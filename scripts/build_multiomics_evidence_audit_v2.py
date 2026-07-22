#!/usr/bin/env python3
"""BTLA multiomics evidence audit v2 — rebuilt directly from raw source tables.

Supersedes audit_v1. Reads the raw assay Excel files directly (not the 08e
long-form derivative), so candidate coverage is no longer limited by the old
`btla_tf_candidate_union.csv` gene gate that 08e iterated over.

Approved specification (see multiomics_source_specification.md):
  * Fold-changes are log2 throughout.
  * Thresholds are the SOURCE-NATIVE values recovered from the study Methods:
      - transcriptomics (DESeq2): padj < 0.05 AND |log2FC| > log2(1.5)
      - proteomics + phosphoproteomics (limma): padj < 0.05 AND |log2FC| > log2(1.3)
    A harmonised sensitivity flag (|log2FC| >= log2(1.5) for all layers) is also
    recorded but is NOT the primary qualifying rule.
  * The non-source "strong" sub-tier from 08e is dropped: each measured item is
    simply qualifying or not-qualifying at the source threshold.
  * Evidence roles:
      - mRNA (transcriptomics): CONTEXTUAL only (regulator's own transcript
        response); never counted as independent corroboration.
      - protein abundance: INDEPENDENT observed molecular evidence; timepoints
        2mn/20mn/4h/24h analysed separately; 4h is the time-matched primary view.
      - phosphosites: INDEPENDENT observed molecular evidence; every site and
        timepoint separate; acute (2/20mn) vs sustained (4/24h) distinguished;
        sites never averaged and phosphorylation direction is NOT read as
        activation without site-specific biology.
      - TF activity: CONTEXTUAL (derived from the transcriptome).
      - coIP: EXCLUDED from corroboration — the available file is an intensity
        matrix with no statistical BTLA contrast; the proper BTLA-IP-vs-IgG
        differential is not in this data tree. Presence is recorded only.
  * Three separated concepts per candidate:
      1. measured_molecular_change (independent assay, any qualifying item)
      2. independent_molecular_association (>=1 qualifying protein/phospho item
         in the BTLA contrast)  -- direction-agnostic
      3. directionally_interpretable_support -- left "unresolved" because model
         rankings are unsigned.
  * Six-state status vocabulary per candidate x assay:
      measured_qualifying | measured_no_qualifying_signal | not_measured |
      absent_from_assay_universe | ambiguous_mapping | failed_qc
  * Missingness is NEVER converted to "no evidence".

Evidence construction only. Manuscript, Table 4, verdict, rankings, CRISPRi and
Paperclip tiers are untouched; rankings are read-only.

Restricted values (effect sizes, adjusted p, phosphosite residues, per-target
direction) go only to a git-ignored local_restricted/ table. Committed tables
carry statuses, counts and provenance hashes.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_DATA_ROOT = "/mnt/R0/Projects/POIAZ/gmorgan/gremln-tcells"
OUT_REL = "results/multiomics/audit_v2"

# Source-native thresholds (recovered from study Methods.docx)
PADJ_MAX = 0.05
RNA_FC_LINEAR = 1.5           # DESeq2: |FC| > 1.5
PROT_FC_LINEAR = 1.3          # limma:  |FC| > 1.3
RNA_LOG2 = math.log2(RNA_FC_LINEAR)     # 0.585
PROT_LOG2 = math.log2(PROT_FC_LINEAR)   # 0.378
HARMONISED_LOG2 = math.log2(1.5)        # labelled sensitivity, all layers

RNA_TIMEPOINTS = ["1h", "4h", "24h"]
MS_TIMEPOINTS = ["2mn", "20mn", "4h", "24h"]
PRIMARY_TP = "4h"
ACUTE_TPS = {"2mn", "20mn"}
SUSTAINED_TPS = {"4h", "24h"}

STATUS = {
    "measured_qualifying", "measured_no_qualifying_signal", "not_measured",
    "absent_from_assay_universe", "ambiguous_mapping", "failed_qc",
    "excluded_no_statistical_contrast",
}


def data_root() -> Path:
    return Path(os.environ.get("BTLA_BENCH_DATA_ROOT", DEFAULT_DATA_ROOT)).expanduser()


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def sha256_file(p: Path) -> str:
    p = Path(p)
    if not p.exists():
        return "missing"
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def utcnow() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _qual(effect, padj, log2_min) -> str:
    """Six-state classification for one measured value at one timepoint."""
    if pd.isna(effect) and pd.isna(padj):
        return "not_measured"
    if pd.isna(effect) or pd.isna(padj):
        return "failed_qc"
    if (padj < PADJ_MAX) and (abs(effect) > log2_min):
        return "measured_qualifying"
    return "measured_no_qualifying_signal"


def _direction(effect) -> str:
    if pd.isna(effect):
        return "unresolved"
    return "increased" if effect > 0 else "decreased" if effect < 0 else "unchanged"


# --------------------------------------------------------------------------- #
def load_union(repo: Path) -> pd.DataFrame:
    u = pd.read_csv(repo / "results/publication_data/top25_union_primary.csv")

    def grp(r):
        if r["in_gremln_top25"] and r["in_genie3_top25"]:
            return "shared"
        return "GREmLN_specific" if r["in_gremln_top25"] else "GENIE3_specific"

    u["candidate_group"] = u.apply(grp, axis=1)
    return u


def build_long(repo: Path, union: pd.DataFrame):
    dr = data_root()
    mm = dr / "data/multiomics"
    cand = list(union["TF"])
    grp = dict(zip(union["TF"], union["candidate_group"]))
    grk = dict(zip(union["TF"], union["gremln_dense_rank"]))
    g3k = dict(zip(union["TF"], union["genie3_dense_rank"]))
    rows = []

    def add(tf, layer, role, contrast, tp, entity, residue, effect, padj,
            log2_min, note):
        status = _qual(effect, padj, log2_min)
        rows.append({
            "TF": tf, "candidate_group": grp[tf],
            "gremln_dense_rank": grk[tf], "genie3_dense_rank": g3k[tf],
            "layer": layer, "evidence_role": role, "contrast": contrast,
            "timepoint": tp, "timepoint_class": ("acute" if tp in ACUTE_TPS
                else "sustained" if tp in SUSTAINED_TPS else "single"),
            "is_primary_timepoint": (tp == PRIMARY_TP),
            "measured_entity": entity, "phosphosite_residue": residue,
            "effect_log2fc": effect, "adjusted_p": padj,
            "qualifying_source": status == "measured_qualifying",
            "qualifying_harmonised": (not pd.isna(effect) and not pd.isna(padj)
                and padj < PADJ_MAX and abs(effect) > HARMONISED_LOG2),
            "status": status, "direction": _direction(effect),
            "note": note,
        })

    # ---- transcriptomics (CONTEXTUAL) ---- BTLA contrast + TCRvsUC context
    rna = pd.read_excel(mm / "Transcriptomics/@Transcriptomics/Tables/DEG_wide_statistics.xlsx")
    rna_u = set(rna["gene"].astype(str))
    rna_by = {str(r["gene"]): r for _, r in rna.iterrows()}
    for tf in cand:
        if tf not in rna_u:
            add(tf, "transcriptomics", "contextual", "BTLA_TCR_vs_TCR", "any",
                tf, "", np.nan, np.nan, RNA_LOG2, "absent_from_assay_universe")
            rows[-1]["status"] = "absent_from_assay_universe"
            continue
        r = rna_by[tf]
        for tp in RNA_TIMEPOINTS:
            add(tf, "transcriptomics", "contextual", "BTLA_TCR_vs_TCR", tp, tf, "",
                r.get(f"BTLAvsTCR_{tp}_log2FoldChange"),
                r.get(f"BTLAvsTCR_{tp}_padj"), RNA_LOG2, "mRNA (regulator's own transcript)")

    # ---- proteomics (INDEPENDENT) ----
    prot = pd.read_excel(mm / "Proteomics/@Proteomics/Tables/WP_DEP_wide_statistics.xlsx")
    prot_u = prot.groupby("Gene").size()
    prot_by = {str(g): sub for g, sub in prot.groupby("Gene")}
    for tf in cand:
        if tf not in prot_by:
            add(tf, "proteomics", "independent", "BTLA_TCR_vs_TCR", "any", tf, "",
                np.nan, np.nan, PROT_LOG2, "absent_from_assay_universe")
            rows[-1]["status"] = "absent_from_assay_universe"
            continue
        sub = prot_by[tf]
        ambiguous = len(sub) > 1
        r = sub.iloc[0]
        for tp in MS_TIMEPOINTS:
            add(tf, "proteomics", "independent", "BTLA_TCR_vs_TCR", tp,
                str(r.get("protein", tf)), "",
                r.get(f"Btla_{tp}_logFC"), r.get(f"Btla_{tp}_adj.P.Val"), PROT_LOG2,
                "multiple protein groups map to gene" if ambiguous else "whole-proteome DEP")
            if ambiguous:
                rows[-1]["status"] = "ambiguous_mapping"

    # ---- phosphoproteomics (INDEPENDENT) ---- per site per timepoint
    phos_path = mm / "Phosphoproteomics/@Phosphoproteomics/@Tables/DEPPs_All_statistics.xlsx"
    xls = pd.ExcelFile(phos_path)
    phos_present = set()
    phos_frames = {}
    for sheet in xls.sheet_names:
        if sheet.startswith("Btla_"):
            phos_frames[sheet.replace("Btla_", "")] = pd.read_excel(phos_path, sheet_name=sheet)
    for df in phos_frames.values():
        phos_present |= set(df["Gene"].astype(str))
    for tf in cand:
        if tf not in phos_present:
            add(tf, "phosphoproteomics", "independent", "BTLA_TCR_vs_TCR", "any", tf, "",
                np.nan, np.nan, PROT_LOG2, "absent_from_assay_universe")
            rows[-1]["status"] = "absent_from_assay_universe"
            continue
        for tp, df in phos_frames.items():
            sites = df[df["Gene"].astype(str) == tf]
            if not len(sites):
                add(tf, "phosphoproteomics", "independent", "BTLA_TCR_vs_TCR", tp, tf, "",
                    np.nan, np.nan, PROT_LOG2, "gene present in assay, no site this timepoint")
                rows[-1]["status"] = "not_measured"
                continue
            for _, s in sites.iterrows():
                site = str(s.get("Gene_STY_canonical", ""))
                residue = site.split("_", 1)[1] if "_" in site else ""
                add(tf, "phosphoproteomics", "independent", "BTLA_TCR_vs_TCR", tp,
                    site, residue, s.get("logFC"), s.get("adj.P.Val"), PROT_LOG2,
                    "phosphosite-level DEP (not averaged)")

    # ---- TF activity (CONTEXTUAL) ----
    tfa = pd.read_excel(mm / "Transcriptomics/@Transcriptomics/Tables/TFactivitiesALL.xlsx")
    tfa_u = set(tfa["source"].astype(str))
    tfa_by = {str(r["source"]): r for _, r in tfa.iterrows()}
    for tf in cand:
        if tf not in tfa_u:
            add(tf, "tf_activity", "contextual", "BTLA_TCR_vs_TCR", "any", tf, "",
                np.nan, np.nan, 0.0, "absent_from_assay_universe")
            rows[-1]["status"] = "absent_from_assay_universe"
            continue
        r = tfa_by[tf]
        for tp in RNA_TIMEPOINTS:
            score = r.get(f"BTLAvsTCR_{tp}_score")
            pval = r.get(f"BTLAvsTCR_{tp}_p_value")
            # activity has no effect-size floor in source; qualify on p only
            status = ("not_measured" if pd.isna(score) and pd.isna(pval)
                      else "failed_qc" if pd.isna(score) or pd.isna(pval)
                      else "measured_qualifying" if pval < PADJ_MAX
                      else "measured_no_qualifying_signal")
            rows.append({
                "TF": tf, "candidate_group": grp[tf],
                "gremln_dense_rank": grk[tf], "genie3_dense_rank": g3k[tf],
                "layer": "tf_activity", "evidence_role": "contextual",
                "contrast": "BTLA_TCR_vs_TCR", "timepoint": tp,
                "timepoint_class": "sustained" if tp in SUSTAINED_TPS else "single",
                "is_primary_timepoint": (tp == PRIMARY_TP),
                "measured_entity": tf, "phosphosite_residue": "",
                "effect_log2fc": score, "adjusted_p": pval,
                "qualifying_source": status == "measured_qualifying",
                "qualifying_harmonised": False,
                "status": status, "direction": _direction(score),
                "note": "decoupleR ULM activity (contextual; p-only, no effect floor)",
            })

    # ---- coIP (EXCLUDED) ---- presence only
    coip = pd.read_excel(mm / "coIP/@coIP/@Tables/@BTLA_coIP_signif_wide.xlsx")
    coip_u = set(coip["Gene"].astype(str))
    for tf in cand:
        present = tf in coip_u
        rows.append({
            "TF": tf, "candidate_group": grp[tf],
            "gremln_dense_rank": grk[tf], "genie3_dense_rank": g3k[tf],
            "layer": "coIP", "evidence_role": "excluded",
            "contrast": "BTLA_coIP_intensity_only", "timepoint": "NA",
            "timepoint_class": "single", "is_primary_timepoint": False,
            "measured_entity": tf, "phosphosite_residue": "",
            "effect_log2fc": np.nan, "adjusted_p": np.nan,
            "qualifying_source": False, "qualifying_harmonised": False,
            "status": "excluded_no_statistical_contrast",
            "direction": "unresolved",
            "note": ("detected in coIP intensity matrix (no IP-vs-IgG stats)"
                     if present else "not detected in coIP intensity matrix"),
        })

    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
def build_summary(long: pd.DataFrame, union: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, u in union.iterrows():
        tf = u["TF"]
        g = long[long["TF"] == tf]

        def layer(l):
            return g[g["layer"] == l]

        prot = layer("proteomics")
        phos = layer("phosphoproteomics")
        rna = layer("transcriptomics")
        tfa = layer("tf_activity")
        coip = layer("coIP")

        prot_qual = prot[prot["qualifying_source"]]
        phos_qual = phos[phos["qualifying_source"]]

        def assay_status(sub):
            if (sub["status"] == "absent_from_assay_universe").all() and len(sub):
                return "absent_from_assay_universe"
            if (sub["qualifying_source"]).any():
                return "measured_qualifying"
            if (sub["status"] == "measured_no_qualifying_signal").any():
                return "measured_no_qualifying_signal"
            if (sub["status"] == "ambiguous_mapping").any():
                return "ambiguous_mapping"
            if (sub["status"] == "not_measured").any():
                return "not_measured"
            return "not_measured"

        indep_assoc = bool(len(prot_qual) or len(phos_qual))
        measured_change = indep_assoc  # any qualifying independent item

        # mRNA context: which timepoints qualify (regulator's own transcript)
        rna_tps = sorted(rna[rna["qualifying_source"]]["timepoint"].tolist())
        tfa_tps = sorted(tfa[tfa["qualifying_source"]]["timepoint"].tolist())

        rows.append({
            "TF": tf, "candidate_group": u["candidate_group"],
            "gremln_dense_rank": u["gremln_dense_rank"],
            "genie3_dense_rank": u["genie3_dense_rank"],
            "in_gremln_top25": bool(u["in_gremln_top25"]),
            "in_genie3_top25": bool(u["in_genie3_top25"]),
            # concept 1 & 2
            "measured_molecular_change": measured_change,
            "independent_molecular_association": indep_assoc,
            # concept 3
            "directionally_interpretable_support": "unresolved",
            # independent detail
            "protein_status": assay_status(prot),
            "protein_qualifying_timepoints": ";".join(sorted(prot_qual["timepoint"].tolist())),
            "protein_primary_4h_qualifying": bool(
                len(prot_qual[prot_qual["timepoint"] == PRIMARY_TP])),
            "protein_acute_qualifying": bool(
                len(prot_qual[prot_qual["timepoint"].isin(ACUTE_TPS)])),
            "phospho_status": assay_status(phos),
            "phospho_qualifying_sites": int(len(phos_qual)),
            "phospho_qualifying_site_timepoints": ";".join(
                f"{r.measured_entity}@{r.timepoint}" for r in phos_qual.itertuples()),
            "phospho_acute_qualifying_sites": int(
                len(phos_qual[phos_qual["timepoint"].isin(ACUTE_TPS)])),
            "phospho_sustained_qualifying_sites": int(
                len(phos_qual[phos_qual["timepoint"].isin(SUSTAINED_TPS)])),
            # contextual detail
            "mrna_status": assay_status(rna),
            "mrna_qualifying_timepoints": ";".join(rna_tps),
            "tf_activity_status": assay_status(tfa),
            "tf_activity_qualifying_timepoints": ";".join(tfa_tps),
            "coip_presence": ("detected_intensity_only"
                if "detected" in coip.iloc[0]["note"] and "not detected" not in coip.iloc[0]["note"]
                else "not_detected"),
            "overall_independent_status": (
                "independent_association_direction_unresolved" if indep_assoc
                else "no_independent_qualifying_signal"),
        })
    return pd.DataFrame(rows)


def model_summaries(summ: pd.DataFrame) -> pd.DataFrame:
    views = {
        "GREmLN_top25": summ[summ["in_gremln_top25"]],
        "GENIE3_top25": summ[summ["in_genie3_top25"]],
        "shared": summ[summ["candidate_group"] == "shared"],
        "GREmLN_specific": summ[summ["candidate_group"] == "GREmLN_specific"],
        "GENIE3_specific": summ[summ["candidate_group"] == "GENIE3_specific"],
        "union_43": summ,
    }
    rows = []
    for name, sub in views.items():
        rows.append({
            "view": name, "n_candidates": int(len(sub)),
            "n_independent_association": int(sub["independent_molecular_association"].sum()),
            "candidates_independent": ";".join(
                sub[sub["independent_molecular_association"]]["TF"].tolist()),
            "n_protein_qualifying": int((sub["protein_status"] == "measured_qualifying").sum()),
            "n_phospho_qualifying": int((sub["phospho_qualifying_sites"] > 0).sum()),
            "n_mrna_context_change": int((sub["mrna_qualifying_timepoints"] != "").sum()),
            "n_tf_activity_context": int((sub["tf_activity_qualifying_timepoints"] != "").sum()),
            "n_protein_absent_universe": int((sub["protein_status"] == "absent_from_assay_universe").sum()),
            "n_phospho_absent_universe": int((sub["phospho_status"] == "absent_from_assay_universe").sum()),
        })
    return pd.DataFrame(rows)


def coverage_table(long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for layer in ["transcriptomics", "proteomics", "phosphoproteomics",
                  "tf_activity", "coIP"]:
        g = long[long["layer"] == layer]
        by_tf = g.groupby("TF")["status"]
        present = by_tf.apply(lambda s: not (s == "absent_from_assay_universe").all())
        rows.append({
            "layer": layer,
            "evidence_role": g["evidence_role"].iloc[0],
            "candidates_in_universe": int(present.sum()),
            "candidates_absent": int((~present).sum()),
            "candidates_qualifying": int(
                g[g["qualifying_source"]]["TF"].nunique()),
        })
    return pd.DataFrame(rows)


def run_checks(union, long, summ, models, coverage) -> pd.DataFrame:
    checks = []

    def chk(name, cond, detail=""):
        checks.append({"check": name, "passed": bool(cond), "detail": str(detail)})

    chk("union == 43", len(union) == 43, len(union))
    chk("summary rows == 43", len(summ) == 43, len(summ))
    chk("GREmLN view == 25", int(summ["in_gremln_top25"].sum()) == 25,
        int(summ["in_gremln_top25"].sum()))
    chk("GENIE3 view == 25", int(summ["in_genie3_top25"].sum()) == 25,
        int(summ["in_genie3_top25"].sum()))
    chk("shared == 7", int((summ["candidate_group"] == "shared").sum()) == 7,
        int((summ["candidate_group"] == "shared").sum()))
    chk("GREmLN_specific + shared == 25",
        int((summ["candidate_group"].isin(["GREmLN_specific", "shared"])).sum()) == 25)
    chk("GENIE3_specific + shared == 25",
        int((summ["candidate_group"].isin(["GENIE3_specific", "shared"])).sum()) == 25)
    # mRNA never independent
    chk("mRNA rows are all contextual role",
        (long[long["layer"] == "transcriptomics"]["evidence_role"] == "contextual").all())
    chk("tf_activity rows are all contextual role",
        (long[long["layer"] == "tf_activity"]["evidence_role"] == "contextual").all())
    chk("coIP excluded from corroboration",
        (long[long["layer"] == "coIP"]["evidence_role"] == "excluded").all())
    # independent association derives only from protein/phospho qualifying
    for _, s in summ[summ["independent_molecular_association"]].iterrows():
        sub = long[(long["TF"] == s["TF"]) &
                   (long["layer"].isin(["proteomics", "phosphoproteomics"])) &
                   (long["qualifying_source"])]
        if not len(sub):
            chk(f"independent association backed by qualifying item ({s['TF']})", False)
            break
    else:
        chk("every independent association backed by a qualifying protein/phospho item", True)
    # not_measured never qualifying
    chk("not_measured never qualifying",
        not long[(long["status"] == "not_measured") & (long["qualifying_source"])].shape[0])
    chk("absent_from_assay_universe never qualifying",
        not long[(long["status"] == "absent_from_assay_universe") & (long["qualifying_source"])].shape[0])
    # thresholds applied match source-native
    chk("RNA threshold = log2(1.5)", abs(RNA_LOG2 - math.log2(1.5)) < 1e-9, RNA_LOG2)
    chk("protein/phospho threshold = log2(1.3)", abs(PROT_LOG2 - math.log2(1.3)) < 1e-9, PROT_LOG2)
    # coverage reconciles with model summaries
    chk("independent count reconciles (union)",
        int(summ["independent_molecular_association"].sum()) ==
        int(models[models["view"] == "union_43"]["n_independent_association"].iloc[0]))
    return pd.DataFrame(checks)


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", default=None)
    args = ap.parse_args(argv)
    repo = Path(args.repo_root).resolve() if args.repo_root else repo_root()
    dr = data_root()
    out = repo / OUT_REL
    (out / "local_restricted").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    union = load_union(repo)
    long = build_long(repo, union)
    summ = build_summary(long, union)
    models = model_summaries(summ)
    coverage = coverage_table(long)
    checks = run_checks(union, long, summ, models, coverage)

    # restricted full
    long.to_csv(out / "local_restricted" / "multiomics_evidence_long_v2.FULL.csv", index=False)
    # redacted committed
    red = long.drop(columns=["effect_log2fc", "adjusted_p", "phosphosite_residue",
                             "direction"]).copy()
    red["effect_log2fc"] = "REDACTED_publication_gate"
    # keep measured_entity for coIP/protein/rna/tfa but strip phospho residue-bearing site id
    red.loc[red["layer"] == "phosphoproteomics", "measured_entity"] = "REDACTED_site"
    red.to_csv(out / "multiomics_evidence_long_v2.csv", index=False)

    # summary: strip phospho site-timepoint detail (residue-bearing) from committed copy
    summ_red = summ.copy()
    summ_red["phospho_qualifying_site_timepoints"] = summ_red[
        "phospho_qualifying_site_timepoints"].apply(
        lambda s: "REDACTED_sites" if s else "")
    summ_red.to_csv(out / "multiomics_candidate_summary_v2.csv", index=False)
    summ.to_csv(out / "local_restricted" / "multiomics_candidate_summary_v2.FULL.csv", index=False)

    models.to_csv(out / "multiomics_model_summaries_v2.csv", index=False)
    coverage.to_csv(out / "multiomics_coverage_v2.csv", index=False)
    checks.to_csv(out / "multiomics_audit_checks_v2.csv", index=False)

    manifest = {
        "generated_utc": utcnow(),
        "supersedes": "results/multiomics/audit_v1",
        "thresholds_source_native": {
            "padj_max": PADJ_MAX, "rna_abs_fc_linear": RNA_FC_LINEAR,
            "prot_phospho_abs_fc_linear": PROT_FC_LINEAR,
            "rna_abs_log2fc": RNA_LOG2, "prot_phospho_abs_log2fc": PROT_LOG2,
            "harmonised_sensitivity_abs_log2fc": HARMONISED_LOG2,
            "recovered_from": "data/multiomics/Documents/Documents/Methods.docx",
        },
        "n_candidates": int(len(union)),
        "n_independent_association": int(summ["independent_molecular_association"].sum()),
        "candidates_independent": summ[summ["independent_molecular_association"]]["TF"].tolist(),
        "coverage": coverage.to_dict(orient="records"),
        "checks_passed": int(checks["passed"].sum()),
        "checks_total": int(len(checks)),
        "raw_source_sha256": {
            "transcriptomics": sha256_file(dr / "data/multiomics/Transcriptomics/@Transcriptomics/Tables/DEG_wide_statistics.xlsx"),
            "proteomics": sha256_file(dr / "data/multiomics/Proteomics/@Proteomics/Tables/WP_DEP_wide_statistics.xlsx"),
            "phosphoproteomics": sha256_file(dr / "data/multiomics/Phosphoproteomics/@Phosphoproteomics/@Tables/DEPPs_All_statistics.xlsx"),
            "coIP": sha256_file(dr / "data/multiomics/coIP/@coIP/@Tables/@BTLA_coIP_signif_wide.xlsx"),
            "tf_activity": sha256_file(dr / "data/multiomics/Transcriptomics/@Transcriptomics/Tables/TFactivitiesALL.xlsx"),
            "methods_doc": sha256_file(dr / "data/multiomics/Documents/Documents/Methods.docx"),
        },
    }
    (out / "run_manifest_v2.json").write_text(json.dumps(manifest, indent=2))

    print(f"[audit_v2] checks {int(checks['passed'].sum())}/{len(checks)} passed")
    for _, c in checks[~checks["passed"]].iterrows():
        print("  FAIL:", c["check"], c["detail"])
    print("[audit_v2] independent association:",
          manifest["candidates_independent"])
    print(coverage.to_string(index=False))
    return 0 if bool(checks["passed"].all()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
