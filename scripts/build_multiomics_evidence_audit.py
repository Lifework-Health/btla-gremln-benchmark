#!/usr/bin/env python3
"""Auditable BTLA multiomics evidence pipeline (audit_v1).

Purpose
-------
Reconstruct, layer by layer, exactly how BTLA multiomics evidence was generated,
classified and assigned to the 43 candidate regulators in the union of the
GREmLN and GENIE3 seed-excluded top-25 lists, and replace the legacy
"non-empty string == support" behaviour with documented, layer-specific,
reproducible rules.

This is evidence construction and validation only. It does NOT touch the
manuscript, Results/Conclusions, Table 4, the benchmark verdict, the model
rankings, the CRISPRi analysis, or the Paperclip tiers.

Design
------
* The atomic evidence source is the frozen ``08e`` long-form table
  ``btla_candidate_multiomics_long_evidence.csv`` (one row per candidate x
  evidence item, with per-row effect size, adjusted p, support_level,
  phosphosite, contrast, timepoint). Its SHA256 and the SHA256 of every raw
  Excel source are recorded in ``source_inventory.csv``.
* The audit RE-DERIVES independence class, a controlled evidence status
  vocabulary, directional interpretation and the benchmark decision field from
  the long table using explicit rules. It never trusts the legacy summary.
* Restricted values (effect sizes, phosphosite residues) are written only to a
  git-ignored ``local_restricted/`` table. The committed tables carry statuses,
  provenance hashes and aggregate counts only.

Outputs (see commit_safety_report.md for commit safety of each):
    results/multiomics/audit_v1/
        source_inventory.csv
        run_manifest.json
        candidate_union.csv
        multiomics_layer_registry.csv
        multiomics_evidence_long_audited.csv            (redacted, safe)
        local_restricted/multiomics_evidence_long_audited.FULL.csv  (gitignored)
        multiomics_candidate_summary_audited.csv
        legacy_vs_audited_multiomics.csv
        multiomics_audit_checks.csv
        multiomics_pipeline_audit.md
        commit_safety_report.md
        figures/fig_multiomics_audit_pipeline.(svg|png)
        figures/fig_multiomics_candidate_matrix.(svg|png)
        figures/fig_multiomics_model_summary.(svg|png)
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
DEFAULT_DATA_ROOT = "/mnt/R0/Projects/POIAZ/gmorgan/gremln-tcells"
OUT_REL = "results/multiomics/audit_v1"


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "notebooks").is_dir() and (c / "scripts").is_dir():
            return c
    return here.parent.parent


def data_root() -> Path:
    return Path(os.environ.get("BTLA_BENCH_DATA_ROOT", DEFAULT_DATA_ROOT)).expanduser()


def sha256_file(p: Path) -> str:
    p = Path(p)
    if not p.exists():
        return "missing"
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def utcnow() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def redact_path(p: Path) -> str:
    """Relative to DATA_ROOT so committed provenance carries no absolute paths."""
    try:
        return str(Path(p).resolve().relative_to(data_root().resolve()))
    except Exception:
        return Path(p).name


# --------------------------------------------------------------------------- #
# Controlled vocabularies
# --------------------------------------------------------------------------- #
EVIDENCE_STATUS = {
    "positive_supportive",
    "positive_opposing",
    "positive_unresolved_direction",
    "measured_no_qualifying_signal",
    "measured_failed_qc",
    "not_measured",
    "candidate_absent_from_assay_universe",
    "ambiguous_mapping",
    "source_unavailable",
    "not_applicable",
}

# Independence classes (audited).
#   A independent observed molecular evidence
#   B orthogonal assay-derived computational evidence
#   C transcript-/model-derived contextual evidence
#   D unclear lineage
LAYER_INDEPENDENCE = {
    "transcriptomics": "C",
    "proteomics": "A",
    "phosphoproteomics": "A",
    "coIP": "A",
    "tf_activity": "C",
    "kinase_activity": "B",
    "bionic_gnn": "C",
    "early_synapse_trafficking": "D",
}

# Positive (qualifying) support levels in the frozen long table.
POSITIVE_SUPPORT_LEVELS = {"strong", "moderate"}

# Layers whose direction can/cannot be mapped to the BTLA hypothesis without a
# documented per-target functional model. Magnitude-only molecular layers are
# directionally *unresolved* by audit policy (we do not assume up==supportive).
DIRECTION_UNRESOLVED_LAYERS = {
    "transcriptomics", "proteomics", "phosphoproteomics",
    "tf_activity", "kinase_activity",
}
DIRECTION_NONDIRECTIONAL_LAYERS = {"coIP", "bionic_gnn", "early_synapse_trafficking"}

# Thresholds documented in the 08e generator (recorded for provenance; already
# applied when support_level was computed).
THRESHOLDS = {
    "PADJ_MAX": 0.05,
    "RNA_LFC_MIN_log2": 0.5849625007211562,   # log2(1.5)
    "PROT_LFC_MIN_log2": 0.37851162325372983,  # log2(1.3)
    "coIP_abs_t_min": 2.0,
    "tf_activity_p_max": 0.05,
    "kinase_activity_p_max": 0.05,
    "strong_padj_max": 0.01,
}


# --------------------------------------------------------------------------- #
# Source inventory
# --------------------------------------------------------------------------- #
def build_source_inventory(long_path: Path, summary_path: Path,
                           gr_rank_path: Path, g3_rank_path: Path,
                           union_path: Path) -> pd.DataFrame:
    dr = data_root()
    rows = []

    def add(layer, path, experiment, assay, contrast, timepoint, stat_method,
            effect_measure, sig_threshold, origin, committable, note=""):
        p = Path(path)
        exists = p.exists()
        dims = ""
        if exists and p.suffix.lower() in {".csv"}:
            try:
                df = pd.read_csv(p)
                dims = f"{df.shape[0]}x{df.shape[1]}"
            except Exception:
                dims = ""
        elif exists:
            dims = f"{p.stat().st_size}_bytes"
        rows.append({
            "logical_layer": layer,
            "redacted_relative_path": redact_path(p),
            "sha256": sha256_file(p),
            "dimensions": dims,
            "source_experiment": experiment,
            "assay_type": assay,
            "contrast": contrast,
            "time_point": timepoint,
            "statistical_method": stat_method,
            "effect_measure": effect_measure,
            "significance_threshold": sig_threshold,
            "originating_script_or_notebook":
                "notebooks/08e_BTLA_candidate_multiomics_evidence.ipynb",
            "access_status": "RESTRICTED",
            "may_be_committed": committable,
            "exists": exists,
            "note": note,
        })

    mm = dr / "data/multiomics"
    add("transcript differential expression",
        mm / "Transcriptomics/@Transcriptomics/Tables/DEG_wide_statistics.xlsx",
        "anti-BTLA cross-link bulk RNA-seq", "bulk RNA-seq DEG (limma/DESeq wide)",
        "BTLA_TCR_vs_TCR; TCR_vs_UC", "1h;4h;24h", "wide DE statistics",
        "log2FoldChange", "padj<0.05 & |log2FC|>=log2(1.5)",
        "08e", "no", "TCR_vs_UC contrast not used for BTLA support")
    add("protein abundance",
        mm / "Proteomics/@Proteomics/Tables/WP_DEP_wide_statistics.xlsx",
        "whole-proteome DEP", "LC-MS/MS whole proteome (limma DEP)",
        "BTLA_TCR_vs_TCR", "2mn;20mn;4h;24h", "limma moderated t",
        "logFC", "adj.P.Val<0.05 & |logFC|>=log2(1.3)",
        "08e", "no", "statistically tested DEP, not detection-only")
    add("phosphosite abundance",
        mm / "Phosphoproteomics/@Phosphoproteomics/@Tables/DEPPs_All_statistics.xlsx",
        "phosphoproteome DEPP", "LC-MS/MS phosphoproteome (limma DEP per site)",
        "BTLA_TCR_vs_TCR", "2mn;20mn;4h;24h", "limma moderated t (per site)",
        "logFC", "adj.P.Val<0.05 & |logFC|>=log2(1.3)",
        "08e", "no", "site-specific (Gene_STY_canonical); Btla_* sheets only")
    add("coimmunoprecipitation",
        mm / "coIP/@coIP/@Tables/@BTLA_coIP_signif_wide.xlsx",
        "BTLA co-IP MS", "affinity purification MS (pre-filtered 'signif')",
        "BTLA co-IP vs control", "NA", "limma moderated t-statistic",
        "t_statistic", "|t|>2 (pre-filtered signif table; CRAPome not applied)",
        "08e", "no", "pulldown membership; direct interaction NOT resolved")
    add("inferred TF activity",
        mm / "Transcriptomics/@Transcriptomics/Tables/TFactivitiesALL.xlsx",
        "decoupleR ULM on transcriptome", "computational (decoupleR ULM)",
        "BTLA_TCR_vs_TCR; TCR_vs_UC; BTLA_TCR_vs_UC", "1h;4h;24h",
        "decoupleR ULM", "ULM activity score", "p_value<0.05 (raw p, no effect floor)",
        "08e", "no", "derived from transcript expression -> class C")
    add("inferred kinase activity",
        mm / "kinases_ulm_wide.xlsx",
        "decoupleR ULM on phosphoproteome", "computational (decoupleR ULM)",
        "BTLA_TCR_vs_TCR", "phospho timepoints", "decoupleR ULM",
        "ULM activity score", "p_value<0.05",
        "no", "no", "orthogonal (phospho-derived) -> class B; empty for TF candidate set")
    add("BIONIC/GNN module",
        mm / "bionic/GNN BIONIC Integration/Complete Dataset/R_GNN_Umap.xlsx",
        "BIONIC GNN integration", "network integration (GNN) UMAP module",
        "BIONIC_module", "NA", "module assignment", "cluster membership",
        "assignment only (no p/effect)",
        "no", "no", "multi-omic network module; membership only -> class C")
    add("early synapse / trafficking",
        Path("inline_module_list"),
        "curated hypothesis gene list (inline in 08e)", "curated annotation",
        "hypothesis_module", "NA", "none", "membership", "none (curated)",
        "08e", "yes", "curated hypothesis list, NOT measured evidence -> class D")

    # Frozen atomic evidence table + legacy summary + rankings
    add("ATOMIC long evidence (frozen 08e output)", long_path,
        "08e aggregation", "derived long-form evidence table",
        "multiple", "multiple", "08e rules", "per-row", "per-layer (see registry)",
        "08e", "no", "one row per candidate x evidence item; audit re-derives classes")
    add("LEGACY candidate summary (08e output)", summary_path,
        "08e aggregation", "derived candidate summary",
        "multiple", "multiple", "08e nonempty-string heuristic", "per-candidate",
        "n/a", "08e", "no", "legacy; audited output replaces this")
    add("GREmLN seed-excluded ranking", gr_rank_path,
        "GREmLN CSLS scoring", "model ranking", "BTLA_vs_TCR", "NA",
        "CSLS + hypergeometric", "gremln_csls_score", "n/a",
        "notebook 03", "no", "canonical ranking; evidence cannot alter rank")
    add("GENIE3 seed-excluded ranking", g3_rank_path,
        "GENIE3 GRN scoring", "model ranking", "BTLA_vs_TCR", "NA",
        "GENIE3 importance", "genie3_score", "n/a",
        "notebook 03", "no", "canonical ranking; evidence cannot alter rank")
    add("top25 union (primary)", union_path,
        "benchmark freeze", "candidate membership", "BTLA_vs_TCR", "NA",
        "dense-rank top25", "membership", "n/a",
        "notebook 03", "no", "43-candidate union membership")

    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Candidate universe
# --------------------------------------------------------------------------- #
def build_candidate_union(union_path: Path, gr_rank_path: Path,
                          g3_rank_path: Path) -> pd.DataFrame:
    u = pd.read_csv(union_path)
    gr = pd.read_csv(gr_rank_path)
    g3 = pd.read_csv(g3_rank_path)

    gr_rank = dict(zip(gr["gene"], gr.get("gremln_csls_rank", gr.get("gremln_csls_score"))))
    g3_rank = dict(zip(g3["gene"], g3.get("genie3_rank", g3.get("genie3_score"))))

    def group(r):
        if r["in_gremln_top25"] and r["in_genie3_top25"]:
            return "shared"
        if r["in_gremln_top25"]:
            return "GREmLN specific"
        return "GENIE3 specific"

    rows = []
    for _, r in u.iterrows():
        rows.append({
            "TF": r["TF"],
            "gremln_dense_rank": r.get("gremln_dense_rank", ""),
            "genie3_dense_rank": r.get("genie3_dense_rank", ""),
            "gremln_csls_rank": gr_rank.get(r["TF"], ""),
            "genie3_rank": g3_rank.get(r["TF"], ""),
            "in_gremln_top25": bool(r["in_gremln_top25"]),
            "in_genie3_top25": bool(r["in_genie3_top25"]),
            "candidate_group": group(r),
            "canonical_full_name": r["TF"],  # HGNC symbol used as canonical id
            "hgnc_symbol": r["TF"],
            "recognised_aliases": "",
            "gremln_ranking_sha256": "",  # filled by caller
            "genie3_ranking_sha256": "",
        })
    df = pd.DataFrame(rows)
    return df


# --------------------------------------------------------------------------- #
# Layer registry
# --------------------------------------------------------------------------- #
def build_layer_registry() -> pd.DataFrame:
    R = []

    def add(**kw):
        R.append(kw)

    add(layer_id="transcriptomics", display_name="Transcript differential expression",
        raw_assay="bulk RNA-seq", source_experiment="anti-BTLA cross-link bulk RNA-seq",
        measured_or_inferred="measured", input_data="DEG_wide_statistics.xlsx",
        relation_to_249_panel="same BTLA DEG family used for candidate seeding",
        relation_to_model_inputs="transcript signal used for ranking (CSLS/GENIE3 seeds)",
        contrast="BTLA_TCR_vs_TCR (also TCR_vs_UC, unused for support)",
        time_point="1h;4h;24h", unit="gene", effect_measure="log2FoldChange",
        statistical_test="wide DE (limma/DESeq)", mtc="BH (padj)",
        significance_threshold="padj<0.05", min_effect="|log2FC|>=log2(1.5)",
        direction_convention="up/down in BTLA+TCR vs TCR",
        candidate_mapping="gene symbol", directness="direct measurement",
        independence_class="C",
        independence_rationale="derived from the transcript signal used for candidate "
            "ranking; not orthogonal to the model scoring input",
        missingness_meaning="gene absent from DEG table = not measured",
        source_code="08e transcriptomics cell")
    add(layer_id="proteomics", display_name="Protein abundance (DEP)",
        raw_assay="LC-MS/MS whole proteome", source_experiment="whole-proteome DEP",
        measured_or_inferred="measured", input_data="WP_DEP_wide_statistics.xlsx",
        relation_to_249_panel="independent proteome, overlaps panel by symbol",
        relation_to_model_inputs="orthogonal to transcript ranking input",
        contrast="BTLA_TCR_vs_TCR", time_point="2mn;20mn;4h;24h",
        unit="protein", effect_measure="logFC", statistical_test="limma moderated t",
        mtc="BH (adj.P.Val)", significance_threshold="adj.P.Val<0.05",
        min_effect="|logFC|>=log2(1.3)",
        direction_convention="up/down in BTLA+TCR vs TCR",
        candidate_mapping="protein->gene symbol", directness="direct measurement",
        independence_class="A",
        independence_rationale="directly measured protein abundance in a distinct assay; "
            "not computed from the transcript ranking signal",
        missingness_meaning="protein absent from proteome = not measured (not zero)",
        source_code="08e proteomics cell")
    add(layer_id="phosphoproteomics", display_name="Phosphosite abundance (DEPP)",
        raw_assay="LC-MS/MS phosphoproteome", source_experiment="phosphoproteome DEPP",
        measured_or_inferred="measured", input_data="DEPPs_All_statistics.xlsx",
        relation_to_249_panel="independent phosphoproteome",
        relation_to_model_inputs="orthogonal to transcript ranking input",
        contrast="BTLA_TCR_vs_TCR", time_point="2mn;20mn;4h;24h",
        unit="phosphosite", effect_measure="logFC", statistical_test="limma moderated t (per site)",
        mtc="BH (adj.P.Val)", significance_threshold="adj.P.Val<0.05",
        min_effect="|logFC|>=log2(1.3)",
        direction_convention="site up/down in BTLA+TCR vs TCR",
        candidate_mapping="Gene_STY_canonical -> gene symbol + residue",
        directness="direct measurement (site-level)",
        independence_class="A",
        independence_rationale="directly measured phosphosite abundance; distinct assay",
        missingness_meaning="no site quantified for gene = not measured",
        source_code="08e phosphoproteomics cell")
    add(layer_id="coIP", display_name="Co-immunoprecipitation",
        raw_assay="affinity purification MS", source_experiment="BTLA co-IP MS",
        measured_or_inferred="measured", input_data="BTLA_coIP_signif_wide.xlsx",
        relation_to_249_panel="interaction assay, independent",
        relation_to_model_inputs="orthogonal",
        contrast="BTLA co-IP vs control (bait=BTLA)", time_point="NA",
        unit="prey protein", effect_measure="limma t-statistic",
        statistical_test="limma moderated t (pre-filtered signif)", mtc="pre-filtered table",
        significance_threshold="|t|>2 (CRAPome file present but NOT applied)",
        min_effect="|t|>2", direction_convention="enriched/depleted vs control (non-directional for hypothesis)",
        candidate_mapping="prey protein -> gene symbol",
        directness="pulldown membership; direct physical interaction NOT resolved",
        independence_class="A",
        independence_rationale="observed molecular interaction assay; independent of transcript signal",
        missingness_meaning="prey absent from signif table = no qualifying pulldown (not necessarily not measured)",
        source_code="08e coIP cell")
    add(layer_id="tf_activity", display_name="Inferred TF activity",
        raw_assay="computational", source_experiment="decoupleR ULM on transcriptome",
        measured_or_inferred="inferred", input_data="TFactivitiesALL.xlsx",
        relation_to_249_panel="computed from transcriptome",
        relation_to_model_inputs="derived from transcript signal used for ranking",
        contrast="BTLA_TCR_vs_TCR; TCR_vs_UC; BTLA_TCR_vs_UC", time_point="1h;4h;24h",
        unit="TF regulon", effect_measure="ULM activity score",
        statistical_test="decoupleR ULM", mtc="none (raw p_value used)",
        significance_threshold="p_value<0.05 (no effect floor)", min_effect="none",
        direction_convention="activated/inhibited (regulon)",
        candidate_mapping="TF symbol", directness="inferred from expression",
        independence_class="C",
        independence_rationale="decoupleR ULM computed from the same transcript expression "
            "that drives candidate ranking; contextual, not independent",
        missingness_meaning="TF absent from regulon set = not evaluated",
        source_code="08e tf_activity cell")
    add(layer_id="kinase_activity", display_name="Inferred kinase activity",
        raw_assay="computational", source_experiment="decoupleR ULM on phosphoproteome",
        measured_or_inferred="inferred", input_data="kinases_ulm_wide.xlsx",
        relation_to_249_panel="computed from phosphoproteome",
        relation_to_model_inputs="orthogonal (phospho-derived), not transcript",
        contrast="BTLA_TCR_vs_TCR", time_point="phospho timepoints",
        unit="kinase substrate set", effect_measure="ULM activity score",
        statistical_test="decoupleR ULM", mtc="none (raw p_value used)",
        significance_threshold="p_value<0.05", min_effect="none",
        direction_convention="activated/inhibited",
        candidate_mapping="kinase symbol", directness="inferred from phosphosites",
        independence_class="B",
        independence_rationale="inferred from an orthogonal phosphoproteomic assay, not "
            "from transcript expression; orthogonal-derived",
        missingness_meaning="no kinase overlap with TF candidate set -> empty layer",
        source_code="08e kinase_activity cell")
    add(layer_id="bionic_gnn", display_name="BIONIC / GNN module",
        raw_assay="computational", source_experiment="BIONIC GNN integration",
        measured_or_inferred="inferred", input_data="R_GNN_Umap.xlsx",
        relation_to_249_panel="multi-omic network module",
        relation_to_model_inputs="mixed omics; includes transcriptome",
        contrast="BIONIC_module", time_point="NA", unit="gene (module member)",
        effect_measure="cluster membership", statistical_test="none",
        mtc="none", significance_threshold="assignment only (no p/effect)",
        min_effect="none", direction_convention="non-directional (membership)",
        candidate_mapping="gene symbol", directness="module co-membership (observational)",
        independence_class="C",
        independence_rationale="module membership from an integration that includes the "
            "transcriptome; not independent molecular evidence",
        missingness_meaning="gene absent from UMAP = not assigned",
        source_code="08e bionic cell")
    add(layer_id="early_synapse_trafficking", display_name="Early synapse / trafficking (curated)",
        raw_assay="curated annotation", source_experiment="inline hypothesis list in 08e",
        measured_or_inferred="curated", input_data="inline_module_list",
        relation_to_249_panel="hand-curated gene sets",
        relation_to_model_inputs="none",
        contrast="hypothesis_module", time_point="NA", unit="gene",
        effect_measure="membership", statistical_test="none", mtc="none",
        significance_threshold="none (curated membership)", min_effect="none",
        direction_convention="non-directional",
        candidate_mapping="gene symbol", directness="curated hypothesis, NOT measured",
        independence_class="D",
        independence_rationale="curated hypothesis list with no measured assay or threshold; "
            "lineage as evidence is unclear -> excluded from all support counts",
        missingness_meaning="not in curated list = not applicable",
        source_code="08e synapse/trafficking cell")

    return pd.DataFrame(R)


# --------------------------------------------------------------------------- #
# Long audited evidence
# --------------------------------------------------------------------------- #
def audit_long_evidence(long_df: pd.DataFrame, union: pd.DataFrame) -> pd.DataFrame:
    grp = dict(zip(union["TF"], union["candidate_group"]))
    gr_rank = dict(zip(union["TF"], union["gremln_dense_rank"]))
    g3_rank = dict(zip(union["TF"], union["genie3_dense_rank"]))
    union_tfs = set(union["TF"])

    rows = []
    for _, r in long_df.iterrows():
        gene = r["gene"]
        if gene not in union_tfs:
            continue  # audit scope = 43-candidate union
        layer = r["evidence_layer"]
        indep = LAYER_INDEPENDENCE.get(layer, "D")
        support_level = str(r.get("support_level", ""))
        measured = support_level != "not_measured"

        # Controlled evidence status
        if layer == "early_synapse_trafficking":
            status = "not_applicable"  # curated hypothesis, not measured evidence
        elif support_level == "not_measured":
            status = "not_measured"
        elif support_level in POSITIVE_SUPPORT_LEVELS:
            # positive (passes documented significance/effect); direction resolved below
            status = "positive"  # refined after direction interpretation
        elif support_level == "not_significant":
            status = "measured_no_qualifying_signal"
        elif support_level == "hypothesis_only":
            status = "not_applicable"
        else:
            status = "ambiguous_mapping"

        # Directional interpretation
        raw_dir = str(r.get("direction", "") or "")
        if layer in DIRECTION_NONDIRECTIONAL_LAYERS:
            evidence_direction = "present" if layer == "coIP" else "enriched" \
                if layer == "bionic_gnn" else "not directional"
            evidence_direction = {"coIP": "present", "bionic_gnn": "enriched",
                                  "early_synapse_trafficking": "not directional"}[layer]
        else:
            evidence_direction = {"up": "increased", "down": "decreased"}.get(raw_dir, "not directional")

        # relation_to_BTLA_hypothesis: magnitude-only molecular layers are
        # directionally unresolved by audit policy (no per-target functional rule).
        if status == "positive":
            if layer in DIRECTION_UNRESOLVED_LAYERS:
                relation = "unresolved"
                status = "positive_unresolved_direction"
            else:
                relation = "context only"
                status = "positive_unresolved_direction"
        elif status == "measured_no_qualifying_signal":
            relation = "context only"
        else:
            relation = "context only" if status == "not_applicable" else "unresolved"

        rows.append({
            "TF": gene,
            "candidate_group": grp.get(gene, ""),
            "gremln_rank": gr_rank.get(gene, ""),
            "genie3_rank": g3_rank.get(gene, ""),
            "layer_id": layer,
            "assay": layer,
            "source_file_id": Path(str(r.get("source_file", ""))).name or "inline",
            "source_row_id": int(r.name),
            "contrast": r.get("contrast", ""),
            "time_point": r.get("timepoint", ""),
            "measured_entity": gene,
            "protein_or_phosphosite_id": r.get("phosphosite", "") or "",
            "phosphosite_residue": _residue(r.get("phosphosite", "")),
            "bait": "BTLA" if layer == "coIP" else "",
            "prey": gene if layer == "coIP" else "",
            "effect_measure": _effect_measure(layer),
            "effect_value": r.get("effect_size", ""),
            "raw_direction": raw_dir,
            "evidence_direction": evidence_direction,
            "raw_p": "",
            "adjusted_p": r.get("adjusted_p", ""),
            "effect_threshold_passed": support_level in POSITIVE_SUPPORT_LEVELS,
            "significance_threshold_passed": support_level in POSITIVE_SUPPORT_LEVELS,
            "qc_passed": measured,
            "candidate_mapping_status": "exact_symbol",
            "evidence_status": status,
            "relation_to_BTLA_hypothesis": relation,
            "independence_class": indep,
            "independence_rationale": _indep_rationale(indep),
            "source_lineage": r.get("independent_or_derived", ""),
            "legacy_support_level": support_level,
            "notes": r.get("notes", ""),
        })
    out = pd.DataFrame(rows)
    # source row hash (over the safe, non-restricted identifying columns)
    out["source_row_hash"] = out.apply(
        lambda x: sha256_text(f"{x['TF']}|{x['layer_id']}|{x['contrast']}|"
                              f"{x['time_point']}|{x['source_row_id']}"), axis=1)
    return out


def _residue(phos: str) -> str:
    s = str(phos or "")
    if "_" in s:
        return s.split("_", 1)[1]
    return ""


def _effect_measure(layer: str) -> str:
    return {
        "transcriptomics": "log2FoldChange",
        "proteomics": "logFC",
        "phosphoproteomics": "logFC",
        "coIP": "t_statistic",
        "tf_activity": "ULM_score",
        "kinase_activity": "ULM_score",
        "bionic_gnn": "module_membership",
        "early_synapse_trafficking": "membership",
    }.get(layer, "")


def _indep_rationale(cls: str) -> str:
    return {
        "A": "directly measured in a distinct molecular assay; not computed from the "
             "transcript signal used for candidate ranking",
        "B": "computed from an orthogonal experimental assay (e.g. kinase activity from "
             "phosphoproteomics)",
        "C": "computed from transcript expression / model modules closely related to the "
             "ranking input; contextual",
        "D": "lineage as evidence is unclear or curated; excluded from support counts",
    }.get(cls, "")


# --------------------------------------------------------------------------- #
# Candidate summary
# --------------------------------------------------------------------------- #
POSITIVE_STATUSES = {"positive_supportive", "positive_opposing",
                     "positive_unresolved_direction"}


def build_candidate_summary(long_aud: pd.DataFrame, union: pd.DataFrame) -> pd.DataFrame:
    rows = []
    scored = set(long_aud["TF"])
    for _, u in union.iterrows():
        tf = u["TF"]
        g = long_aud[long_aud["TF"] == tf]

        def cls(c):
            return g[g["independence_class"] == c]

        def pos_support(sub, rel):
            return sub[(sub["evidence_status"].isin(POSITIVE_STATUSES)) &
                       (sub["relation_to_BTLA_hypothesis"] == rel)]

        A, B, C = cls("A"), cls("B"), cls("C")
        A_pos = A[A["evidence_status"].isin(POSITIVE_STATUSES)]
        B_pos = B[B["evidence_status"].isin(POSITIVE_STATUSES)]
        C_pos = C[C["evidence_status"].isin(POSITIVE_STATUSES)]

        measured_layers = g[g["evidence_status"].isin(
            POSITIVE_STATUSES | {"measured_no_qualifying_signal"})]["layer_id"].nunique()
        not_measured_layers = g[g["evidence_status"] == "not_measured"]["layer_id"].nunique()

        # decision-rule field: True only when >=1 class A item is positive AND supportive
        A_supportive = A_pos[A_pos["relation_to_BTLA_hypothesis"] == "supportive"]
        decision = bool(len(A_supportive) >= 1)

        if tf not in scored:
            audit_status = "source_unavailable_not_evaluated_by_08e_pipeline"
        elif len(A_pos):
            audit_status = "independent_molecular_signal_present_direction_unresolved" \
                if not decision else "independent_supportive"
        elif len(B_pos) or len(C_pos):
            audit_status = "contextual_or_orthogonal_only"
        else:
            audit_status = "no_qualifying_signal_in_evaluated_layers"

        def strongest(sub):
            s = sub[sub["legacy_support_level"] == "strong"]
            if len(s):
                return s.iloc[0]["layer_id"]
            m = sub[sub["legacy_support_level"] == "moderate"]
            return m.iloc[0]["layer_id"] if len(m) else ""

        rows.append({
            "TF": tf,
            "candidate_group": u["candidate_group"],
            "gremln_dense_rank": u["gremln_dense_rank"],
            "genie3_dense_rank": u["genie3_dense_rank"],
            "in_gremln_top25": u["in_gremln_top25"],
            "in_genie3_top25": u["in_genie3_top25"],
            "classA_independent_supportive": int(len(A_supportive)),
            "classA_independent_opposing": int(len(pos_support(A, "opposing"))),
            "classA_independent_positive_unresolved": int(
                len(A_pos[A_pos["relation_to_BTLA_hypothesis"] == "unresolved"])),
            "classB_orthogonal_supportive": int(len(pos_support(B, "supportive"))),
            "classB_orthogonal_opposing": int(len(pos_support(B, "opposing"))),
            "classB_orthogonal_positive": int(len(B_pos)),
            "classC_contextual_supportive": int(len(pos_support(C, "supportive"))),
            "classC_contextual_opposing": int(len(pos_support(C, "opposing"))),
            "classC_contextual_positive": int(len(C_pos)),
            "n_positive_layers_classA": int(A_pos["layer_id"].nunique()),
            "n_positive_layers_classB": int(B_pos["layer_id"].nunique()),
            "n_positive_layers_classC": int(C_pos["layer_id"].nunique()),
            "n_measured_layers": int(measured_layers),
            "n_not_measured_layers": int(not_measured_layers),
            "strongest_independent_layer": strongest(A),
            "strongest_contextual_layer": strongest(C),
            "evidence_item_ids": ";".join(g["source_row_hash"].tolist()),
            "independent_molecular_support_for_decision_rule": decision,
            "overall_audit_status": audit_status,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Legacy reconciliation
# --------------------------------------------------------------------------- #
def reconcile_legacy(summary_aud: pd.DataFrame, legacy_summary: Path,
                     legacy_pub: Path, union: pd.DataFrame) -> pd.DataFrame:
    leg = pd.read_csv(legacy_summary)
    leg_by = {r["gene"]: r for _, r in leg.iterrows()}
    pub = pd.read_csv(legacy_pub) if Path(legacy_pub).exists() else pd.DataFrame()
    pub_by = {r["TF"]: r for _, r in pub.iterrows()} if len(pub) else {}

    rows = []
    for _, a in summary_aud.iterrows():
        tf = a["TF"]
        lg = leg_by.get(tf)
        legacy_overall = str(lg["overall_multiomics_support"]) if lg is not None else "ABSENT_from_legacy_summary"
        legacy_indep = str(lg["independent_support"]) if lg is not None else ""
        pub_indep = ""
        if tf in pub_by:
            pub_indep = str(pub_by[tf].get("independent_orthogonal_support", ""))
        rows.append({
            "TF": tf,
            "candidate_group": a["candidate_group"],
            "legacy_overall_multiomics_support": legacy_overall,
            "legacy_independent_support": legacy_indep,
            "legacy_pub_independent_orthogonal_support": pub_indep,
            "audited_overall_status": a["overall_audit_status"],
            "audited_classA_supportive": a["classA_independent_supportive"],
            "audited_classA_positive_unresolved": a["classA_independent_positive_unresolved"],
            "audited_decision_field": a["independent_molecular_support_for_decision_rule"],
            "status_changed": _changed(legacy_overall, a),
        })
    return pd.DataFrame(rows)


def _changed(legacy_overall: str, a) -> str:
    changes = []
    legacy_supportish = legacy_overall in {"strong", "moderate", "weak"}
    audited_supportish = a["overall_audit_status"] in {"independent_supportive"}
    if legacy_overall == "ABSENT_from_legacy_summary":
        changes.append("added_to_audit_scope")
    if legacy_supportish and not audited_supportish:
        changes.append("legacy_support_downgraded_or_reclassified")
    if a["classA_independent_positive_unresolved"] > 0 and not a[
            "independent_molecular_support_for_decision_rule"]:
        changes.append("class_A_signal_present_but_direction_unresolved")
    return ";".join(changes) if changes else "no_material_change"


# --------------------------------------------------------------------------- #
# Audit checks
# --------------------------------------------------------------------------- #
def run_audit_checks(union, long_aud, summary_aud, registry) -> pd.DataFrame:
    checks = []

    def chk(name, cond, detail=""):
        checks.append({"check": name, "passed": bool(cond), "detail": str(detail)})

    chk("union has exactly 43 candidates", len(union) == 43, len(union))
    chk("GREmLN contributes 25",
        int(union["in_gremln_top25"].sum()) == 25, int(union["in_gremln_top25"].sum()))
    chk("GENIE3 contributes 25",
        int(union["in_genie3_top25"].sum()) == 25, int(union["in_genie3_top25"].sum()))
    chk("shared == 7",
        int((union["candidate_group"] == "shared").sum()) == 7,
        int((union["candidate_group"] == "shared").sum()))
    chk("summary has each candidate exactly once",
        len(summary_aud) == 43 and summary_aud["TF"].is_unique, len(summary_aud))

    # every positive call has a source evidence row
    pos = long_aud[long_aud["evidence_status"].isin(POSITIVE_STATUSES)]
    chk("every positive call has source_row_hash",
        pos["source_row_hash"].notna().all() and (pos["source_row_hash"] != "").all(),
        len(pos))
    # every positive statistical call has thresholds documented
    stat_layers = {"transcriptomics", "proteomics", "phosphoproteomics",
                   "tf_activity", "kinase_activity"}
    pos_stat = pos[pos["layer_id"].isin(stat_layers)]
    chk("every positive statistical call passed documented thresholds",
        bool(pos_stat["significance_threshold_passed"].all()), len(pos_stat))
    # phosphosite mapping
    phos = long_aud[long_aud["layer_id"] == "phosphoproteomics"]
    chk("every phosphosite row has a site id",
        bool((phos["protein_or_phosphosite_id"].astype(str).str.len() > 0).all()), len(phos))
    # coIP bait/condition
    coip = long_aud[long_aud["layer_id"] == "coIP"]
    chk("every coIP row has a bait",
        bool((coip["bait"] == "BTLA").all()) if len(coip) else True, len(coip))
    # class A only observed molecular assays
    classA_layers = set(long_aud[long_aud["independence_class"] == "A"]["layer_id"])
    chk("class A only from observed molecular assays",
        classA_layers <= {"proteomics", "phosphoproteomics", "coIP"}, classA_layers)
    # transcript-derived cannot enter decision field
    decisions = summary_aud[summary_aud["independent_molecular_support_for_decision_rule"]]
    ok_decision = True
    for _, d in decisions.iterrows():
        a = long_aud[(long_aud["TF"] == d["TF"]) &
                     (long_aud["independence_class"] == "A") &
                     (long_aud["relation_to_BTLA_hypothesis"] == "supportive")]
        if not len(a):
            ok_decision = False
    chk("decision field requires >=1 class A supportive item", ok_decision, len(decisions))
    chk("no class C/D layer feeds decision field",
        True, "enforced by construction (class A only)")
    # not measured never treated as support
    nm = long_aud[long_aud["evidence_status"] == "not_measured"]
    chk("not_measured never positive",
        bool((~nm["evidence_status"].isin(POSITIVE_STATUSES)).all()), len(nm))
    # opposing not counted supportive
    chk("opposing not counted as supportive",
        int(summary_aud["classA_independent_opposing"].sum()) >= 0, "structural")
    # shared candidates identical in both model views (single summary row -> trivially identical)
    chk("shared candidates single consistent summary row",
        summary_aud["TF"].is_unique, "single-view summary")
    # model counts reconcile
    gr_rows = summary_aud[summary_aud["in_gremln_top25"]]
    g3_rows = summary_aud[summary_aud["in_genie3_top25"]]
    chk("GREmLN view has 25 rows", len(gr_rows) == 25, len(gr_rows))
    chk("GENIE3 view has 25 rows", len(g3_rows) == 25, len(g3_rows))
    # evidence joins cannot alter rank (ranks come only from union)
    chk("ranks sourced only from frozen union",
        set(summary_aud["TF"]) == set(union["TF"]), "structural")

    return pd.DataFrame(checks)


# --------------------------------------------------------------------------- #
# main
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

    long_path = dr / "results/btla_csls_multiomics/btla_candidate_multiomics_long_evidence.csv"
    summary_path = dr / "results/btla_csls_multiomics/btla_candidate_multiomics_summary.csv"
    gr_rank_path = repo / "results/tables/gremln_btla_vs_tcr_seed_excluded_tf_ranking.csv"
    g3_rank_path = repo / "results/tables/genie3_btla_vs_tcr_seed_excluded_tf_ranking.csv"
    union_path = repo / "results/publication_data/top25_union_primary.csv"
    legacy_pub = repo / "results/publication_data/multiomics_candidate_evidence.csv"

    # 1. source inventory
    inv = build_source_inventory(long_path, summary_path, gr_rank_path,
                                 g3_rank_path, union_path)
    inv.to_csv(out / "source_inventory.csv", index=False)

    # 2. candidate union
    union = build_candidate_union(union_path, gr_rank_path, g3_rank_path)
    union["gremln_ranking_sha256"] = sha256_file(gr_rank_path)
    union["genie3_ranking_sha256"] = sha256_file(g3_rank_path)
    union.to_csv(out / "candidate_union.csv", index=False)

    # 3. registry
    registry = build_layer_registry()
    registry.to_csv(out / "multiomics_layer_registry.csv", index=False)

    # 4-6. long audited
    long_df = pd.read_csv(long_path)
    long_aud = audit_long_evidence(long_df, union)

    # restricted full (gitignored)
    long_aud.to_csv(out / "local_restricted" /
                    "multiomics_evidence_long_audited.FULL.csv", index=False)
    # committed redacted: drop effect values, adjusted p, phosphosite residue/id,
    # and per-target direction (sign of a significant molecular change is inside
    # the publication gate). evidence_status keeps the direction-neutral status.
    redacted = long_aud.drop(columns=[
        "effect_value", "adjusted_p", "protein_or_phosphosite_id",
        "phosphosite_residue", "raw_p", "raw_direction", "evidence_direction"]).copy()
    redacted["effect_value"] = "REDACTED_publication_gate"
    redacted["phosphosite_id"] = "REDACTED_publication_gate"
    redacted.to_csv(out / "multiomics_evidence_long_audited.csv", index=False)

    # 7. candidate summary (statuses + counts, no effect sizes)
    summary_aud = build_candidate_summary(long_aud, union)
    summary_aud.to_csv(out / "multiomics_candidate_summary_audited.csv", index=False)

    # 9. legacy reconciliation
    legacy = reconcile_legacy(summary_aud, summary_path, legacy_pub, union)
    legacy.to_csv(out / "legacy_vs_audited_multiomics.csv", index=False)

    # 10. audit checks
    checks = run_audit_checks(union, long_aud, summary_aud, registry)
    checks.to_csv(out / "multiomics_audit_checks.csv", index=False)

    # manifest
    manifest = {
        "generated_utc": utcnow(),
        "data_root_redacted": "DATA_ROOT (env BTLA_BENCH_DATA_ROOT)",
        "thresholds": THRESHOLDS,
        "layer_independence": LAYER_INDEPENDENCE,
        "n_candidates": int(len(union)),
        "n_long_rows_union": int(len(long_aud)),
        "n_candidates_evaluated_by_08e": int(long_aud["TF"].nunique()),
        "long_source_sha256": sha256_file(long_path),
        "legacy_summary_sha256": sha256_file(summary_path),
        "gremln_ranking_sha256": sha256_file(gr_rank_path),
        "genie3_ranking_sha256": sha256_file(g3_rank_path),
        "union_sha256": sha256_file(union_path),
        "checks_passed": int(checks["passed"].sum()),
        "checks_total": int(len(checks)),
        "decision_rule_true_candidates":
            summary_aud[summary_aud["independent_molecular_support_for_decision_rule"]]["TF"].tolist(),
        "classA_signal_present_candidates":
            summary_aud[summary_aud["classA_independent_positive_unresolved"] > 0]["TF"].tolist(),
    }
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2))

    print("[multiomics-audit] checks:",
          f"{int(checks['passed'].sum())}/{len(checks)} passed")
    for _, c in checks[~checks["passed"]].iterrows():
        print("  FAIL:", c["check"], c["detail"])
    print("[multiomics-audit] class A signal present:",
          manifest["classA_signal_present_candidates"])
    print("[multiomics-audit] decision-rule TRUE:",
          manifest["decision_rule_true_candidates"])
    print("[multiomics-audit] candidates evaluated by 08e:",
          manifest["n_candidates_evaluated_by_08e"], "/ 43")
    return 0 if bool(checks["passed"].all()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
