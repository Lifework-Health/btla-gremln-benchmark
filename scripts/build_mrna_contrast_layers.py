#!/usr/bin/env python3
"""Candidate mRNA layer across the full experimental sequence (3 contrasts).

Refinement of the contextual mRNA layer. Instead of only the incremental
BTLA+TCR-vs-TCR contrast, the transcriptional responsiveness of each of the 43
candidates is described across three prespecified contrasts and reported
separately (never collapsed into one strongest result):

  C1  TCR_activation            = TCR vs unstimulated control  (source: TCRvsUC)
  C2  Incremental_BTLA          = BTLA+TCR vs TCR              (source: BTLAvsTCR)
  C3  Combined_BTLA_TCR_vs_UC   = BTLA+TCR vs unstimulated control

Contrast 3 does NOT exist in the source differential-expression output. Both
source workbooks (DEG_statistics.xlsx long-form and DEG_wide_statistics.xlsx)
contain only BTLAvsTCR and TCRvsUC, each at 1h/4h/24h, and the Methods confirm
only pairwise per-timepoint designs were fitted. It CANNOT be regenerated from
the fitted model because the DESeq2 dds object / raw counts / coldata are not
present in either data root. It is therefore represented by the explicit,
non-fabricated state `source_unavailable_pending_model`. Contrast 3 is NEVER
approximated by adding fold changes or combining p-values from C1 and C2.

Source-native transcriptomics threshold (Methods.docx):
    adjusted p < 0.05  AND  |log2 fold change| > log2(1.5)

Candidate mRNA remains CONTEXTUAL evidence: it does not enter the independent
molecular corroboration count (protein + phosphosite only).

Redaction: statuses, qualifying contrasts/timepoints and direction sign (up/down)
are committed; numeric log2FC / adjusted p go only to a git-ignored restricted
file.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

PADJ_MAX = 0.05
RNA_LOG2 = math.log2(1.5)  # 0.585 source-native transcriptomics effect threshold
TPS = ["1h", "4h", "24h"]

# prespecified contrasts, in experimental-sequence order
CONTRASTS = [
    {"key": "TCR_activation", "label": "TCR activation",
     "biological": "TCR vs unstimulated control", "source_prefix": "TCRvsUC",
     "available": True},
    {"key": "Incremental_BTLA", "label": "Incremental BTLA",
     "biological": "BTLA+TCR vs TCR", "source_prefix": "BTLAvsTCR",
     "available": True},
    {"key": "Combined_BTLA_TCR_vs_UC", "label": "Combined BTLA+TCR state",
     "biological": "BTLA+TCR vs unstimulated control", "source_prefix": None,
     "available": False},
]
PENDING = "source_unavailable_pending_model"


def data_root() -> Path:
    return Path(os.environ.get("BTLA_BENCH_DATA_ROOT",
                               "/mnt/R0/Projects/POIAZ/gmorgan/gremln-tcells")).expanduser()


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def status_of(effect, padj) -> str:
    """Six-state source-native call for one measured value."""
    if pd.isna(effect) and pd.isna(padj):
        return "not_measured"
    if pd.isna(effect) or pd.isna(padj):
        return "failed_qc"
    if (padj < PADJ_MAX) and (abs(effect) > RNA_LOG2):
        return "measured_qualifying"
    return "measured_no_qualifying_signal"


def direction_of(effect) -> str:
    if pd.isna(effect):
        return ""
    return "up" if effect > 0 else "down"


def main() -> int:
    dr = data_root()
    repo = repo_root()
    out = repo / "results/multiomics/audit_v2"
    rst = out / "local_restricted"
    rst.mkdir(parents=True, exist_ok=True)

    union = pd.read_csv(repo / "results/publication_data/top25_union_primary.csv")
    assert len(union) == 43, f"expected 43 candidates, got {len(union)}"
    group = dict(zip(union["TF"], union["status"]))

    deg = pd.read_excel(
        dr / "data/multiomics/Transcriptomics/@Transcriptomics/Tables/DEG_wide_statistics.xlsx")
    deg_by = {str(r["gene"]): r for _, r in deg.iterrows()}
    universe = set(deg["gene"].astype(str))

    long_rows = []          # committed (statuses + direction sign)
    restricted_rows = []    # git-ignored (numeric)
    for tf in union["TF"]:
        in_universe = tf in universe
        r = deg_by.get(tf)
        for c in CONTRASTS:
            for tp in TPS:
                if not c["available"]:
                    long_rows.append({
                        "TF": tf, "candidate_group": group[tf],
                        "contrast_key": c["key"], "contrast_label": c["label"],
                        "contrast_biological": c["biological"], "timepoint": tp,
                        "status": PENDING, "direction": "",
                        "qualifying_source": False,
                        "source_prefix": "", "in_assay_universe": in_universe})
                    continue
                if not in_universe:
                    long_rows.append({
                        "TF": tf, "candidate_group": group[tf],
                        "contrast_key": c["key"], "contrast_label": c["label"],
                        "contrast_biological": c["biological"], "timepoint": tp,
                        "status": "absent_from_assay_universe", "direction": "",
                        "qualifying_source": False,
                        "source_prefix": c["source_prefix"], "in_assay_universe": False})
                    continue
                eff = r.get(f"{c['source_prefix']}_{tp}_log2FoldChange")
                padj = r.get(f"{c['source_prefix']}_{tp}_padj")
                st = status_of(eff, padj)
                qual = st == "measured_qualifying"
                long_rows.append({
                    "TF": tf, "candidate_group": group[tf],
                    "contrast_key": c["key"], "contrast_label": c["label"],
                    "contrast_biological": c["biological"], "timepoint": tp,
                    "status": st,
                    "direction": direction_of(eff) if qual else "",
                    "qualifying_source": qual,
                    "source_prefix": c["source_prefix"], "in_assay_universe": True})
                restricted_rows.append({
                    "TF": tf, "contrast_key": c["key"], "timepoint": tp,
                    "log2FoldChange": eff, "padj": padj, "status": st})

    long = pd.DataFrame(long_rows)
    long.to_csv(out / "mrna_contrast_layers_long_v2.csv", index=False)
    pd.DataFrame(restricted_rows).to_csv(
        rst / "mrna_contrast_layers_values.RESTRICTED.csv", index=False)

    # ---- candidate-level summary ----
    srows = []
    for tf in union["TF"]:
        rec = {"TF": tf, "candidate_group": group[tf]}
        qual_contexts = []
        for c in CONTRASTS:
            sub = long[(long["TF"] == tf) & (long["contrast_key"] == c["key"])]
            if not c["available"]:
                rec[f"{c['key']}__status"] = PENDING
                rec[f"{c['key']}__qualifies"] = False
                rec[f"{c['key']}__qual_timepoints"] = ""
                rec[f"{c['key']}__directions"] = ""
                continue
            qh = sub[sub["qualifying_source"]]
            qualifies = len(qh) > 0
            tps = sorted(qh["timepoint"], key=TPS.index)
            dirs = sorted(set(qh["direction"]))
            rec[f"{c['key']}__qualifies"] = qualifies
            rec[f"{c['key']}__qual_timepoints"] = ";".join(tps)
            rec[f"{c['key']}__directions"] = ";".join(dirs)
            # layer status when not qualifying
            if qualifies:
                rec[f"{c['key']}__status"] = "measured_qualifying"
                for tp in tps:
                    d = qh[qh["timepoint"] == tp]["direction"].iloc[0]
                    qual_contexts.append(f"{c['key']}@{tp}:{d}")
            elif (sub["status"] == "absent_from_assay_universe").all():
                rec[f"{c['key']}__status"] = "absent_from_assay_universe"
            elif (sub["status"] == "not_measured").all():
                rec[f"{c['key']}__status"] = "not_measured"
            elif (sub["status"] == "failed_qc").any() and \
                 not (sub["status"] == "measured_no_qualifying_signal").any():
                rec[f"{c['key']}__status"] = "failed_qc"
            else:
                rec[f"{c['key']}__status"] = "measured_no_qualifying_signal"

        n_qual_contrasts = sum(bool(rec[f"{c['key']}__qualifies"])
                               for c in CONTRASTS if c["available"])
        rec["n_qualifying_contrasts"] = n_qual_contrasts
        rec["transcriptionally_responsive_any_context"] = n_qual_contrasts > 0
        # the summary flag must always be accompanied by source contrast + timepoint
        rec["responsive_contexts"] = ";".join(qual_contexts)
        srows.append(rec)
    summ = pd.DataFrame(srows)
    summ.to_csv(out / "mrna_contrast_layers_summary_v2.csv", index=False)

    # ---- counts ----
    def qn(key):
        return int(summ[f"{key}__qualifies"].sum())
    counts = {
        "n_candidates": len(summ),
        "qualify_C1_TCR_activation": qn("TCR_activation"),
        "qualify_C2_Incremental_BTLA": qn("Incremental_BTLA"),
        "qualify_C3_Combined_pending": PENDING,  # not measurable
        "qualify_in_more_than_one_available_contrast": int(
            (summ["n_qualifying_contrasts"] > 1).sum()),
        "qualify_in_any_available_context": int(
            summ["transcriptionally_responsive_any_context"].sum()),
        "available_contrasts": ["TCR_activation", "Incremental_BTLA"],
        "unavailable_contrasts": {"Combined_BTLA_TCR_vs_UC": PENDING},
    }
    by_group = {}
    for g in ["GREmLN_specific", "shared", "GENIE3_specific"]:
        s = summ[summ["candidate_group"] == g]
        by_group[g] = {
            "n": len(s),
            "C1": int(s["TCR_activation__qualifies"].sum()),
            "C2": int(s["Incremental_BTLA__qualifies"].sum()),
            "any": int(s["transcriptionally_responsive_any_context"].sum()),
        }
    counts["by_group"] = by_group
    (rst / "mrna_contrast_layers_counts.json").write_text(json.dumps(counts, indent=2))

    # lists (committed via report elsewhere; printed here)
    c1_tfs = sorted(summ[summ["TCR_activation__qualifies"]]["TF"])
    c2_tfs = sorted(summ[summ["Incremental_BTLA__qualifies"]]["TF"])
    both = sorted(set(c1_tfs) & set(c2_tfs))
    any_tfs = sorted(summ[summ["transcriptionally_responsive_any_context"]]["TF"])

    print("=== candidate mRNA across the experimental sequence (43 candidates) ===")
    print("C1 TCR activation (TCR vs UC) qualifying:", counts["qualify_C1_TCR_activation"])
    print("   ", c1_tfs)
    print("C2 Incremental BTLA (BTLA+TCR vs TCR) qualifying:",
          counts["qualify_C2_Incremental_BTLA"])
    print("   ", c2_tfs)
    print("C3 Combined BTLA+TCR vs UC:", PENDING, "(not in source; not regenerable)")
    print("qualify in >1 available contrast:",
          counts["qualify_in_more_than_one_available_contrast"], "->", both)
    print("qualify in ANY available context:",
          counts["qualify_in_any_available_context"], "->", any_tfs)
    print("by group:", json.dumps(by_group))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
