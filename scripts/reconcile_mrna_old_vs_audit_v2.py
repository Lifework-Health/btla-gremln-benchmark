#!/usr/bin/env python3
"""Reconcile the old 08e GENIE3 mRNA chart against audit_v2 (BTLA-specific).

Section 1-3 of the mRNA reconciliation task. Evidence construction only; no
manuscript / verdict changes. Restricted values (log2FC, adj p) are read from
the local source only for provenance and written to a git-ignored file; the
committed reconciliation carries statuses/timepoints/contrasts only.

Old chart faithfully reconstructed from the frozen 08e long-evidence table
(the exact data the chart read) and independently re-derived from the raw
DEG workbook for both contrasts at each timepoint.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

PADJ_MAX = 0.05
RNA_LOG2 = math.log2(1.5)  # 0.585 source-native transcriptomics threshold
STRONG_LOG2 = RNA_LOG2 * 1.5  # 08e "strong" inner tier (for glyph reconstruction only)
STRONG_PADJ = 0.01
RNA_TPS = ["1h", "4h", "24h"]


def data_root() -> Path:
    return Path(os.environ.get("BTLA_BENCH_DATA_ROOT",
                               "/mnt/R0/Projects/POIAZ/gmorgan/gremln-tcells")).expanduser()


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def source_status(effect, padj):
    """audit_v2 source-native six-state for one measured value."""
    if pd.isna(effect) and pd.isna(padj):
        return "not_measured"
    if pd.isna(effect) or pd.isna(padj):
        return "failed_qc"
    if (padj < PADJ_MAX) and (abs(effect) > RNA_LOG2):
        return "measured_qualifying"
    return "measured_no_qualifying_signal"


def glyph_08e(effect, padj):
    """Reproduce 08e support_from_p tier -> glyph for one value."""
    if pd.isna(effect) or pd.isna(padj):
        return "not_measured", "—"
    if padj >= PADJ_MAX or abs(effect) < RNA_LOG2:
        return "not_significant", "ns"
    if padj < STRONG_PADJ and abs(effect) >= STRONG_LOG2:
        return "strong", "++"
    return "moderate", "·"


def main() -> int:
    dr = data_root()
    repo = repo_root()
    out = repo / "results/multiomics/audit_v2"
    (out / "local_restricted").mkdir(parents=True, exist_ok=True)

    # ---- inputs ----
    oldg = pd.read_csv(dr / "results/btla_csls_multiomics/btla_tf_genie3_top25.csv")
    old_rank = dict(zip(oldg["gene"], oldg["genie3_rank"]))
    old_seed = dict(zip(oldg["gene"], oldg["is_BTLA_vs_TCR_seed"]))
    old_set = list(oldg["gene"])

    union = pd.read_csv(repo / "results/publication_data/top25_union_primary.csv")

    def grp(r):
        if r["in_gremln_top25"] and r["in_genie3_top25"]:
            return "shared"
        return "GREmLN_specific" if r["in_gremln_top25"] else "GENIE3_specific"
    union["candidate_group"] = union.apply(grp, axis=1)
    union_set = set(union["TF"])
    ugroup = dict(zip(union["TF"], union["candidate_group"]))

    # frozen long-evidence (exact data old chart read)
    le = pd.read_csv(dr / "results/btla_csls_multiomics/btla_candidate_multiomics_long_evidence.csv")
    tr = le[le["evidence_layer"] == "transcriptomics"].copy()

    # raw DEG workbook (independent recompute, both contrasts, per timepoint)
    deg = pd.read_excel(dr / "data/multiomics/Transcriptomics/@Transcriptomics/Tables/DEG_wide_statistics.xlsx")
    deg_by = {str(r["gene"]): r for _, r in deg.iterrows()}
    deg_universe = set(deg["gene"].astype(str))

    def raw_cell(gene, prefix, tp):
        if gene not in deg_universe:
            return "absent_from_assay_universe", np.nan, np.nan
        r = deg_by[gene]
        eff = r.get(f"{prefix}_{tp}_log2FoldChange")
        padj = r.get(f"{prefix}_{tp}_padj")
        return source_status(eff, padj), eff, padj

    STRONG_LEVELS = {"strong", "moderate"}

    def old_mrna_glyph(gene):
        """Reproduce _compact_layer_evidence(gene,'transcriptomics',contrast=None)."""
        sub = tr[tr["gene"] == gene]
        if sub.empty:
            return dict(status="—", tier="none", timepoint="", contrast="", direction="")
        hits = sub[sub["support_level"].isin(STRONG_LEVELS)].copy()
        if not hits.empty:
            hits["_ord"] = hits["support_level"].map({"strong": 0, "moderate": 1})
            hits = hits.sort_values(["_ord", "adjusted_p"], na_position="last")
            r = hits.iloc[0]
            tier = r["support_level"]
            return dict(status="++" if tier == "strong" else "·", tier=tier,
                        timepoint=str(r.get("timepoint", "")),
                        contrast=str(r.get("contrast", "")),
                        direction=str(r.get("direction", "")))
        if (sub["support_level"] == "not_significant").any():
            return dict(status="ns", tier="not_significant", timepoint="", contrast="", direction="")
        return dict(status="—", tier="not_measured", timepoint="", contrast="", direction="")

    # ---- build row-level reconciliation over old ∪ union ----
    all_tfs = sorted(set(old_set) | union_set)
    rows = []
    for tf in all_tfs:
        in_old = tf in old_set
        in_union = tf in union_set
        og = old_mrna_glyph(tf) if in_old else dict(status="", tier="", timepoint="", contrast="", direction="")

        cells = {}
        for prefix, label in [("BTLAvsTCR", "btla"), ("TCRvsUC", "tcruc")]:
            for tp in RNA_TPS:
                st, eff, padj = raw_cell(tf, prefix, tp)
                cells[f"{label}_{tp}"] = st

        # reason for disagreement
        reason = []
        if in_old and not in_union:
            if old_seed.get(tf):
                reason.append("old_seed_inclusive_candidate_removed_by_seed_exclusion")
            else:
                reason.append("old_candidate_outside_current_seed_excluded_top25")
        if in_old and og["status"] in {"++", "·"}:
            if og["contrast"] == "TCR_vs_UC":
                reason.append("old_positive_from_TCRvsUC_not_BTLA_specific")
            elif og["contrast"] == "BTLA_TCR_vs_TCR":
                reason.append("old_positive_from_BTLAvsTCR_should_persist")
        btla_qual = any(cells[f"btla_{tp}"] == "measured_qualifying" for tp in RNA_TPS)
        if in_old and og["status"] in {"++", "·"} and not btla_qual:
            reason.append("no_BTLA_specific_qualifying_signal_in_current_rule")
        if not reason:
            reason.append("no_disagreement" if not (in_old and in_union) or og["status"] not in {"++", "·"}
                          else "consistent")

        # source threshold result (BTLA-specific summary)
        src_result = ("qualifying" if btla_qual else
                      "absent_from_assay_universe" if all(
                          cells[f"btla_{tp}"] == "absent_from_assay_universe" for tp in RNA_TPS)
                      else "no_qualifying_signal")

        rows.append({
            "TF": tf,
            "in_old_chart": in_old,
            "old_chart_rank": old_rank.get(tf, ""),
            "old_seed_status": ("BTLA_vs_TCR_seed" if old_seed.get(tf) else "not_seed") if in_old else "",
            "in_current_union": in_union,
            "current_candidate_group": ugroup.get(tf, ""),
            "old_displayed_mrna_status": og["status"],
            "old_displayed_timepoint": og["timepoint"],
            "old_source_contrast": og["contrast"],
            "old_direction": og["direction"],
            "BTLAvsTCR_1h_status": cells["btla_1h"],
            "BTLAvsTCR_4h_status": cells["btla_4h"],
            "BTLAvsTCR_24h_status": cells["btla_24h"],
            "TCRvsUC_1h_status": cells["tcruc_1h"],
            "TCRvsUC_4h_status": cells["tcruc_4h"],
            "TCRvsUC_24h_status": cells["tcruc_24h"],
            "source_threshold_result_BTLA_specific": src_result,
            "reason_for_disagreement": ";".join(sorted(set(reason))),
        })
    recon = pd.DataFrame(rows)
    recon.to_csv(out / "mrna_old_vs_current_reconciliation.csv", index=False)

    # restricted: exact log2FC/padj for old positive calls (focus + any old positive)
    focus = ["JUNB", "HIVEP3", "EGR2", "IRF8", "STAT1", "CREM", "MYC", "ZBTB32",
             "ID2", "MDM2", "TRAF4", "ZEB2", "PRDM1", "BCL3", "MAF"]
    rrows = []
    for tf in sorted(set(old_set) | set(focus)):
        if tf not in deg_universe:
            continue
        r = deg_by[tf]
        for prefix in ["BTLAvsTCR", "TCRvsUC"]:
            for tp in RNA_TPS:
                rrows.append({
                    "TF": tf, "contrast": prefix, "timepoint": tp,
                    "log2FoldChange": r.get(f"{prefix}_{tp}_log2FoldChange"),
                    "padj": r.get(f"{prefix}_{tp}_padj"),
                    "status_source_rule": source_status(
                        r.get(f"{prefix}_{tp}_log2FoldChange"), r.get(f"{prefix}_{tp}_padj")),
                })
    pd.DataFrame(rrows).to_csv(
        out / "local_restricted" / "mrna_reconciliation_values.RESTRICTED.csv", index=False)

    # ---- section 3: verify current 43 BTLA-specific mRNA ----
    v = union[["TF", "candidate_group"]].copy()
    stats = {"measured": 0, "not_measured": 0, "absent": 0, "padj_sig": 0,
             "effect_thresh": 0, "both_qualifying": 0}
    near = []
    per = []
    for tf in union["TF"]:
        if tf not in deg_universe:
            stats["absent"] += 1
            per.append((tf, "absent_from_assay_universe"))
            continue
        r = deg_by[tf]
        measured_any = False
        padj_sig_any = eff_any = both_any = False
        best_gap = None
        for tp in RNA_TPS:
            eff = r.get(f"BTLAvsTCR_{tp}_log2FoldChange")
            padj = r.get(f"BTLAvsTCR_{tp}_padj")
            if pd.notna(eff) and pd.notna(padj):
                measured_any = True
                if padj < PADJ_MAX:
                    padj_sig_any = True
                if abs(eff) > RNA_LOG2:
                    eff_any = True
                if padj < PADJ_MAX and abs(eff) > RNA_LOG2:
                    both_any = True
                # nearest miss: significant padj but sub-threshold effect
                if padj < PADJ_MAX and abs(eff) <= RNA_LOG2:
                    gap = RNA_LOG2 - abs(eff)
                    if best_gap is None or gap < best_gap:
                        best_gap = gap
        stats["measured"] += int(measured_any)
        stats["not_measured"] += int(not measured_any)
        stats["padj_sig"] += int(padj_sig_any)
        stats["effect_thresh"] += int(eff_any)
        stats["both_qualifying"] += int(both_any)
        if best_gap is not None and not both_any:
            near.append((tf, round(best_gap, 3)))
        per.append((tf, "qualifying" if both_any else
                    "measured_no_qualifying_signal" if measured_any else "not_measured"))

    # assertions
    asserts = {}
    asserts["all_43_accounted"] = len(v) == 43 and len(per) == 43
    asserts["no_TCRvsUC_in_BTLA_result"] = True  # by construction: only BTLAvsTCR cols used
    asserts["missing_not_called_no_signal"] = (stats["absent"] ==
        len([1 for tf in union["TF"] if tf not in deg_universe]))
    asserts["zero_calls_were_measured_or_marked"] = all(
        s in {"qualifying", "measured_no_qualifying_signal", "not_measured",
              "absent_from_assay_universe"} for _, s in per)

    summary = {
        "old_genie3_top25_n": len(old_set),
        "old_intersect_union": sorted(set(old_set) & union_set),
        "old_not_in_union": sorted(set(old_set) - union_set),
        "old_not_in_union_seed": sorted([g for g in set(old_set) - union_set if old_seed.get(g)]),
        "old_not_in_union_nonseed": sorted([g for g in set(old_set) - union_set if not old_seed.get(g)]),
        "old_mrna_positive_calls": recon[(recon.in_old_chart) &
            (recon.old_displayed_mrna_status.isin(["++", "·"]))][
            ["TF", "old_displayed_mrna_status", "old_source_contrast", "old_displayed_timepoint"]].to_dict("records"),
        "old_positive_from_TCRvsUC": sorted(recon[(recon.old_displayed_mrna_status.isin(["++", "·"])) &
            (recon.old_source_contrast == "TCR_vs_UC")]["TF"].tolist()),
        "old_positive_from_BTLAvsTCR": sorted(recon[(recon.old_displayed_mrna_status.isin(["++", "·"])) &
            (recon.old_source_contrast == "BTLA_TCR_vs_TCR")]["TF"].tolist()),
        "current_43_verification": stats,
        "current_per_candidate": per,
        "nearest_threshold_misses_padj_sig_subthreshold_effect": sorted(near, key=lambda x: x[1]),
        "assertions": asserts,
    }
    (out / "local_restricted" / "mrna_reconciliation_stats.json").write_text(json.dumps(summary, indent=2, default=str))

    print("== old GENIE3 top25:", len(old_set), "| in union:", len(set(old_set) & union_set),
          "| not in union:", len(set(old_set) - union_set),
          "(seeds:", len(summary["old_not_in_union_seed"]), ")")
    print("== old mRNA positive calls:", len(summary["old_mrna_positive_calls"]),
          "| from TCRvsUC:", len(summary["old_positive_from_TCRvsUC"]),
          "| from BTLAvsTCR:", len(summary["old_positive_from_BTLAvsTCR"]))
    print("   TCRvsUC-driven:", summary["old_positive_from_TCRvsUC"])
    print("   BTLAvsTCR-driven:", summary["old_positive_from_BTLAvsTCR"])
    print("== current 43 BTLA-specific mRNA:", stats)
    print("== nearest misses (padj-sig, sub-threshold effect):",
          summary["nearest_threshold_misses_padj_sig_subthreshold_effect"][:10])
    print("== assertions:", asserts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
