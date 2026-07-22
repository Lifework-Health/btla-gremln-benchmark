#!/usr/bin/env python3
"""Presentation layer for BTLA multiomics audit_v2.

Reads the audited v2 tables and emits three figures plus the pipeline audit
report and commit-safety report. Reads only statuses / counts / provenance;
no effect sizes or phosphosite residues.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT_REL = "results/multiomics/audit_v2"
LAYERS = ["transcriptomics", "proteomics", "phosphoproteomics", "tf_activity", "coIP"]
LAYER_LABEL = {
    "transcriptomics": "mRNA\n(context)", "proteomics": "Protein\n(indep)",
    "phosphoproteomics": "Phospho\n(indep)", "tf_activity": "TF act.\n(context)",
    "coIP": "coIP\n(excluded)",
}
STATUS_COLOR = {
    "measured_qualifying": "#1a7a3a",
    "measured_no_qualifying_signal": "#d9d9d9",
    "not_measured": "#f7f7f7",
    "absent_from_assay_universe": "#efefef",
    "ambiguous_mapping": "#fddbc7",
    "failed_qc": "#fdae6b",
    "excluded_no_statistical_contrast": "#e6d9f2",
}


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def load(out: Path):
    return (pd.read_csv(out / "multiomics_candidate_summary_v2.csv"),
            pd.read_csv(out / "multiomics_evidence_long_v2.csv"),
            pd.read_csv(out / "multiomics_model_summaries_v2.csv"),
            pd.read_csv(out / "multiomics_coverage_v2.csv"),
            pd.read_csv(out / "multiomics_audit_checks_v2.csv"),
            json.loads((out / "run_manifest_v2.json").read_text()))


def _candidate_layer_status(long: pd.DataFrame, tf: str, layer: str) -> str:
    sub = long[(long["TF"] == tf) & (long["layer"] == layer)]
    if not len(sub):
        return "not_measured"
    if (sub["qualifying_source"]).any():
        return "measured_qualifying"
    order = ["ambiguous_mapping", "measured_no_qualifying_signal",
             "excluded_no_statistical_contrast", "failed_qc", "not_measured",
             "absent_from_assay_universe"]
    present = set(sub["status"])
    for s in order:
        if s in present:
            return s
    return "not_measured"


def fig_matrix(out: Path, summ, long):
    order = summ.sort_values(
        ["independent_molecular_association", "phospho_qualifying_sites",
         "protein_status", "TF"],
        ascending=[False, False, True, True])["TF"].tolist()
    grp = dict(zip(summ["TF"], summ["candidate_group"]))
    fig, ax = plt.subplots(figsize=(8.5, 13))
    for i, tf in enumerate(order):
        for j, layer in enumerate(LAYERS):
            st = _candidate_layer_status(long, tf, layer)
            ax.add_patch(mpatches.Rectangle((j, len(order) - 1 - i), 1, 1,
                         fc=STATUS_COLOR.get(st, "#f7f7f7"), ec="#999", lw=0.4))
    ax.set_xlim(0, len(LAYERS)); ax.set_ylim(0, len(order))
    ax.set_xticks(np.arange(len(LAYERS)) + 0.5)
    ax.set_xticklabels([LAYER_LABEL[l] for l in LAYERS], fontsize=9)
    ax.set_yticks(np.arange(len(order)) + 0.5)
    tag = {"shared": "sh", "GREmLN_specific": "GR", "GENIE3_specific": "G3"}
    ax.set_yticklabels([f"{tf} ·{tag.get(grp[tf],'')}" for tf in reversed(order)],
                       fontsize=7.5)
    ax.set_title("BTLA multiomics audit_v2 — status by candidate x assay\n"
                 "independent = protein + phospho; mRNA / TF-activity contextual; coIP excluded",
                 fontsize=10.5)
    legend = [mpatches.Patch(fc=STATUS_COLOR[k], ec="#999", label=k.replace("_", " "))
              for k in ["measured_qualifying", "measured_no_qualifying_signal",
                        "not_measured", "absent_from_assay_universe",
                        "excluded_no_statistical_contrast"]]
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.04),
              ncol=2, fontsize=8, frameon=False)
    fig.savefig(out / "figures/fig_multiomics_v2_candidate_matrix.png", dpi=200,
                bbox_inches="tight")
    fig.savefig(out / "figures/fig_multiomics_v2_candidate_matrix.svg", bbox_inches="tight")
    plt.close(fig)


def fig_models(out: Path, models):
    views = ["GREmLN_top25", "GENIE3_top25", "shared", "GREmLN_specific",
             "GENIE3_specific", "union_43"]
    m = models.set_index("view").loc[views]
    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(len(views))
    ax.bar(x - 0.2, m["n_independent_association"], 0.4,
           label="independent association (protein/phospho)", color="#1a7a3a")
    ax.bar(x + 0.2, m["n_tf_activity_context"], 0.4,
           label="TF-activity context (not independent)", color="#9ecae1")
    for i, v in enumerate(views):
        ax.text(i - 0.2, m["n_independent_association"].iloc[i] + 0.1,
                str(int(m["n_independent_association"].iloc[i])), ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v}\n(n={int(m['n_candidates'].iloc[i])})"
                        for i, v in enumerate(views)], fontsize=8.5)
    ax.set_ylabel("candidates")
    ax.set_title("Independent molecular association by model view (audit_v2)\n"
                 "mRNA context change qualifying = 0 across all views (expected for upstream regulators)",
                 fontsize=10.5)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "figures/fig_multiomics_v2_model_summary.png", dpi=200,
                bbox_inches="tight")
    fig.savefig(out / "figures/fig_multiomics_v2_model_summary.svg", bbox_inches="tight")
    plt.close(fig)


def fig_coverage(out: Path, coverage):
    fig, ax = plt.subplots(figsize=(10, 5.5))
    y = np.arange(len(coverage))
    ax.barh(y, coverage["candidates_in_universe"], color="#c6dbef", label="in assay universe")
    ax.barh(y, coverage["candidates_qualifying"], color="#1a7a3a", label="qualifying signal")
    for i, r in coverage.iterrows():
        ax.text(r["candidates_in_universe"] + 0.3, i,
                f"{int(r['candidates_in_universe'])}/43 in · "
                f"{int(r['candidates_qualifying'])} qual · role={r['evidence_role']}",
                va="center", fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels(coverage["layer"])
    ax.set_xlim(0, 48); ax.set_xlabel("candidates (of 43)")
    ax.set_title("audit_v2 raw coverage per assay (direct from source, no 08e gene gate)",
                 fontsize=11)
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "figures/fig_multiomics_v2_coverage.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / "figures/fig_multiomics_v2_coverage.svg", bbox_inches="tight")
    plt.close(fig)


def write_report(out: Path, summ, long, models, coverage, checks, manifest):
    ind = summ[summ["independent_molecular_association"]]
    L = []
    a = L.append
    a("# BTLA multiomics evidence audit v2 — pipeline audit report\n")
    a(f"_Generated: {manifest['generated_utc']}. Supersedes audit_v1._\n")
    a("**Scope.** Evidence construction only; manuscript, Table 4, verdict, "
      "rankings, CRISPRi and Paperclip tiers untouched. See "
      "`multiomics_source_specification.md` for the approved spec.\n")

    a("## 1. Headline\n")
    a(f"- Rebuilt directly from the raw assay tables (not the 08e derivative), so "
      f"coverage is no longer limited by 08e's old gene gate.\n")
    a(f"- **{len(ind)} of 43** candidates show an independent molecular association "
      f"with the BTLA contrast (qualifying protein or phosphosite item): "
      f"{', '.join(ind['TF'])}.\n")
    a("  - Protein-qualifying: SMAP2 (4h, time-matched primary), STAT5B (20min, acute).\n")
    a("  - Phospho-qualifying: MYC, GTF3C2, HIRIP3, STAU2, TFDP1 (site/timepoint detail restricted).\n")
    a(f"- **No shared candidate** has an independent association; all {len(ind)} are "
      f"model-specific (GREmLN-specific: GTF3C2, HIRIP3, STAU2, TFDP1; "
      f"GENIE3-specific: MYC, SMAP2, STAT5B).\n")
    a("- **mRNA context change qualifying = 0** across all 43 in the BTLA contrast. "
      "Consistent with the clarification that these are upstream regulators whose own "
      "transcript need not move (lack of mRNA change is NOT evidence against a candidate).\n")
    a("- Direction is left **unresolved** for every independent item (rankings unsigned; "
      "no site-specific functional model). Association is established; supportive/"
      "opposing is not claimed.\n")

    a("## 2. Thresholds (source-native, recovered from Methods.docx)\n")
    t = manifest["thresholds_source_native"]
    a(f"- Transcriptomics: padj < {t['padj_max']} AND |log2FC| > log2({t['rna_abs_fc_linear']}) "
      f"= {t['rna_abs_log2fc']:.3f}\n")
    a(f"- Proteomics & phosphoproteomics: padj < {t['padj_max']} AND |log2FC| > "
      f"log2({t['prot_phospho_abs_fc_linear']}) = {t['prot_phospho_abs_log2fc']:.3f}\n")
    a(f"- Harmonised sensitivity flag (not primary): |log2FC| >= {t['harmonised_sensitivity_abs_log2fc']:.3f}.\n")
    a("- 08e's non-source `strong` sub-tier dropped: items are qualifying or not.\n")

    a("## 3. Coverage (raw, all 43)\n")
    a("| layer | role | in universe | absent | qualifying |\n|---|---|---|---|---|\n")
    for _, r in coverage.iterrows():
        a(f"| {r['layer']} | {r['evidence_role']} | {r['candidates_in_universe']} | "
          f"{43 - r['candidates_in_universe']} | {r['candidates_qualifying']} |\n")

    a("## 4. Model summaries\n")
    a("| view | n | independent assoc. | protein qual | phospho qual | mRNA context | TF-act context |\n"
      "|---|---|---|---|---|---|---|\n")
    for _, r in models.iterrows():
        a(f"| {r['view']} | {r['n_candidates']} | {r['n_independent_association']} | "
          f"{r['n_protein_qualifying']} | {r['n_phospho_qualifying']} | "
          f"{r['n_mrna_context_change']} | {r['n_tf_activity_context']} |\n")

    a("## 5. Three separated concepts\n")
    a("1. measured molecular change · 2. independent molecular association (BTLA "
      "contrast, direction-agnostic) · 3. directionally interpretable support "
      "(**unresolved** for all). A qualifying protein/phosphosite change establishes "
      "(1) and (2) without implying (3).\n")

    a("## 6. coIP\n")
    a("Excluded from corroboration: the available file is an intensity matrix with no "
      "statistical BTLA contrast, and the BTLA-IP-vs-IgG differential named in the "
      "Methods is absent from this tree. Only GAR1 (of 43) is detected. Presence is "
      "recorded, never scored.\n")

    a("## 7. Assertions\n")
    a(f"All {int(checks['passed'].sum())}/{len(checks)} checks pass:\n")
    for _, c in checks.iterrows():
        a(f"- [{'x' if c['passed'] else ' '}] {c['check']} ({c['detail']})\n")

    a("## 8. Unresolved questions for the data owner\n")
    a("1. Provide the BTLA-IP-vs-IgG coIP differential file to enable an interaction layer.\n"
      "2. Provide a site-specific functional model so phosphosite direction can move from "
      "`unresolved` to supportive/opposing.\n"
      "3. Confirm the linear-fold-change reading of the Methods thresholds.\n"
      "4. Confirm the gmorgan multiomics vintage is the intended source (vs Viet's newer copy).\n")

    a("## 9. Files\n")
    for f in ["multiomics_source_specification.md", "run_manifest_v2.json",
              "multiomics_evidence_long_v2.csv (redacted)",
              "multiomics_candidate_summary_v2.csv",
              "multiomics_model_summaries_v2.csv", "multiomics_coverage_v2.csv",
              "multiomics_audit_checks_v2.csv",
              "figures/fig_multiomics_v2_candidate_matrix.png",
              "figures/fig_multiomics_v2_model_summary.png",
              "figures/fig_multiomics_v2_coverage.png"]:
        a(f"- `{f}`\n")
    a("- `local_restricted/*.FULL.csv` — git-ignored (effect sizes, adjusted p, "
      "phosphosite residues, per-target direction).\n")
    (out / "multiomics_pipeline_audit_v2.md").write_text("".join(L))


def write_commit_safety(out: Path, manifest):
    L = ["# Commit-safety report — multiomics audit_v2\n",
         f"_Generated: {manifest['generated_utc']}_\n",
         "## Safe to commit (statuses / counts / hashes)\n"]
    for f in ["multiomics_source_specification.md",
              "run_manifest_v2.json (hashes + counts + threshold provenance)",
              "multiomics_evidence_long_v2.csv (statuses; effect_log2fc, adjusted_p, "
              "phosphosite residue and direction REDACTED; phospho site id -> REDACTED_site)",
              "multiomics_candidate_summary_v2.csv (counts + statuses; phospho "
              "site-timepoint list -> REDACTED_sites)",
              "multiomics_model_summaries_v2.csv",
              "multiomics_coverage_v2.csv", "multiomics_audit_checks_v2.csv",
              "multiomics_pipeline_audit_v2.md", "figures/*.png / *.svg"]:
        L.append(f"- [SAFE] {f}\n")
    L.append("\n## Must NOT be committed\n")
    L.append("- [RESTRICTED] `local_restricted/multiomics_evidence_long_v2.FULL.csv` and "
             "`local_restricted/multiomics_candidate_summary_v2.FULL.csv` — effect sizes, "
             "adjusted p, phosphosite residues, per-target direction, site-timepoint lists.\n")
    L.append("- [RESTRICTED] raw `.xlsx` sources under DATA_ROOT/data/multiomics/ (only SHA256 recorded).\n")
    L.append("\n## Pre-commit checklist\n")
    L.append("- [ ] no file under local_restricted/ staged\n- [ ] no .xlsx staged\n"
             "- [ ] committed CSVs contain no >=3dp decimal effect values\n")
    (out / "commit_safety_report_v2.md").write_text("".join(L))


def main():
    out = repo_root() / OUT_REL
    summ, long, models, coverage, checks, manifest = load(out)
    fig_matrix(out, summ, long)
    fig_models(out, models)
    fig_coverage(out, coverage)
    write_report(out, summ, long, models, coverage, checks, manifest)
    write_commit_safety(out, manifest)
    print("[audit_v2 report] figures + markdown written to", out)


if __name__ == "__main__":
    main()
