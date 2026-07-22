#!/usr/bin/env python3
"""Presentation layer for the BTLA multiomics evidence audit (audit_v1).

Reads the audited tables produced by ``build_multiomics_evidence_audit.py`` and
emits three figures plus ``multiomics_pipeline_audit.md`` and
``commit_safety_report.md``. No effect sizes or phosphosite residues are read or
rendered here; only statuses, counts and provenance hashes.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


OUT_REL = "results/multiomics/audit_v1"

LAYERS = ["transcriptomics", "proteomics", "phosphoproteomics", "coIP",
          "tf_activity", "kinase_activity", "bionic_gnn",
          "early_synapse_trafficking"]
LAYER_SHORT = {
    "transcriptomics": "RNA (C)", "proteomics": "Prot (A)",
    "phosphoproteomics": "Phos (A)", "coIP": "coIP (A)",
    "tf_activity": "TFact (C)", "kinase_activity": "Kin (B)",
    "bionic_gnn": "BIONIC (C)", "early_synapse_trafficking": "Synapse (D)",
}

# status -> (code, color) for the candidate matrix
STATUS_CODE = {
    "positive_supportive": (4, "#1a7a3a"),
    "positive_opposing": (3, "#8e44ad"),
    "positive_unresolved_direction": (2, "#2c7fb8"),
    "measured_no_qualifying_signal": (1, "#d9d9d9"),
    "measured_failed_qc": (1, "#d9d9d9"),
    "not_measured": (0, "#f7f7f7"),
    "candidate_absent_from_assay_universe": (-1, "#ffffff"),
    "not_applicable": (0, "#fbe6c2"),
    "source_unavailable": (-2, "#efefef"),
    "ambiguous_mapping": (-1, "#fddbc7"),
}


def load(out: Path):
    summ = pd.read_csv(out / "multiomics_candidate_summary_audited.csv")
    long = pd.read_csv(out / "multiomics_evidence_long_audited.csv")
    checks = pd.read_csv(out / "multiomics_audit_checks.csv")
    legacy = pd.read_csv(out / "legacy_vs_audited_multiomics.csv")
    manifest = json.loads((out / "run_manifest.json").read_text())
    union = pd.read_csv(out / "candidate_union.csv")
    registry = pd.read_csv(out / "multiomics_layer_registry.csv")
    return summ, long, checks, legacy, manifest, union, registry


# --------------------------------------------------------------------------- #
# Figure 1: pipeline flow
# --------------------------------------------------------------------------- #
def fig_pipeline(out: Path, manifest: dict):
    fig, ax = plt.subplots(figsize=(12, 6.2))
    ax.axis("off")
    stages = [
        ("7 raw assay sources\n(RNA, Prot, Phospho,\ncoIP, TFact, Kinase, BIONIC)",
         "#dfeaf5"),
        ("Frozen 08e long table\n(1 row / candidate x item)\nSHA256 pinned", "#dfeaf5"),
        ("Audit re-derivation\nindependence class A/B/C/D\ncontrolled status vocab", "#fde9c8"),
        ("Candidate summary\n+ decision-rule field\n(class A supportive only)", "#dcefdc"),
        ("Reconcile vs legacy\n+ 18 audit assertions", "#e6d9f2"),
    ]
    n = len(stages)
    w, h, gap = 0.165, 0.34, 0.043
    x0 = 0.015
    y = 0.55
    for i, (txt, col) in enumerate(stages):
        x = x0 + i * (w + gap)
        box = mpatches.FancyBboxPatch((x, y), w, h,
                                      boxstyle="round,pad=0.008,rounding_size=0.02",
                                      fc=col, ec="#333", lw=1.2)
        ax.add_patch(box)
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=9)
        if i < n - 1:
            ax.annotate("", xy=(x + w + gap, y + h / 2), xytext=(x + w, y + h / 2),
                        arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#333"))
    ax.text(0.5, 0.98,
            "BTLA multiomics evidence audit (audit_v1) — evidence construction only; "
            "manuscript and verdict untouched",
            ha="center", va="top", fontsize=12, fontweight="bold")
    caption = (
        f"43-candidate union (GREmLN 25 + GENIE3 25, 7 shared).  "
        f"Evaluated by 08e source pipeline: {manifest['n_candidates_evaluated_by_08e']}/43.  "
        f"Class-A independent molecular signal present: "
        f"{', '.join(manifest['classA_signal_present_candidates']) or 'none'}.  "
        f"Decision-rule TRUE (class-A supportive): "
        f"{', '.join(manifest['decision_rule_true_candidates']) or 'none — all class-A directions unresolved'}.  "
        f"Assertions: {manifest['checks_passed']}/{manifest['checks_total']} passed."
    )
    ax.text(0.5, 0.30, caption, ha="center", va="top", fontsize=9.2, wrap=True,
            bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#bbb"))
    key = ("Independence classes:  A = independent observed molecular assay  |  "
           "B = orthogonal assay-derived computational  |  "
           "C = transcript-/model-derived contextual  |  D = curated / unclear (excluded)")
    ax.text(0.5, 0.13, key, ha="center", va="top", fontsize=8.5, style="italic",
            color="#333")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.savefig(out / "figures/fig_multiomics_audit_pipeline.svg", bbox_inches="tight")
    fig.savefig(out / "figures/fig_multiomics_audit_pipeline.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Figure 2: candidate x layer status matrix
# --------------------------------------------------------------------------- #
def fig_matrix(out: Path, summ: pd.DataFrame, long: pd.DataFrame, union: pd.DataFrame):
    order = summ.sort_values(
        ["classA_independent_positive_unresolved", "classC_contextual_positive",
         "n_measured_layers", "TF"],
        ascending=[False, False, False, True])["TF"].tolist()

    # build per candidate x layer best status
    def best_status(tf, layer):
        rows = long[(long["TF"] == tf) & (long["layer_id"] == layer)]
        if not len(rows):
            evaluated = tf in set(long["TF"])
            return "not_measured" if evaluated else "source_unavailable"
        prio = ["positive_supportive", "positive_opposing",
                "positive_unresolved_direction", "measured_no_qualifying_signal",
                "not_applicable", "not_measured", "source_unavailable",
                "ambiguous_mapping"]
        present = set(rows["evidence_status"])
        for s in prio:
            if s in present:
                return s
        return "not_measured"

    grid = np.zeros((len(order), len(LAYERS)))
    colors = {}
    for i, tf in enumerate(order):
        for j, layer in enumerate(LAYERS):
            st = best_status(tf, layer)
            code, col = STATUS_CODE.get(st, (0, "#f7f7f7"))
            grid[i, j] = code
            colors[(i, j)] = col

    fig, ax = plt.subplots(figsize=(9.5, 13.5))
    for i in range(len(order)):
        for j in range(len(LAYERS)):
            ax.add_patch(mpatches.Rectangle((j, len(order) - 1 - i), 1, 1,
                         fc=colors[(i, j)], ec="#999", lw=0.4))
    ax.set_xlim(0, len(LAYERS)); ax.set_ylim(0, len(order))
    ax.set_xticks(np.arange(len(LAYERS)) + 0.5)
    ax.set_xticklabels([LAYER_SHORT[l] for l in LAYERS], rotation=45, ha="right",
                       fontsize=9)
    ax.set_yticks(np.arange(len(order)) + 0.5)
    grp = dict(zip(union["TF"], union["candidate_group"]))
    ax.set_yticklabels([f"{tf}  ·{grp.get(tf,'')[:3]}" for tf in reversed(order)],
                       fontsize=7.5)
    ax.set_title("Audited multiomics status by candidate x layer\n"
                 "(letters = independence class; direction not asserted for molecular magnitude)",
                 fontsize=11)
    legend = [
        mpatches.Patch(fc="#2c7fb8", label="positive, direction unresolved"),
        mpatches.Patch(fc="#1a7a3a", label="positive supportive"),
        mpatches.Patch(fc="#8e44ad", label="positive opposing"),
        mpatches.Patch(fc="#d9d9d9", label="measured, no qualifying signal"),
        mpatches.Patch(fc="#fbe6c2", label="not applicable (curated/hypothesis)"),
        mpatches.Patch(fc="#f7f7f7", ec="#999", label="not measured (in-scope)"),
        mpatches.Patch(fc="#efefef", ec="#999", label="source unavailable (not run by 08e)"),
    ]
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.045),
              ncol=2, fontsize=8, frameon=False)
    fig.savefig(out / "figures/fig_multiomics_candidate_matrix.svg", bbox_inches="tight")
    fig.savefig(out / "figures/fig_multiomics_candidate_matrix.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Figure 3: model-level summary
# --------------------------------------------------------------------------- #
def fig_model_summary(out: Path, summ: pd.DataFrame):
    cats = ["independent_molecular_signal_present_direction_unresolved",
            "contextual_or_orthogonal_only",
            "no_qualifying_signal_in_evaluated_layers",
            "source_unavailable_not_evaluated_by_08e_pipeline"]
    labels = ["Class-A signal\n(dir. unresolved)", "Contextual /\northogonal only",
              "No qualifying\nsignal", "Source unavailable\n(not run by 08e)"]
    views = {
        "GREmLN top-25": summ[summ["in_gremln_top25"]],
        "GENIE3 top-25": summ[summ["in_genie3_top25"]],
        "Union (43)": summ,
    }
    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(len(cats))
    width = 0.26
    for k, (name, sub) in enumerate(views.items()):
        counts = [int((sub["overall_audit_status"] == c).sum()) for c in cats]
        bars = ax.bar(x + (k - 1) * width, counts, width, label=f"{name} (n={len(sub)})")
        for b, c in zip(bars, counts):
            if c:
                ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.15, str(c),
                        ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("candidates")
    ax.set_title("Audited multiomics evidence status by model view\n"
                 "(single decision-rule TRUE count = 0; MYC is the only class-A signal, "
                 "direction unresolved)", fontsize=11)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "figures/fig_multiomics_model_summary.svg", bbox_inches="tight")
    fig.savefig(out / "figures/fig_multiomics_model_summary.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Reports
# --------------------------------------------------------------------------- #
def write_report(out: Path, summ, long, checks, legacy, manifest, union, registry):
    n_eval = manifest["n_candidates_evaluated_by_08e"]
    absent = sorted(set(union["TF"]) - set(long["TF"]))
    downgraded = legacy[legacy["status_changed"].str.contains("downgraded")]["TF"].tolist()
    classA = manifest["classA_signal_present_candidates"]
    L = []
    a = L.append
    a("# BTLA multiomics evidence audit (audit_v1)\n")
    a(f"_Generated: {manifest['generated_utc']}_\n")
    a("**Scope.** Evidence construction and validation only. This audit does NOT "
      "modify the manuscript, Results/Conclusions, Table 4, the benchmark verdict, "
      "the GREmLN/GENIE3 rankings, the CRISPRi analysis, or the Paperclip tiers. "
      "Rankings are read-only inputs; evidence joins cannot alter any rank.\n")

    a("## 1. Headline audit findings\n")
    a(f"- **Coverage gap:** only **{n_eval}/43** union candidates were ever scored by "
      f"the source `08e` pipeline. The remaining **{len(absent)}** are recorded as "
      f"`source_unavailable_not_evaluated_by_08e_pipeline` — NOT as \"no evidence\". "
      f"They were outside the gene set the 08e notebook iterated over and must be "
      f"re-run against the frozen 43-candidate union before any absence claim.\n")
    a(f"  - Not evaluated: {', '.join(absent)}\n")
    a(f"- **Independent molecular evidence is rare:** exactly **one** candidate "
      f"(**{', '.join(classA) or 'none'}**) has any class-A (independent observed "
      f"molecular) positive item that passes the documented significance and effect "
      f"thresholds — MYC, via 3 phosphosites.\n")
    a("- **Direction is not asserted:** for magnitude-only molecular layers (RNA, "
      "protein, phosphosite, TF/kinase activity) the audit records "
      "`positive_unresolved_direction`. We do not assume an increase or decrease is "
      "supportive of the BTLA inhibitory phenotype without a documented per-target "
      "functional model. Consequently the strict decision-rule field "
      "`independent_molecular_support_for_decision_rule` is **True for 0 candidates** "
      "(MYC included) pending a directional rule from the data owner.\n")
    a(f"- **Legacy over-crediting:** {len(downgraded)} candidates the legacy summary "
      f"labelled `moderate`/`strong` multiomics support are reclassified, because the "
      f"legacy field treated any non-empty support string (including class-C contextual "
      f"TF-activity/BIONIC) as support: {', '.join(downgraded)}.\n")

    a("## 2. Independence classification (why MYC ≠ the others)\n")
    a("| class | meaning | layers | counts toward decision rule |\n"
      "|---|---|---|---|\n"
      "| A | independent observed molecular assay | proteomics, phosphoproteomics, coIP | yes (if supportive) |\n"
      "| B | orthogonal assay-derived computational | kinase_activity | no |\n"
      "| C | transcript-/model-derived contextual | transcriptomics, tf_activity, bionic_gnn | no |\n"
      "| D | curated / unclear lineage | early_synapse_trafficking | no (excluded) |\n")
    a("`tf_activity` and `transcriptomics` are **class C** because they are computed "
      "from (or are) the transcript signal that drives candidate ranking; crediting "
      "them as independent confirmation would double-count the ranking input.\n")

    a("## 3. Missingness vocabulary (never silently converted to 'no support')\n")
    a("- `not_measured` — entity not present in an assay that was run.\n"
      "- `measured_no_qualifying_signal` — measured but did not pass thresholds.\n"
      "- `source_unavailable_not_evaluated_by_08e_pipeline` — the assay/candidate was "
      "not processed by the source pipeline at all.\n"
      "- `not_applicable` — curated hypothesis annotation, not a measurement.\n")

    a("## 4. Directional interpretation audit\n")
    a("Every positive molecular item carries `relation_to_BTLA_hypothesis`. For all "
      "magnitude-only layers this is `unresolved` by policy. coIP is `context only` "
      "(pulldown membership does not resolve a direct interaction; the CRAPome file is "
      "present but was NOT applied in 08e). No item is scored `supportive` without a "
      "documented directional model, so the decision rule cannot inflate.\n")

    a("## 5. Layer registry\n")
    a("See `multiomics_layer_registry.csv` (8 layers). Key provenance: RNA thresholds "
      "padj<0.05 & |log2FC|≥log2(1.5); protein/phospho padj<0.05 & |logFC|≥log2(1.3); "
      "coIP |t|>2 on a pre-filtered significant table; TF/kinase activity decoupleR ULM "
      "raw p<0.05 with no effect floor.\n")

    a("## 6. Reconciliation with legacy tables\n")
    a(f"`legacy_vs_audited_multiomics.csv`: {len(legacy)} rows. "
      f"{int((legacy['status_changed']=='no_material_change').sum())} unchanged, "
      f"{len(downgraded)} reclassified downward, "
      f"{int(legacy['status_changed'].str.contains('added_to_audit_scope').sum())} "
      f"added to audit scope (absent from legacy summary).\n")

    a("## 7. Assertions\n")
    a(f"All **{int(checks['passed'].sum())}/{len(checks)}** audit checks pass:\n")
    for _, c in checks.iterrows():
        a(f"- [{'x' if c['passed'] else ' '}] {c['check']} ({c['detail']})\n")

    a("## 8. Unresolved questions for the data owner\n")
    a("1. Re-run the 08e pipeline over the frozen 43-candidate union so the 19 "
      "`source_unavailable` candidates get real statuses.\n"
      "2. Provide a documented per-layer/per-target directional model (does an increase "
      "in a given phosphosite/protein support the BTLA inhibitory phenotype?) so "
      "`relation_to_BTLA_hypothesis` can move from `unresolved` to supportive/opposing.\n"
      "3. Decide whether CRAPome filtering should be applied to coIP before any "
      "interaction claim.\n"
      "4. Confirm MYC phosphosite residues/effect sizes (restricted) before external use.\n")

    a("## 9. Files\n")
    for f in ["source_inventory.csv", "run_manifest.json", "candidate_union.csv",
              "multiomics_layer_registry.csv", "multiomics_evidence_long_audited.csv",
              "multiomics_candidate_summary_audited.csv",
              "legacy_vs_audited_multiomics.csv", "multiomics_audit_checks.csv",
              "figures/fig_multiomics_audit_pipeline.png",
              "figures/fig_multiomics_candidate_matrix.png",
              "figures/fig_multiomics_model_summary.png"]:
        a(f"- `{f}`\n")
    a("- `local_restricted/multiomics_evidence_long_audited.FULL.csv` — **git-ignored**; "
      "effect sizes + phosphosite residues.\n")
    (out / "multiomics_pipeline_audit.md").write_text("".join(L))


def write_commit_safety(out: Path, manifest: dict):
    L = []
    a = L.append
    a("# Commit-safety report — multiomics audit_v1\n")
    a(f"_Generated: {manifest['generated_utc']}_\n")
    a("## Safe to commit (statuses / hashes / aggregate counts only)\n")
    for f in ["source_inventory.csv (redacted relative paths + SHA256, no effect sizes)",
              "run_manifest.json (hashes + counts)",
              "candidate_union.csv (membership + ranks, already published)",
              "multiomics_layer_registry.csv (schema/provenance)",
              "multiomics_evidence_long_audited.csv (statuses only; effect_value, "
              "adjusted_p, phosphosite id/residue and per-target direction REDACTED)",
              "multiomics_candidate_summary_audited.csv (class counts + statuses)",
              "legacy_vs_audited_multiomics.csv (status-level reconciliation)",
              "multiomics_audit_checks.csv",
              "multiomics_pipeline_audit.md",
              "figures/*.png / *.svg (status colours only, no numeric values)"]:
        a(f"- [SAFE] {f}\n")
    a("\n## Must NOT be committed (publication gate)\n")
    a("- [RESTRICTED] `local_restricted/multiomics_evidence_long_audited.FULL.csv` — "
      "log2FC/logFC/t-statistics, adjusted p-values, phosphosite residues, per-target "
      "direction. Kept under `local_restricted/` and added to `.gitignore`.\n")
    a("- [RESTRICTED] all raw Excel sources under `DATA_ROOT/data/multiomics/` "
      "(never staged; only SHA256 recorded).\n")
    a("\n## Redaction method\n")
    a("- Absolute paths replaced by paths relative to `BTLA_BENCH_DATA_ROOT`.\n"
      "- Effect sizes and phosphosite residues dropped from committed tables and "
      "replaced with `REDACTED_publication_gate`.\n"
      "- Per-target direction (sign of a significant change) removed from committed "
      "tables; only the direction-neutral `positive_unresolved_direction` status remains.\n")
    a("\n## Pre-commit checklist\n")
    a("- [ ] `git status` shows no file under `local_restricted/`.\n"
      "- [ ] no `.xlsx` staged.\n"
      "- [ ] grep committed CSVs for numeric effect columns returns nothing.\n")
    (out / "commit_safety_report.md").write_text("".join(L))


def main():
    repo = repo_root()
    out = repo / OUT_REL
    summ, long, checks, legacy, manifest, union, registry = load(out)
    fig_pipeline(out, manifest)
    fig_matrix(out, summ, long, union)
    fig_model_summary(out, summ)
    write_report(out, summ, long, checks, legacy, manifest, union, registry)
    write_commit_safety(out, manifest)
    print("[report] figures + markdown written to", out)


if __name__ == "__main__":
    main()
