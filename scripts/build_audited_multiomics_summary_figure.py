#!/usr/bin/env python3
"""Section 6-7: two-panel audited BTLA multiomics summary figure.

Panel A: 43-candidate x 3-layer status matrix (mRNA contextual | protein,
phospho independent), grouped GREmLN-specific / shared / GENIE3-specific,
ordered by canonical model rank, qualifying cells annotated with timepoint(s).
Panel B: coverage-aware layer summary per group (qualifying / measured-no-signal
/ not-measured, with the measured denominator).

Source-specific primary thresholds only. No effect sizes, adjusted p or
phosphosite residues are read or drawn. Reads the committed audit_v2 tables.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT_REL = "results/multiomics/audit_v2"
LAYERS = [("transcriptomics", "Candidate\nmRNA"), ("proteomics", "Protein\nabundance"),
          ("phosphoproteomics", "Phosphosite\nchange")]
GROUPS = [("GREmLN_specific", "GREmLN-specific"), ("shared", "shared"),
          ("GENIE3_specific", "GENIE3-specific")]

STATUS_COLOR = {
    "qualifying association": "#1a7a3a",
    "measured, no qualifying signal": "#d9d9d9",
    "not measured": "#f7f7f7",
    "absent from assay universe": "#ececec",
    "ambiguous mapping": "#fddbc7",
    "failed quality control": "#fdae6b",
}
STATUS_SHORT = {
    "qualifying association": "qual",
    "measured, no qualifying signal": "ns",
    "not measured": "nm",
    "absent from assay universe": "absent",
    "ambiguous mapping": "amb",
    "failed quality control": "qc",
}


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def cell_status(long, tf, layer):
    sub = long[(long["TF"] == tf) & (long["layer"] == layer)]
    if not len(sub):
        return "not measured", []
    if (sub["qualifying_source"]).any():
        tps = sorted(set(sub[sub["qualifying_source"]]["timepoint"].astype(str)),
                     key=lambda t: ["2mn", "20mn", "1h", "4h", "24h"].index(t)
                     if t in ["2mn", "20mn", "1h", "4h", "24h"] else 99)
        return "qualifying association", tps
    st = set(sub["status"])
    if "ambiguous_mapping" in st:
        return "ambiguous mapping", []
    if "measured_no_qualifying_signal" in st:
        return "measured, no qualifying signal", []
    if "failed_qc" in st:
        return "failed quality control", []
    if (sub["status"] == "absent_from_assay_universe").all():
        return "absent from assay universe", []
    return "not measured", []


def main():
    repo = repo_root()
    out = repo / OUT_REL
    summ = pd.read_csv(out / "multiomics_candidate_summary_v2.csv")
    long = pd.read_csv(out / "multiomics_evidence_long_v2.csv")

    # ----- assemble figure data (ordered) -----
    rows = []
    ordered = []
    for gkey, _ in GROUPS:
        g = summ[summ["candidate_group"] == gkey].copy()
        rank_col = "gremln_dense_rank" if gkey != "GENIE3_specific" else "genie3_dense_rank"
        g = g.sort_values(rank_col)
        for _, r in g.iterrows():
            tf = r["TF"]
            ordered.append((tf, gkey))
            rec = {"TF": tf, "candidate_group": gkey,
                   "gremln_dense_rank": r["gremln_dense_rank"],
                   "genie3_dense_rank": r["genie3_dense_rank"]}
            for lkey, _ in LAYERS:
                st, tps = cell_status(long, tf, lkey)
                rec[f"{lkey}_status"] = st
                rec[f"{lkey}_qual_timepoints"] = ";".join(tps)
            rows.append(rec)
    figdata = pd.DataFrame(rows)
    figdata.to_csv(out / "audited_multiomics_figure_data.csv", index=False)

    # ----- Panel B counts (coverage-aware) -----
    panelB = []
    for gkey, glabel in GROUPS:
        sub = figdata[figdata["candidate_group"] == gkey]
        for lkey, llabel in LAYERS:
            s = sub[f"{lkey}_status"]
            qual = int((s == "qualifying association").sum())
            nosig = int((s == "measured, no qualifying signal").sum())
            failed = int((s == "failed quality control").sum())
            measured = qual + nosig + failed
            notmeas = int((s == "not measured").sum())
            absent = int((s == "absent from assay universe").sum())
            panelB.append({"group": gkey, "layer": lkey, "n_group": len(sub),
                           "qualifying": qual, "measured_no_signal": nosig + failed,
                           "measured_denominator": measured,
                           "not_measured_or_absent": notmeas + absent})
    panelB = pd.DataFrame(panelB)

    # ================= FIGURE =================
    n = len(ordered)
    fig = plt.figure(figsize=(14, 12.5))
    gsA = fig.add_axes([0.06, 0.10, 0.46, 0.78])   # Panel A
    gsB = fig.add_axes([0.61, 0.30, 0.36, 0.50])   # Panel B

    # ----- Panel A -----
    group_bounds = {}
    y = 0
    yticks, ylabels = [], []
    # draw from top: reverse so GREmLN-specific at top
    seq = []
    for gkey, glabel in GROUPS:
        members = [tf for tf, gg in ordered if gg == gkey]
        seq.append((gkey, glabel, members))
    total = sum(len(m) for _, _, m in seq)
    row_i = total
    for gkey, glabel, members in seq:
        top = row_i
        for tf in members:
            row_i -= 1
            yticks.append(row_i + 0.5)
            ylabels.append(tf)
            rec = figdata[figdata["TF"] == tf].iloc[0]
            for j, (lkey, _) in enumerate(LAYERS):
                st = rec[f"{lkey}_status"]
                # visual separation: gap between mRNA (col0) and protein/phospho (col1,2)
                xj = j + (0.18 if j >= 1 else 0)
                ax_fc = STATUS_COLOR[st]
                gsA.add_patch(mpatches.Rectangle((xj, row_i), 0.95, 0.92,
                              fc=ax_fc, ec="#8a8a8a", lw=0.5))
                if st == "qualifying association":
                    tps = rec[f"{lkey}_qual_timepoints"]
                    gsA.text(xj + 0.475, row_i + 0.46, tps, ha="center", va="center",
                             fontsize=6.2, color="white", fontweight="bold")
        group_bounds[gkey] = (row_i, top, glabel)

    gsA.set_xlim(-0.05, len(LAYERS) + 0.18 + 0.1)
    gsA.set_ylim(0, total)
    gsA.set_yticks(yticks)
    gsA.set_yticklabels(ylabels, fontsize=7.4)
    gsA.set_xticks([0.475, 1 + 0.18 + 0.475, 2 + 0.18 + 0.475])
    gsA.set_xticklabels([lab for _, lab in LAYERS], fontsize=9)
    gsA.xaxis.set_ticks_position("top")
    gsA.xaxis.set_label_position("top")
    for spine in gsA.spines.values():
        spine.set_visible(False)
    gsA.tick_params(length=0)
    # group separators + labels
    for gkey, (lo, hi, glabel) in group_bounds.items():
        gsA.plot([-0.05, len(LAYERS) + 0.28], [hi, hi], color="#333", lw=0.8)
        gsA.text(-0.55, (lo + hi) / 2, glabel, rotation=90, va="center", ha="center",
                 fontsize=9, fontweight="bold")
    # contextual vs independent separator
    gsA.axvline(1 + 0.09, color="#444", lw=1.4, ls="--")
    gsA.set_title("A  Candidate-level evidence matrix (BTLA+TCR vs TCR)\n"
                  "left of dashed line = transcript-derived contextual   |   "
                  "right = independent observed molecular",
                  fontsize=10.5, fontweight="bold", loc="left", pad=22)

    # ----- Panel B -----
    x = np.arange(len(LAYERS))
    width = 0.26
    offsets = {"GREmLN_specific": -width, "shared": 0, "GENIE3_specific": width}
    gcolor = {"GREmLN_specific": "#4477aa", "shared": "#999933", "GENIE3_specific": "#cc6677"}
    for gkey, glabel in GROUPS:
        sub = panelB[panelB["group"] == gkey].set_index("layer").loc[[l for l, _ in LAYERS]]
        qual = sub["qualifying"].values
        nosig = sub["measured_no_signal"].values
        notm = sub["not_measured_or_absent"].values
        xb = x + offsets[gkey]
        gsB.bar(xb, qual, width, color="#1a7a3a", edgecolor="w", label="_q")
        gsB.bar(xb, nosig, width, bottom=qual, color="#d9d9d9", edgecolor="w", label="_n")
        gsB.bar(xb, notm, width, bottom=qual + nosig, color="#f2f2f2",
                edgecolor="#bbb", hatch="//", label="_m")
        n_group = sub["n_group"].iloc[0]
        for xi, q, ns, nm in zip(xb, qual, nosig, notm):
            meas = q + ns
            gsB.text(xi, q + ns + nm + 0.3, f"{q}/{meas}", ha="center", va="bottom",
                     fontsize=6.6, color=gcolor[gkey], fontweight="bold")
    gsB.set_xticks(x)
    gsB.set_xticklabels([lab.replace("\n", " ") for _, lab in LAYERS], fontsize=9)
    gsB.set_ylabel("candidates")
    gsB.set_title("B  Coverage-aware layer summary\n"
                  "bars per group; label = qualifying / measured denominator",
                  fontsize=10.5, fontweight="bold", loc="left")
    gsB.spines[["top", "right"]].set_visible(False)
    gsB.set_ylim(0, 20)
    # legends placed below Panel B to avoid bars
    grp_handles = [mpatches.Patch(color=gcolor[g], label=l) for g, l in GROUPS]
    st_handles = [mpatches.Patch(fc="#1a7a3a", label="qualifying association"),
                  mpatches.Patch(fc="#d9d9d9", label="measured, no qualifying signal"),
                  mpatches.Patch(fc="#f2f2f2", ec="#bbb", hatch="//", label="not measured / absent")]
    leg1 = gsB.legend(handles=grp_handles, loc="upper left", fontsize=7.5,
                      frameon=False, title="group (left/centre/right bar)",
                      bbox_to_anchor=(0.0, -0.14))
    gsB.add_artist(leg1)
    gsB.legend(handles=st_handles, loc="upper left", fontsize=7.5, frameon=False,
               title="stack (bottom to top)", bbox_to_anchor=(0.52, -0.14))

    # global status legend for Panel A
    a_handles = [mpatches.Patch(fc=STATUS_COLOR[k], ec="#8a8a8a", label=k)
                 for k in ["qualifying association", "measured, no qualifying signal",
                           "not measured", "absent from assay universe",
                           "ambiguous mapping", "failed quality control"]]
    fig.legend(handles=a_handles, loc="lower center", ncol=3, fontsize=8,
               frameon=False, bbox_to_anchor=(0.29, 0.02))

    fig.suptitle("Audited candidate-level BTLA multiomics evidence (audit_v2, source-specific thresholds)",
                 fontsize=12.5, fontweight="bold", y=0.995)
    fig.savefig(out / "fig_audited_multiomics_summary.png", dpi=210, bbox_inches="tight")
    fig.savefig(out / "fig_audited_multiomics_summary.svg", bbox_inches="tight")
    plt.close(fig)

    # ----- caption -----
    ind = summ[summ["independent_molecular_association"]]["TF"].tolist()
    cov = {l: {"in_universe": int((figdata[f"{l}_status"] != "absent from assay universe").sum()),
               "qualifying": int((figdata[f"{l}_status"] == "qualifying association").sum())}
           for l, _ in LAYERS}
    caption = (
        "Figure | Audited candidate-level BTLA multiomics evidence.\n\n"
        "The matrix summarises candidate regulator measurements in the BTLA+TCR versus "
        "TCR contrast across transcriptomics, whole proteomics and phosphoproteomics. "
        "Candidate mRNA is shown as transcript-derived contextual evidence because the "
        "ranked regulators were nominated as upstream explanations of the BTLA "
        "transcriptional response rather than selected for changes in their own "
        "expression. Protein abundance and phosphosite changes are classified as "
        "independent observed molecular associations. Primary calls used the "
        "source-defined thresholds: adjusted p<0.05 and |log2 fold change|>log2(1.5) for "
        "transcriptomics, and adjusted p<0.05 and |log2 fold change|>log2(1.3) for protein "
        "and phosphosite measurements. Rows are grouped by model-specific or shared "
        "nomination. Grey cells indicate measured candidates without a qualifying signal; "
        "blank or separately marked cells indicate candidates not measured or outside the "
        "assay universe. Molecular direction was not interpreted as supportive or opposing "
        "because both model rankings were unsigned.\n\n"
        f"Coverage and qualifying-association summary. Of 43 candidates, "
        f"{cov['transcriptomics']['in_universe']} were measured for mRNA "
        f"(BTLA-specific qualifying = {cov['transcriptomics']['qualifying']}), "
        f"{cov['proteomics']['in_universe']} for protein "
        f"(qualifying = {cov['proteomics']['qualifying']}: SMAP2 at 4h, STAT5B at 20min), "
        f"and {cov['phosphoproteomics']['in_universe']} for phosphosites "
        f"(qualifying = {cov['phosphoproteomics']['qualifying']}: MYC, GTF3C2, HIRIP3, STAU2, TFDP1). "
        f"In total {len(ind)} of 43 candidates showed an independent molecular association "
        f"(protein or phosphosite): {', '.join(ind)}; none were shared candidates "
        f"(GREmLN-specific 4, shared 0, GENIE3-specific 3). No candidate showed a "
        f"BTLA-specific qualifying mRNA change; the ~10-15 mRNA signals in the earlier "
        f"08e chart arose from the TCR-vs-unstimulated contrast and a seed-inclusive "
        f"candidate set (see mrna_discrepancy_report.md). Under the stricter harmonised "
        f"threshold (|log2FC|>=0.585) the seven independent associations are unchanged."
    )
    (out / "fig_audited_multiomics_summary_caption.txt").write_text(caption)

    print("[figure] wrote fig_audited_multiomics_summary.{png,svg}, caption, data csv")
    print("coverage:", cov, "| independent:", ind)
    print(panelB.to_string(index=False))


if __name__ == "__main__":
    main()
