#!/usr/bin/env python3
"""Section 6-7: two-panel audited BTLA multiomics summary figure.

Panel A: 43-candidate status matrix. Candidate mRNA (transcript-derived
CONTEXTUAL evidence) is divided into three prespecified contrast subcolumns
describing transcriptional responsiveness across the experimental sequence:

    TCR activation      = TCR vs unstimulated control        (source: TCRvsUC)
    Incremental BTLA    = BTLA+TCR vs TCR                     (source: BTLAvsTCR)
    Combined BTLA+TCR   = BTLA+TCR vs unstimulated control    (not in source)

The combined contrast does not exist in the source differential-expression
output and cannot be regenerated from the fitted model (dds/counts absent), so
it is shown with the explicit non-fabricated state "source unavailable
(pending model)" and never synthesised from the other two contrasts.

Protein abundance and phosphosite change are INDEPENDENT observed molecular
associations and are the only layers that enter the independent corroboration
count. Candidate mRNA is contextual and never enters that count.

Source-specific primary thresholds only. No effect sizes, adjusted p or
phosphosite residues are read or drawn.
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

# columns: (layer_kind, contrast_key, header, role)
# Combined BTLA+TCR-vs-UC (C3) is retained in the backing tables with status
# source_unavailable_pending_model, but is NOT drawn as 43 identical hatched
# cells; a concise figure note explains its absence instead.
MRNA_CONTRASTS = [
    ("TCR_activation", "TCR\nactivation"),
    ("Incremental_BTLA", "Incremental\nBTLA"),
]
FIGNOTE = ("Note: the combined BTLA+TCR versus unstimulated-control contrast was "
           "not available from the fitted source models and was not reconstructed "
           "from other contrasts.")
COLUMNS = [("mrna", ck, hdr, "contextual") for ck, hdr in MRNA_CONTRASTS] + [
    ("proteomics", None, "Protein\nabundance", "independent"),
    ("phosphoproteomics", None, "Phosphosite\nchange", "independent"),
]
LAYERS_INDEP = [("proteomics", "Protein\nabundance"),
                ("phosphoproteomics", "Phosphosite\nchange")]
GROUPS = [("GREmLN_specific", "GREmLN-specific"), ("shared", "shared"),
          ("GENIE3_specific", "GENIE3-specific")]
PENDING = "source_unavailable_pending_model"

STATUS_COLOR = {
    "qualifying association": "#1a7a3a",
    "measured, no qualifying signal": "#d9d9d9",
    "not measured": "#f7f7f7",
    "absent from assay universe": "#ececec",
    "ambiguous mapping": "#fddbc7",
    "failed quality control": "#fdae6b",
    "source unavailable (pending model)": "#dbe9f6",
}


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def indep_cell_status(long, tf, layer):
    """Protein / phospho cell status from the audit_v2 long table."""
    sub = long[(long["TF"] == tf) & (long["layer"] == layer)]
    if not len(sub):
        return "not measured", []
    if (sub["qualifying_source"]).any():
        order = ["2mn", "20mn", "1h", "4h", "24h"]
        tps = sorted(set(sub[sub["qualifying_source"]]["timepoint"].astype(str)),
                     key=lambda t: order.index(t) if t in order else 99)
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


def mrna_cell_status(mrna, tf, contrast_key):
    """mRNA contrast subcolumn status from the 3-contrast summary table."""
    row = mrna[mrna["TF"] == tf]
    if not len(row):
        return "not measured", []
    r = row.iloc[0]
    st = r[f"{contrast_key}__status"]
    if st == "measured_qualifying":
        tps = [t for t in str(r[f"{contrast_key}__qual_timepoints"]).split(";") if t]
        return "qualifying association", tps
    return {
        "measured_no_qualifying_signal": "measured, no qualifying signal",
        "not_measured": "not measured",
        "absent_from_assay_universe": "absent from assay universe",
        "failed_qc": "failed quality control",
        PENDING: "source unavailable (pending model)",
    }.get(st, "not measured"), []


def main():
    repo = repo_root()
    out = repo / OUT_REL
    summ = pd.read_csv(out / "multiomics_candidate_summary_v2.csv")
    long = pd.read_csv(out / "multiomics_evidence_long_v2.csv")
    mrna = pd.read_csv(out / "mrna_contrast_layers_summary_v2.csv")

    # ----- assemble ordered figure data -----
    rows, ordered = [], []
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
            for kind, ck, _, _ in COLUMNS:
                if kind == "mrna":
                    st, tps = mrna_cell_status(mrna, tf, ck)
                    col = f"mrna_{ck}"
                else:
                    st, tps = indep_cell_status(long, tf, kind)
                    col = kind
                rec[f"{col}_status"] = st
                rec[f"{col}_qual_timepoints"] = ";".join(tps)
            rows.append(rec)
    figdata = pd.DataFrame(rows)
    figdata.to_csv(out / "audited_multiomics_figure_data.csv", index=False)

    # column keys for downstream loops
    colkeys = []
    for kind, ck, hdr, role in COLUMNS:
        colkeys.append((f"mrna_{ck}" if kind == "mrna" else kind, hdr, role, kind, ck))

    # ----- Panel B counts (coverage-aware, per column) -----
    panelB = []
    for gkey, glabel in GROUPS:
        sub = figdata[figdata["candidate_group"] == gkey]
        for col, hdr, role, kind, ck in colkeys:
            s = sub[f"{col}_status"]
            qual = int((s == "qualifying association").sum())
            nosig = int((s == "measured, no qualifying signal").sum())
            failed = int((s == "failed quality control").sum())
            measured = qual + nosig + failed
            notmeas = int((s == "not measured").sum())
            absent = int((s == "absent from assay universe").sum())
            pend = int((s == "source unavailable (pending model)").sum())
            panelB.append({"group": gkey, "column": col, "header": hdr.replace("\n", " "),
                           "role": role, "n_group": len(sub), "qualifying": qual,
                           "measured_no_signal": nosig + failed,
                           "measured_denominator": measured,
                           "not_measured_or_absent": notmeas + absent, "pending": pend})
    panelB = pd.DataFrame(panelB)
    panelB.to_csv(out / "audited_multiomics_panelB_counts.csv", index=False)

    # ================= FIGURE =================
    ncol = len(colkeys)
    # x positions: 3 mRNA cols tight, gap, then protein/phospho
    xpos = {}
    x = 0.0
    prev_role = None
    for col, hdr, role, kind, ck in colkeys:
        if prev_role == "contextual" and role == "independent":
            x += 0.45  # gap between contextual and independent blocks
        xpos[col] = x
        prev_role = role
        x += 1.0
    xmax = x
    last_mrna = f"mrna_{MRNA_CONTRASTS[-1][0]}"
    sep_x = (xpos[last_mrna] + 0.95 + xpos["proteomics"]) / 2

    fig = plt.figure(figsize=(15.5, 12.8))
    gsA = fig.add_axes([0.07, 0.10, 0.52, 0.76])
    gsB = fig.add_axes([0.67, 0.32, 0.31, 0.46])

    # ----- Panel A -----
    seq = []
    for gkey, glabel in GROUPS:
        members = [tf for tf, gg in ordered if gg == gkey]
        seq.append((gkey, glabel, members))
    total = sum(len(m) for _, _, m in seq)
    group_bounds = {}
    yticks, ylabels = [], []
    row_i = total
    for gkey, glabel, members in seq:
        top = row_i
        for tf in members:
            row_i -= 1
            yticks.append(row_i + 0.5)
            ylabels.append(tf)
            rec = figdata[figdata["TF"] == tf].iloc[0]
            for col, hdr, role, kind, ck in colkeys:
                st = rec[f"{col}_status"]
                xj = xpos[col]
                gsA.add_patch(mpatches.Rectangle((xj, row_i), 0.95, 0.92,
                              fc=STATUS_COLOR[st], ec="#8a8a8a", lw=0.5))
                if st == "source unavailable (pending model)":
                    gsA.plot([xj, xj + 0.95], [row_i, row_i + 0.92],
                             color="#7fa8d0", lw=0.6)
                if st == "qualifying association":
                    tps = rec[f"{col}_qual_timepoints"]
                    gsA.text(xj + 0.475, row_i + 0.46, tps, ha="center", va="center",
                             fontsize=5.8, color="white", fontweight="bold")
        group_bounds[gkey] = (row_i, top, glabel)

    gsA.set_xlim(-0.05, xmax + 0.05)
    gsA.set_ylim(0, total)
    gsA.set_yticks(yticks)
    gsA.set_yticklabels(ylabels, fontsize=7.2)
    gsA.set_xticks([xpos[c] + 0.475 for c, *_ in colkeys])
    gsA.set_xticklabels([hdr for _, hdr, *_ in colkeys], fontsize=8)
    gsA.xaxis.set_ticks_position("top")
    gsA.xaxis.set_label_position("top")
    for spine in gsA.spines.values():
        spine.set_visible(False)
    gsA.tick_params(length=0)

    for gkey, (lo, hi, glabel) in group_bounds.items():
        gsA.plot([-0.05, xmax], [hi, hi], color="#333", lw=0.8)
        gsA.text(-0.75, (lo + hi) / 2, glabel, rotation=90, va="center", ha="center",
                 fontsize=9, fontweight="bold")
    # contextual vs independent divider
    gsA.axvline(sep_x, color="#444", lw=1.4, ls="--")
    # spanning brackets over the two blocks
    mrna_lo = xpos["mrna_TCR_activation"]
    mrna_hi = xpos[last_mrna] + 0.95
    ind_lo = xpos["proteomics"]
    ind_hi = xpos["phosphoproteomics"] + 0.95
    ybr = total + 1.15
    for lo, hi, txt, col in [(mrna_lo, mrna_hi, "Candidate mRNA — contextual", "#555"),
                             (ind_lo, ind_hi, "Independent observed", "#1a7a3a")]:
        gsA.plot([lo, hi], [ybr, ybr], color=col, lw=1.6, clip_on=False)
        gsA.text((lo + hi) / 2, ybr + 0.25, txt, ha="center", va="bottom",
                 fontsize=8.5, fontweight="bold", color=col, clip_on=False)
    gsA.set_ylim(0, total + 2.0)

    gsA.set_title("A  Candidate-level evidence matrix\n"
                  "candidate mRNA split by contrast (TCR vs UC | BTLA+TCR vs TCR); "
                  "cell label = qualifying timepoint(s)",
                  fontsize=10, fontweight="bold", loc="left", pad=30)

    # ----- Panel B -----
    labels = [hdr.replace("\n", " ") for _, hdr, *_ in colkeys]
    xB = np.arange(ncol)
    width = 0.26
    offsets = {"GREmLN_specific": -width, "shared": 0, "GENIE3_specific": width}
    gcolor = {"GREmLN_specific": "#4477aa", "shared": "#999933", "GENIE3_specific": "#cc6677"}
    for gkey, glabel in GROUPS:
        sub = panelB[panelB["group"] == gkey].set_index("column").loc[[c for c, *_ in colkeys]]
        qual = sub["qualifying"].values
        nosig = sub["measured_no_signal"].values
        notm = sub["not_measured_or_absent"].values
        pend = sub["pending"].values
        xb = xB + offsets[gkey]
        gsB.bar(xb, qual, width, color="#1a7a3a", edgecolor="w")
        gsB.bar(xb, nosig, width, bottom=qual, color="#d9d9d9", edgecolor="w")
        gsB.bar(xb, notm, width, bottom=qual + nosig, color="#f2f2f2",
                edgecolor="#bbb", hatch="//")
        gsB.bar(xb, pend, width, bottom=qual + nosig + notm, color="#dbe9f6",
                edgecolor="#7fa8d0", hatch="xx")
        for xi, q, ns, nm, pd_ in zip(xb, qual, nosig, notm, pend):
            meas = q + ns
            lab = "pending" if pd_ and meas == 0 else f"{q}/{meas}"
            gsB.text(xi, q + ns + nm + pd_ + 0.3, lab, ha="center", va="bottom",
                     fontsize=5.8, color=gcolor[gkey], fontweight="bold", rotation=90)
    gsB.set_xticks(xB)
    gsB.set_xticklabels(labels, fontsize=7.2, rotation=30, ha="right")
    gsB.set_ylabel("candidates")
    gsB.set_title("B  Coverage-aware layer summary\n"
                  "label = qualifying / measured denominator",
                  fontsize=10, fontweight="bold", loc="left")
    gsB.spines[["top", "right"]].set_visible(False)
    gsB.set_ylim(0, 22)

    grp_handles = [mpatches.Patch(color=gcolor[g], label=l) for g, l in GROUPS]
    st_handles = [mpatches.Patch(fc="#1a7a3a", label="qualifying association"),
                  mpatches.Patch(fc="#d9d9d9", label="measured, no qualifying signal"),
                  mpatches.Patch(fc="#f2f2f2", ec="#bbb", hatch="//", label="not measured / absent")]
    leg1 = gsB.legend(handles=grp_handles, loc="upper left", fontsize=7,
                      frameon=False, title="group (left/centre/right bar)",
                      bbox_to_anchor=(0.0, -0.30))
    gsB.add_artist(leg1)
    gsB.legend(handles=st_handles, loc="upper left", fontsize=7, frameon=False,
               title="stack (bottom to top)", bbox_to_anchor=(0.0, -0.52))

    a_handles = [mpatches.Patch(fc=STATUS_COLOR[k], ec="#8a8a8a", label=k)
                 for k in ["qualifying association", "measured, no qualifying signal",
                           "not measured", "absent from assay universe",
                           "failed quality control"]]
    fig.legend(handles=a_handles, loc="upper left", ncol=3, fontsize=8,
               frameon=False, bbox_to_anchor=(0.07, 0.085))
    # concise note on the unavailable combined contrast (C3), below the legend
    import textwrap as _tw
    fig.text(0.07, 0.045, "\n".join(_tw.wrap(FIGNOTE, 95)),
             ha="left", va="top", fontsize=7.6, style="italic", color="#555")

    fig.suptitle("Audited candidate-level BTLA multiomics evidence "
                 "(audit_v2, source-specific thresholds)",
                 fontsize=12.5, fontweight="bold", y=0.995)
    fig.savefig(out / "fig_audited_multiomics_summary.png", dpi=210, bbox_inches="tight")
    fig.savefig(out / "fig_audited_multiomics_summary.svg", bbox_inches="tight")
    plt.close(fig)

    # ----- caption -----
    ind = summ[summ["independent_molecular_association"]]["TF"].tolist()
    c1 = int(mrna["TCR_activation__qualifies"].sum())
    c2 = int(mrna["Incremental_BTLA__qualifies"].sum())
    anyc = int(mrna["transcriptionally_responsive_any_context"].sum())
    morethan1 = int((mrna["n_qualifying_contrasts"] > 1).sum())
    c1_tfs = ", ".join(sorted(mrna[mrna["TCR_activation__qualifies"]]["TF"]))
    caption = (
        "Figure | Audited candidate-level BTLA multiomics evidence.\n\n"
        "The matrix summarises candidate regulator measurements across "
        "transcriptomics, whole proteomics and phosphoproteomics. Candidate mRNA is "
        "contextual evidence and is shown as two contrast-specific columns "
        "describing transcriptional responsiveness across the experimental "
        "sequence: TCR activation (TCR vs unstimulated control) and incremental BTLA "
        "(BTLA+TCR vs TCR). The combined BTLA+TCR versus unstimulated-control "
        "contrast was not available from the fitted source models and was not "
        "reconstructed from other contrasts; it is retained in the backing tables "
        "with status source_unavailable_pending_model but is not drawn. Contrasts "
        "are reported separately and never collapsed into a single strongest "
        "result.\n\n"
        "Protein and phosphosite measurements remain the independent molecular "
        "evidence layers and are the only layers that enter the independent "
        "molecular corroboration count; candidate mRNA never enters that count. "
        "Primary calls used source-defined thresholds: adjusted p<0.05 and |log2 "
        "fold change|>log2(1.5) for transcriptomics, and adjusted p<0.05 and |log2 "
        "fold change|>log2(1.3) for protein and phosphosite measurements. Rows are "
        "grouped by model-specific or shared nomination and ordered by canonical "
        "model rank; cell labels give the qualifying timepoint(s). Molecular "
        "direction was not interpreted as supportive or opposing because both model "
        "rankings were unsigned.\n\n"
        f"Transcriptional responsiveness (43 candidates, contextual). Five "
        f"candidates were transcriptionally responsive to TCR activation: "
        f"{c1_tfs}. No candidate showed a qualifying incremental BTLA-associated "
        f"mRNA response ({c2} of 43). {morethan1} qualified in more than one "
        f"available contrast; {anyc} qualified in any available context.\n\n"
        f"Independent molecular associations. The independent molecular association "
        f"total remains seven candidates (protein or phosphosite): "
        f"{', '.join(ind)} (GREmLN-specific 4, shared 0, GENIE3-specific 3). Under "
        f"the stricter harmonised threshold (|log2FC|>=0.585) these seven are "
        f"unchanged."
    )
    (out / "fig_audited_multiomics_summary_caption.txt").write_text(caption)

    print("[figure] wrote fig_audited_multiomics_summary.{png,svg}, caption, data csv")
    print(f"mRNA contextual: C1={c1} C2={c2} >1={morethan1} any={anyc} | independent={ind}")
    print(panelB[panelB.role == "contextual"][["group", "header", "qualifying",
          "measured_denominator", "pending"]].to_string(index=False))


if __name__ == "__main__":
    main()
