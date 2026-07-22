#!/usr/bin/env python3
"""Paperclip-style multiomics evidence tables (GREmLN + GENIE3).

Two matched publication figures, each covering one model's seed-excluded top 25,
styled like the Paperclip candidate-evidence summaries.

Ordering is taken EXACTLY from Table 4A / 4B (canonical ordinal ranking).
Multiomics evidence is joined by TF symbol and never used to re-sort rows.

Visual encoding distinguishes evidence class:
  - contextual candidate mRNA  -> blue pills
  - independent protein/phosphosite -> green pills
Written labels remain self-explanatory without colour.

Built only from audited Table S7 / audit_v2; evidence is not rerun.
"""
from __future__ import annotations

import re
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd

DASH = "\u2014"

# Distinct encodings for the two evidence classes
BLUE_CTX = "#2c7bb6"       # contextual mRNA (qualifying)
GREEN_IND = "#1b7837"      # independent protein / phosphosite (qualifying)
GREY_NONE = "#d9d9d9"
YELLOW_WEAK = "#fee08b"
ZEBRA = "#f7f7f7"
HEADER_BG = "#f2f2f2"
TEXT_DARK = "#222222"

OUT_STEM = {
    "GREmLN": "fig_gremln_multiomics_evidence",
    "GENIE3": "fig_genie3_multiomics_evidence",
}
TABLE4 = {
    "GREmLN": "results/publication_assets/tables/table4a_gremln_top25_evidence.csv",
    "GENIE3": "results/publication_assets/tables/table4b_genie3_top25_evidence.csv",
}

# Expected reconciliation (assertions, not hardcoded replacements)
EXP = {
    "GREmLN": {"mrna": 2, "incr": 0, "indep": 4},
    "GENIE3": {"mrna": 5, "incr": 0, "indep": 3},
}


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def strip_shared_mark(name: str) -> str:
    """Remove Table 4 shared superscript (ˢ, U+02E2) from the candidate symbol."""
    return re.sub(r"\u02e2$", "", str(name)).strip()


def shorten_tps(s: str) -> str:
    s = str(s)
    return (s.replace(" min", "min").replace(" h", "h")
             .replace("Increased: ", "↑ ").replace("Decreased: ", "↓ ")
             .replace("Qualifying association: ", "")
             .replace(" qualifying sites", " sites")
             .replace(" qualifying site", " site"))


def classify_cell(text: str, kind: str):
    """Return (display_label, style).

    style in {ctx, indep, grey, dash, qc}:
      ctx   = qualifying contextual mRNA (blue)
      indep = qualifying independent protein/phosphosite (green)
      grey  = measured, no qualifying change
      dash  = absent / not assessable
      qc    = failed QC
    """
    t = str(text)
    if kind in ("tcr", "incr"):
        if t.startswith("Increased") or t.startswith("Decreased"):
            return shorten_tps(t), "ctx"
        if t == "Measured, no qualifying change":
            return "No qualifying change", "grey"
        if t == "Failed quality control":
            return "Failed QC", "qc"
        if t in ("Absent from assay universe",):
            return DASH, "dash"
        return "No qualifying change", "grey"

    if kind in ("protein", "phospho"):
        if t.startswith("Qualifying"):
            return shorten_tps(t), "indep"
        if t == "Measured, no qualifying change":
            return "No qualifying change", "grey"
        if t == "Failed quality control":
            return "Failed QC", "qc"
        return DASH, "dash"

    # independent association summary
    if t.startswith("Yes"):
        if "protein and phosphosite" in t:
            return "Protein + phospho", "indep"
        if "protein" in t:
            return "Protein", "indep"
        if "phosphosite" in t:
            return "Phosphosite", "indep"
        return "Yes", "indep"
    if t == "No qualifying association":
        return "No qualifying change", "grey"
    return DASH, "dash"


def pill(ax, cx, cy, text, style, w=0.78, h=0.50, fs=6.6):
    if style == "dash":
        ax.text(cx, cy, DASH, ha="center", va="center", fontsize=10, color="#666")
        return
    if style == "ctx":
        fc, tc = BLUE_CTX, "white"
    elif style == "indep":
        fc, tc = GREEN_IND, "white"
    elif style == "qc":
        fc, tc = YELLOW_WEAK, "black"
    else:
        fc, tc = GREY_NONE, "#333"
    # grey "No qualifying change" needs a wider pill
    tw = max(w, min(1.70, 0.20 + 0.078 * len(text)))
    box = FancyBboxPatch((cx - tw / 2, cy - h / 2), tw, h,
                         boxstyle="round,pad=0.012,rounding_size=0.18",
                         fc=fc, ec="none", clip_on=False)
    ax.add_patch(box)
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if style in ("ctx", "indep") else "normal",
            clip_on=False)


def canonical_order(repo: Path, model: str) -> list[str]:
    """Exact Table 4A/4B ordinal list (Position 1..25). Never re-sort by evidence."""
    t4 = pd.read_csv(repo / TABLE4[model])
    # Use file order as written; Position must already be 1..25 sequentially.
    assert list(t4["Position"]) == list(range(1, 26)), (
        f"{model}: Table 4 Position column is not already 1..25 in file order")
    tfs = [strip_shared_mark(x) for x in t4["Candidate regulator"]]
    assert len(tfs) == 25 and len(set(tfs)) == 25, f"{model}: not 25 unique TFs"
    return tfs


def render_model(model: str, by_tf: dict, order: list[str], shared: set, out: Path):
    # Join evidence onto the canonical order WITHOUT re-sorting
    rows = []
    for i, tf in enumerate(order):
        assert tf in by_tf, f"{model}: {tf} missing from audited multiomics table"
        ev = by_tf[tf]
        rows.append({
            "Position": i + 1,
            "TF": tf,
            "tcr": ev["tcr"], "incr": ev["incr"],
            "protein": ev["protein"], "phospho": ev["phospho"], "indep": ev["indep"],
        })
    sub = pd.DataFrame(rows)
    # ASSERT: displayed order exactly matches the canonical ranking table
    assert list(sub["TF"]) == order, f"{model}: displayed order != canonical Table 4 order"

    cols = [
        ("Rank", 0.55),
        ("Candidate regulator", 1.55),
        ("TCR activation\nmRNA", 1.65),
        ("Incremental\nBTLA mRNA", 1.55),
        ("Protein\nabundance", 1.45),
        ("Phosphosite\nabundance", 1.70),
        ("Independent\nassociation", 1.45),
    ]
    gap = 0.10
    x0 = 0.35
    xpos, x = [], x0
    for _, w in cols:
        xpos.append(x + w / 2)
        x += w + gap
    width = x + 0.2

    nrows = 25
    row_h = 0.78
    top = 2.35
    header_h = 1.15
    y_header = top
    y0 = top + header_h
    y_end = y0 + nrows * row_h
    legend_y = y_end + 0.55
    fig_h = legend_y + 3.6

    fig, ax = plt.subplots(figsize=(width * 0.95, fig_h * 0.55))
    ax.set_xlim(0, width)
    ax.set_ylim(fig_h, 0)
    ax.axis("off")

    ax.text(width / 2, 0.55, f"{model} candidate multiomics evidence",
            ha="center", va="center", fontsize=14, fontweight="bold", color=TEXT_DARK)

    ax.add_patch(Rectangle((0.15, y_header), width - 0.3, header_h, fc=HEADER_BG, ec="none"))
    for (lab, _), cx in zip(cols, xpos):
        ax.text(cx, y_header + header_h / 2, lab, ha="center", va="center",
                fontsize=8.0, fontweight="bold", color=TEXT_DARK)
    for yy, lw in ((y_header, 1.35), (y0, 0.9), (y_end, 1.35)):
        ax.plot([0.15, width - 0.15], [yy, yy], color="#222", lw=lw, solid_capstyle="butt")

    n_mrna = n_incr = n_indep = 0
    for i, r in sub.iterrows():
        cy = y0 + i * row_h + row_h / 2
        if i % 2 == 1:
            ax.add_patch(Rectangle((0.15, y0 + i * row_h), width - 0.3, row_h,
                                   fc=ZEBRA, ec="none"))
        ax.text(xpos[0], cy, str(int(r["Position"])), ha="center", va="center",
                fontsize=9, color=TEXT_DARK)
        tf = r["TF"]
        is_shared = tf in shared
        ax.text(xpos[1] - (0.18 if is_shared else 0), cy, tf, ha="center",
                va="center", fontsize=9.2, fontweight="bold", color=TEXT_DARK)
        if is_shared:
            ax.plot(xpos[1] + 0.42, cy, marker="o", markersize=4.5,
                    color="black", linestyle="None")

        for j, key in enumerate(["tcr", "incr", "protein", "phospho", "indep"]):
            lab, style = classify_cell(r[key], key)
            if key == "tcr" and style == "ctx":
                n_mrna += 1
            if key == "incr" and style == "ctx":
                n_incr += 1
            if key == "indep" and style == "indep":
                n_indep += 1
            pill(ax, xpos[j + 2], cy, lab, style,
                 w=0.95 if style == "grey" else 0.72, h=0.50,
                 fs=6.2 if style == "grey" else 6.6)

    # reconciliation assertions
    exp = EXP[model]
    assert n_mrna == exp["mrna"], f"{model}: mRNA {n_mrna} != {exp['mrna']}"
    assert n_incr == exp["incr"], f"{model}: incr {n_incr} != {exp['incr']}"
    assert n_indep == exp["indep"], f"{model}: indep {n_indep} != {exp['indep']}"
    shared_indep = [r["TF"] for _, r in sub.iterrows()
                    if r["TF"] in shared and classify_cell(r["indep"], "indep")[1] == "indep"]
    assert shared_indep == [], f"{model}: shared independent associations {shared_indep}"

    # legend
    ly = legend_y + 0.15
    ax.plot([0.15, width - 0.15], [legend_y - 0.15, legend_y - 0.15], color="#222", lw=1.0)

    def _leg_pill(x, fc, label, tc="white"):
        box = FancyBboxPatch((x, ly), 0.55, 0.38,
                             boxstyle="round,pad=0.01,rounding_size=0.12", fc=fc, ec="none")
        ax.add_patch(box)
        ax.text(x + 0.70, ly + 0.19, label, ha="left", va="center", fontsize=7.3)

    _leg_pill(0.35, BLUE_CTX, "Contextual mRNA (qualifying)")
    _leg_pill(3.55, GREEN_IND, "Independent protein / phosphosite")
    _leg_pill(7.15, GREY_NONE, "No qualifying change", tc="#333")
    ax.text(9.55, ly + 0.19, f"{DASH}  absent / not assessable",
            ha="left", va="center", fontsize=7.3)
    ax.plot(0.55, ly + 0.85, marker="o", markersize=5, color="black")
    ax.text(0.75, ly + 0.85, "shared with the other model's top 25",
            ha="left", va="center", fontsize=7.3)
    ax.text(4.55, ly + 0.85, "Failed QC  =  quality-control failure",
            ha="left", va="center", fontsize=7.3)

    note = (
        "Candidate mRNA is contextual (TCR vs unstimulated control; BTLA+TCR vs TCR) and is shown "
        "in blue. Protein and phosphosite are independent observed molecular layers "
        "(BTLA+TCR vs TCR; source thresholds) and are shown in green. "
        "The combined BTLA+TCR vs unstimulated-control mRNA contrast was not available "
        "from the fitted source models and was not reconstructed. "
        "Direction is not interpreted as supportive or opposing (unsigned rankings). "
        f"Row order matches Table {'4A' if model == 'GREmLN' else '4B'} exactly."
    )
    ax.text(0.25, ly + 1.35, textwrap.fill(note, 118), ha="left", va="top",
            fontsize=6.8, color="#444", linespacing=1.35)
    ax.text(0.25, ly + 2.55,
            f"Figure | {model} top 25 candidate regulator multiomics evidence summary.",
            ha="left", va="top", fontsize=8.5, fontweight="bold", color=TEXT_DARK)

    stem = OUT_STEM[model]
    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(out / f"{stem}{ext}", dpi=400 if ext == ".png" else 300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[multiomics table] {stem}  mRNA={n_mrna} incr={n_incr} indep={n_indep}  OK")


def main() -> int:
    repo = repo_root()
    out = repo / "results/multiomics/audit_v2"
    s7 = pd.read_csv(out / "table_S7_audited_multiomics_43_candidates.csv")
    union = pd.read_csv(repo / "results/publication_data/top25_union_primary.csv")

    by_tf = {}
    for _, r in s7.iterrows():
        by_tf[r["Candidate regulator"]] = {
            "tcr": r["TCR activation mRNA"],
            "incr": r["Incremental BTLA mRNA"],
            "protein": r["Protein abundance"],
            "phospho": r["Phosphosite abundance"],
            "indep": r["Independent molecular association"],
        }
    shared = set(union.loc[union["in_gremln_top25"] & union["in_genie3_top25"], "TF"])

    for model in ("GREmLN", "GENIE3"):
        order = canonical_order(repo, model)
        render_model(model, by_tf, order, shared, out)

    # cross-model: zero shared independent associations
    shared_with_indep = [
        tf for tf in shared
        if classify_cell(by_tf[tf]["indep"], "indep")[1] == "indep"
    ]
    assert shared_with_indep == [], f"shared independent associations: {shared_with_indep}"
    print("[assert] shared independent associations = 0  OK")

    cap = (
        "Figure | Model-specific audited multiomics evidence summaries.\n\n"
        "Matched Paperclip-style tables for the GREmLN and GENIE3 seed-excluded "
        "top-25 candidate regulators. Row order matches Table 4A / 4B exactly. "
        "Blue pills mark qualifying contextual candidate mRNA; green pills mark "
        "qualifying independent protein or phosphosite associations. Grey pills "
        "read 'No qualifying change'; an em dash marks absence from the assay "
        "universe. Black dots mark candidates shared between the two models' "
        "top-25 lists. Direction is not interpreted as supportive or opposing "
        "because both rankings are unsigned.\n"
    )
    (out / "fig_model_multiomics_evidence_caption.txt").write_text(cap)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
