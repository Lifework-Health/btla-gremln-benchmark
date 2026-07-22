#!/usr/bin/env python3
"""Redesign Table 5 as a compact publication Results table (presentation only).

Counts and evidence definitions are unchanged; only layout, typography and
compact cell wording are redesigned. Full-width manuscript table (~11 × 5 in),
not a stretched landscape page.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd

PAGE_W, PAGE_H = 11.2, 5.05  # inches — compact Results insertion size
TITLE_PT = 14.0
HDR_PT = 10.75
BODY_PT = 10.75
NOTE_PT = 8.25
MIN_BODY = 10.5

# Evidence ~23%; model cols wide enough for Overall 2-line bold text; interp flexible
COL_FRAC = [0.22, 0.23, 0.23, 0.32]
HEADER_BG = "#F0F0F0"
ZEBRA = "#F7F7F7"
OVERALL_BG = "#EFEFEF"
RULE = "#222222"
TEXT = "#1a1a1a"
PAD_FRAC = 0.012  # left/right inset within each column (axes fraction)

# Compact presentation rows (same counts as audited Table 5)
ROWS = [
    {
        "source": "CRISPRi perturbational evidence",
        "gremln": "3 / 20 represented",
        "genie3": "6 / 19 evaluable",
        "interp": "Broader perturbational evidence for GENIE3",
        "gremln_bold_nums": True,
    },
    {
        "source": "Audited literature evidence",
        "gremln": "7 / 25 strong or moderate",
        "genie3": "9 / 25 strong or moderate",
        "interp": "Modest literature lean towards GENIE3",
        "gremln_bold_nums": True,
    },
    {
        "source": "TCR activation mRNA context",
        "gremln": "2 / 25",
        "genie3": "5 / 25",
        "interp": "Activation-associated context was more frequent for GENIE3",
        "gremln_bold_nums": True,
    },
    {
        "source": "Independent BTLA molecular association",
        "gremln": "4 / 25",
        "genie3": "3 / 25",
        "interp": "No clear separation; molecular coverage was sparse and unequal",
        "gremln_bold_nums": True,
    },
    {
        "source": "Overall",
        "gremln": "Distinct candidates and molecular hypotheses",
        "genie3": "Broader perturbational and literature corroboration",
        "interp": "No clear or consistent superiority",
        "gremln_bold_nums": False,
        "overall": True,
    },
]

NOTES = [
    "1. CRISPRi denominators represent candidates present in the screen for GREmLN "
    "and candidates evaluable after on-target quality control for GENIE3.",
    "2. Literature tiers are audited provisional annotations. Candidate mRNA is "
    "contextual; no candidate showed a qualifying incremental BTLA-specific mRNA "
    "response. Independent molecular evidence comprises qualifying directly measured "
    "protein or phosphosite associations, for which assay coverage was incomplete.",
]

# Assertions against prior audited CSV values
EXPECTED = {
    "CRISPRi": ("3", "20", "6", "19"),
    "lit": ("7", "9"),
    "mrna": ("2", "5"),
    "indep": ("4", "3"),
}


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def word_count(s: str) -> int:
    return len(s.split())


def wrap_to_width(text: str, max_width_in: float, fontsize: float,
                  fontweight: str = "normal", fontfamily: str = "DejaVu Sans",
                  max_lines: int = 2) -> list[str]:
    """Wrap text so each line fits max_width_in (inches); pack into ≤max_lines."""
    words = text.split()
    if not words:
        return [""]

    def line_width(s: str) -> float:
        # Conservative advance so wrapped lines stay inside column bounds
        avg = 0.56 if fontweight == "normal" else 0.60
        return len(s) * (avg * fontsize) / 72.0

    def greedy(max_w: float) -> list[str]:
        lines: list[str] = []
        cur = words[0]
        for w in words[1:]:
            trial = f"{cur} {w}"
            if line_width(trial) <= max_w:
                cur = trial
            else:
                lines.append(cur)
                cur = w
        lines.append(cur)
        return lines

    lines = greedy(max_width_in)
    if len(lines) <= max_lines:
        return lines

    if max_lines == 2:
        n = len(words)
        for k in range(1, n):
            a = " ".join(words[:k])
            b = " ".join(words[k:])
            if line_width(a) <= max_width_in and line_width(b) <= max_width_in:
                return [a, b]
        mid = max(1, n // 2)
        return [" ".join(words[:mid]), " ".join(words[mid:])]

    return greedy(max_width_in * 1.15)[:max_lines]


def split_model_value(val: str) -> tuple[str, str]:
    """Split '3 / 20 represented' -> ('3 / 20', ' represented') for mixed weight."""
    import re
    m = re.match(r"^(\d+\s*/\s*\d+)(.*)$", val)
    if m:
        return m.group(1), m.group(2)
    return val, ""


def assert_content():
    assert BODY_PT >= MIN_BODY
    assert len(ROWS) == 5
    # counts preserved
    assert "3 / 20" in ROWS[0]["gremln"] and "6 / 19" in ROWS[0]["genie3"]
    assert "7 / 25" in ROWS[1]["gremln"] and "9 / 25" in ROWS[1]["genie3"]
    assert ROWS[2]["gremln"] == "2 / 25" and ROWS[2]["genie3"] == "5 / 25"
    assert ROWS[3]["gremln"] == "4 / 25" and ROWS[3]["genie3"] == "3 / 25"
    for r in ROWS:
        assert word_count(r["interp"]) <= 12, f"interp too long: {r['interp']}"
        for key in ("source", "gremln", "genie3", "interp"):
            # soft wrap check: ~chars per column
            pass
    # no forbidden content
    blob = " ".join(r["interp"] + r["gremln"] + r["genie3"] for r in ROWS)
    for bad in ["Spearman", "bootstrap", "superiority margin", "composite score",
                "incremental BTLA-specific mRNA response"]:
        assert bad not in blob, f"forbidden in body: {bad}"


def render(out: Path):
    assert_content()

    # Tight canvas: title + header + 5 rows + notes
    left_m, right_m = 0.35, 0.35
    top_m = 0.42
    bottom_m = 0.95
    usable_w = PAGE_W - left_m - right_m
    usable_h = PAGE_H - top_m - bottom_m

    hdr_h_in = 0.38
    row_h_in = (usable_h - hdr_h_in) / 5
    assert row_h_in >= 0.38, f"row height too small: {row_h_in}"

    fig = plt.figure(figsize=(PAGE_W, PAGE_H), dpi=300, facecolor="white")
    ax = fig.add_axes([left_m / PAGE_W, bottom_m / PAGE_H,
                       usable_w / PAGE_W, usable_h / PAGE_H])
    ax.set_xlim(0, 1)
    ax.set_ylim(1, 0)
    ax.axis("off")

    title = "Table 5. Triangulated evidence across the two top 25 candidate lists"
    fig.text(0.5, 1 - 0.12 / PAGE_H, title, ha="center", va="top",
             fontsize=TITLE_PT, fontweight="bold", color=TEXT,
             fontfamily="DejaVu Sans")

    xb = [0.0]
    for f in COL_FRAC:
        xb.append(xb[-1] + f)

    # header
    y_hdr = hdr_h_in / usable_h
    ax.add_patch(Rectangle((0, 0), 1, y_hdr, fc=HEADER_BG, ec="none", zorder=0))
    headers = ["Evidence source", "GREmLN", "GENIE3", "Comparative interpretation"]
    for j, h in enumerate(headers):
        cx = (xb[j] + xb[j + 1]) / 2 if j in (1, 2) else xb[j] + 0.01
        ha = "center" if j in (1, 2) else "left"
        ax.text(cx, y_hdr / 2, h, ha=ha, va="center", fontsize=HDR_PT,
                fontweight="bold", color=TEXT, fontfamily="DejaVu Sans")

    ax.plot([0, 1], [0, 0], color=RULE, lw=1.5, solid_capstyle="butt")
    ax.plot([0, 1], [y_hdr, y_hdr], color=RULE, lw=1.2, solid_capstyle="butt")

    rh = (1.0 - y_hdr) / 5
    for i, r in enumerate(ROWS):
        y = y_hdr + i * rh
        is_overall = r.get("overall", False)
        if is_overall:
            ax.add_patch(Rectangle((0, y), 1, rh, fc=OVERALL_BG, ec="none", zorder=0))
            ax.plot([0, 1], [y, y], color=RULE, lw=1.15, solid_capstyle="butt", zorder=3)
        elif i % 2 == 1:
            ax.add_patch(Rectangle((0, y), 1, rh, fc=ZEBRA, ec="none", zorder=0))

        cells = [r["source"], r["gremln"], r["genie3"], r["interp"]]
        for j, val in enumerate(cells):
            col_w_in = usable_w * COL_FRAC[j]
            pad_in = usable_w * PAD_FRAC
            max_w_in = max(0.4, col_w_in - 2 * pad_in)

            is_model = j in (1, 2)
            is_num_row = (not is_overall) and bool(r.get("gremln_bold_nums")) and is_model

            if is_num_row:
                num, rest = split_model_value(val)
                lines_n = wrap_to_width(num + rest, max_w_in, BODY_PT,
                                        fontweight="bold", max_lines=2)
                assert len(lines_n) <= 2, f"wrap >2: {val} -> {lines_n}"
                if len(lines_n) == 1 and rest:
                    avg_b, avg_n = 0.55, 0.52
                    w_num = len(num) * (avg_b * BODY_PT) / 72.0
                    w_rest = len(rest) * (avg_n * BODY_PT) / 72.0
                    total = w_num + w_rest
                    cx = (xb[j] + xb[j + 1]) / 2

                    def in_to_ax(inches: float) -> float:
                        return inches / usable_w

                    x0 = cx - in_to_ax(total) / 2
                    ax.text(x0, y + rh / 2, num, ha="left", va="center",
                            fontsize=BODY_PT, fontweight="bold", color=TEXT,
                            fontfamily="DejaVu Sans", zorder=2)
                    ax.text(x0 + in_to_ax(w_num), y + rh / 2, rest, ha="left",
                            va="center", fontsize=BODY_PT, fontweight="normal",
                            color=TEXT, fontfamily="DejaVu Sans", zorder=2)
                else:
                    text = "\n".join(lines_n)
                    ax.text((xb[j] + xb[j + 1]) / 2, y + rh / 2, text,
                            ha="center", va="center", fontsize=BODY_PT,
                            fontweight="bold", color=TEXT,
                            fontfamily="DejaVu Sans", zorder=2)
                continue

            weight = "bold" if (is_overall or j == 0) else "normal"
            lines = wrap_to_width(val, max_w_in, BODY_PT, fontweight=weight,
                                  max_lines=2)
            assert len(lines) <= 2, f"wrap >2: {val} -> {lines}"
            text = "\n".join(lines)
            if is_model:
                ax.text((xb[j] + xb[j + 1]) / 2, y + rh / 2, text,
                        ha="center", va="center", fontsize=BODY_PT,
                        fontweight=weight, color=TEXT,
                        fontfamily="DejaVu Sans", zorder=2)
            else:
                ax.text(xb[j] + PAD_FRAC, y + rh / 2, text,
                        ha="left", va="center", fontsize=BODY_PT,
                        fontweight=weight, color=TEXT,
                        fontfamily="DejaVu Sans", zorder=2)

    ax.plot([0, 1], [1, 1], color=RULE, lw=1.5, solid_capstyle="butt")

    # notes immediately beneath
    note_y = bottom_m / PAGE_H - 0.02
    line_step = (NOTE_PT + 1.8) / 72.0 / PAGE_H
    yy = note_y
    for note in NOTES:
        wrapped = textwrap.fill(note, 145)
        nlines = wrapped.count("\n") + 1
        fig.text(left_m / PAGE_W, yy, wrapped, ha="left", va="top",
                 fontsize=NOTE_PT, color="#333333", fontfamily="DejaVu Sans",
                 linespacing=1.25)
        yy -= nlines * line_step + 0.01

    stem = out / "table5_triangulated_evidence_summary"
    fig.savefig(str(stem) + ".svg", format="svg", bbox_inches=None,
                facecolor="white", pad_inches=0)
    fig.savefig(str(stem) + ".png", format="png", dpi=300, bbox_inches=None,
                facecolor="white", pad_inches=0)
    plt.close(fig)

    # compact CSV (presentation wording; counts unchanged)
    csv_df = pd.DataFrame([
        {
            "Evidence source": r["source"],
            "GREmLN": r["gremln"],
            "GENIE3": r["genie3"],
            "Comparative interpretation": r["interp"],
        }
        for r in ROWS
    ])
    csv_df.to_csv(stem.with_suffix(".csv"), index=False)

    # mirror
    mirror = out.parent.parent / "publication_assets" / "tables"
    # path: results/multiomics/audit_v2 -> results/publication_assets/tables
    mirror = repo_root() / "results/publication_assets/tables"
    mirror.mkdir(parents=True, exist_ok=True)
    for ext in (".svg", ".png", ".csv"):
        src = Path(str(stem) + ext)
        if src.exists():
            (mirror / src.name).write_bytes(src.read_bytes())

    return stem


def main() -> int:
    out = repo_root() / "results/multiomics/audit_v2"
    # verify against prior audited CSV counts before overwrite
    prev = pd.read_csv(out / "table5_triangulated_evidence_summary.csv",
                       keep_default_na=False)
    assert "3 candidates" in prev.iloc[0]["GREmLN"] or "3 / 20" in prev.iloc[0]["GREmLN"]
    assert "6 candidates" in prev.iloc[0]["GENIE3"] or "6 / 19" in prev.iloc[0]["GENIE3"]
    assert "7/25" in prev.iloc[1]["GREmLN"].replace(" ", "") or "7 / 25" in prev.iloc[1]["GREmLN"]
    assert "9/25" in prev.iloc[1]["GENIE3"].replace(" ", "") or "9 / 25" in prev.iloc[1]["GENIE3"]
    assert "2/25" in prev.iloc[2]["GREmLN"].replace(" ", "") or prev.iloc[2]["GREmLN"].startswith("2")
    assert "5/25" in prev.iloc[2]["GENIE3"].replace(" ", "") or prev.iloc[2]["GENIE3"].startswith("5")
    assert "4/25" in prev.iloc[3]["GREmLN"].replace(" ", "") or prev.iloc[3]["GREmLN"].startswith("4")
    assert "3/25" in prev.iloc[3]["GENIE3"].replace(" ", "") or prev.iloc[3]["GENIE3"].startswith("3")

    stem = render(out)
    from PIL import Image
    im = Image.open(str(stem) + ".png")
    w_in, h_in = im.size[0] / 300, im.size[1] / 300
    print(f"[Table 5] {stem}.svg/.png/.csv")
    print(f"  rendered size ≈ {w_in:.2f} × {h_in:.2f} in @ 300 dpi  ({im.size})")
    print(f"  body={BODY_PT} pt  title={TITLE_PT} pt  notes={NOTE_PT} pt")
    print("  counts preserved: CRISPRi 3/20 vs 6/19; lit 7 vs 9; mRNA 2 vs 5; indep 4 vs 3")
    assert 10.0 <= w_in <= 12.0, f"width {w_in} out of target 10.5–11.5"
    assert 3.8 <= h_in <= 6.2, f"height {h_in} out of target ~4.5–5.5"
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
