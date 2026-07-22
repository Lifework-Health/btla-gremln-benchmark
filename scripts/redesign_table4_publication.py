#!/usr/bin/env python3
"""Redesign Tables 4A/4B for publication readability (presentation only).

Reads the audited source CSVs produced by build_tables_4_and_5_evidence.py and
renders landscape journal-style SVG/PNG. Evidence calls are not recomputed or
altered — only display labels, typography and layout change.

Primary vector: SVG. Raster: PNG @ 300 dpi. Page: 14 × 8.5 in landscape.
If 25 rows cannot fit at >=9.5 pt with <=2-line wraps, split into pages
1–13 and 14–25 rather than shrinking the font.
"""
from __future__ import annotations

import re
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

SUP_S = "\u02e2"
PAGE_W, PAGE_H = 14.0, 8.5  # inches

# Colours (restrained journal highlighting)
PG = "#E3F1E1"          # pale green — detected CRISPRi / qualifying protein-phospho
PG_STRONG = "#EDF7ED"   # very light green — Strong literature only
PYEL = "#FFF8E7"        # pale yellow — qualifying contextual mRNA
GREY = "#EEEEEE"        # not screened / failed QC / not measured
ZEBRA = "#F7F7F7"
HEADER_BG = "#F2F2F2"
RULE = "#222222"
TEXT = "#1a1a1a"

# Typography (pt). Body must never drop below 9.5.
TITLE_PT = 15.5
HDR_PT = 10.75
BODY_PT = 10.25
NOTE_PT = 8.5
MIN_BODY_PT = 9.5

# Column fractions
COL_FRAC = [0.05, 0.12, 0.24, 0.11, 0.20, 0.28]
COL_KEYS = ["rank", "cand", "crispri", "lit", "mrna", "indep"]

NOTES = [
    "1. S, shared between both top-25 lists.",
    "2. CRISPRi values show responsive predicted targets / tested predicted targets; "
    "absent and failed-QC perturbations were not converted to zero.",
    "3. Literature tiers are audited provisional Paperclip v2 annotations; strong and "
    "moderate tiers remain pending documented human verification.",
    "4. Candidate mRNA is contextual and shows qualifying TCR-versus-unstimulated-control "
    "responses. No candidate had a qualifying incremental BTLA+TCR-versus-TCR mRNA response.",
    "5. Protein and phosphosite entries are qualifying directly measured molecular "
    "associations in BTLA+TCR versus TCR. Candidates not measured were not treated as negative evidence.",
    "6. Molecular direction was not interpreted as supportive or opposing because both "
    "model rankings were unsigned.",
]

EXP = {
    "GREmLN": {"crispri": 3, "lit_sm": 7, "mrna": 2, "indep": 4},
    "GENIE3": {"crispri": 6, "lit_sm": 9, "mrna": 5, "indep": 3},
}


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def strip_s(name: str) -> str:
    return re.sub(r"\u02e2$", "", str(name)).strip()


# ---------- compact display mapping (presentation only) ----------

def compact_crispri(v: str) -> str:
    v = str(v)
    m = re.search(r"\((\d+)/(\d+)\)", v)
    frac = f" ({m.group(1)}/{m.group(2)})" if m else ""
    if "anti-concordant" in v:
        return f"Anti-concordant{frac}"
    if "concordant" in v and "anti" not in v:
        return f"Concordant{frac}"
    if v.startswith("Detected") and "mixed" in v.lower():
        return f"Mixed{frac}"
    if v.startswith("No detected") or v.startswith("No response"):
        return f"No response{frac}" if frac else "No response"
    if "Not represented" in v or v == "Not screened":
        return "Not screened"
    if "Failed" in v:
        return "Failed QC"
    return v


def compact_mrna(v: str) -> str:
    v = str(v)
    if v.startswith("TCR activation"):
        # "TCR activation ↑: 1 h, 24 h" -> "TCR ↑ 1 h, 24 h"
        v = v.replace("TCR activation ", "TCR ").replace(": ", " ")
        return v
    if "Absent" in v or v == "Not measured":
        return "Not measured"
    if "Failed" in v:
        return "Failed QC"
    return "No qualifying change"


def compact_indep(v: str) -> str:
    v = str(v)
    if v.startswith("Protein association:"):
        return "Protein: " + v.split(":", 1)[1].strip()
    if v.startswith("Phosphosite association:"):
        return "Phospho: " + v.split(":", 1)[1].strip()
    if "Protein association:" in v and "Phosphosite" in v:
        # combined — keep both short
        parts = []
        for bit in v.split(";"):
            bit = bit.strip()
            if bit.startswith("Protein"):
                parts.append("Protein: " + bit.split(":", 1)[1].strip())
            elif bit.startswith("Phosphosite"):
                parts.append("Phospho: " + bit.split(":", 1)[1].strip())
        return "; ".join(parts)
    if "Not measured" in v:
        return "Not measured"
    if "Failed" in v:
        return "Failed QC"
    return "No qualifying signal"


def to_display(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        raw = str(r["Candidate regulator"])
        shared = SUP_S in raw
        tf = strip_s(raw)
        rows.append({
            "rank": int(r["Position"]),
            "cand": tf,
            "shared": shared,
            "crispri": compact_crispri(r["CRISPRi evidence"]),
            "lit": str(r["Audited literature evidence"]),
            "mrna": compact_mrna(r["Candidate mRNA context"]),
            "indep": compact_indep(r["Independent BTLA molecular evidence"]),
            # keep source for assertions
            "_src_crispri": str(r["CRISPRi evidence"]),
            "_src_lit": str(r["Audited literature evidence"]),
            "_src_mrna": str(r["Candidate mRNA context"]),
            "_src_indep": str(r["Independent BTLA molecular evidence"]),
        })
    return pd.DataFrame(rows)


def cell_bg(key: str, val: str) -> str | None:
    if key == "crispri":
        if val.startswith("Concordant") or val.startswith("Anti-concordant") or val.startswith("Mixed"):
            return PG
        if val in ("Not screened", "Failed QC"):
            return GREY
    if key == "lit":
        if val == "Strong":
            return PG_STRONG
        # Moderate: no fill; None/Weak: no fill
    if key == "mrna":
        if val.startswith("TCR"):
            return PYEL
        if val in ("Not measured", "Failed QC"):
            return GREY
    if key == "indep":
        if val.startswith("Protein:") or val.startswith("Phospho:"):
            return PG
        if val in ("Not measured", "Failed QC"):
            return GREY
    return None


def wrap_lines(text: str, width_chars: int) -> list[str]:
    if not text:
        return [""]
    return textwrap.wrap(text, width=max(4, width_chars), break_long_words=False,
                         break_on_hyphens=False) or [text]


def estimate_wrap(text: str, col_frac: float, page_w_in: float, fontsize: float,
                  usable_width_frac: float = 0.92) -> int:
    """Estimate number of wrapped lines for a cell given font size and column fraction."""
    # approximate char width ≈ 0.5 * fontsize points; 1 inch = 72 pt
    usable_in = page_w_in * usable_width_frac * col_frac
    char_w_in = (0.50 * fontsize) / 72.0
    nchars = max(4, int(usable_in / char_w_in) - 1)
    return len(wrap_lines(text, nchars)), nchars


def assert_wrap_limits(disp: pd.DataFrame, body_pt: float, page_w: float):
    """Fail if any body/header cell would wrap to >2 lines at the chosen font."""
    headers_l1 = ["Rank", "Candidate", "CRISPRi", "Literature", "mRNA context", "Protein / phosphosite"]
    # grouped header "BTLA experimental evidence" spans cols 4+5
    issues = []
    for j, h in enumerate(headers_l1):
        n, _ = estimate_wrap(h, COL_FRAC[j], page_w, HDR_PT)
        if n > 2:
            issues.append(f"header '{h}' wraps to {n} lines")
    n, _ = estimate_wrap("BTLA experimental evidence", COL_FRAC[4] + COL_FRAC[5], page_w, HDR_PT)
    if n > 2:
        issues.append(f"grouped header wraps to {n} lines")

    for _, r in disp.iterrows():
        for j, key in enumerate(COL_KEYS):
            text = str(r[key])
            if key == "cand" and r["shared"]:
                text = text + "S"
            n, _ = estimate_wrap(text, COL_FRAC[j], page_w, body_pt)
            if n > 2:
                issues.append(f"row {r['rank']} {key}='{text}' wraps to {n} lines")
    if body_pt < MIN_BODY_PT:
        issues.append(f"body font {body_pt} pt < {MIN_BODY_PT} pt minimum")
    if issues:
        raise AssertionError("wrap/typography assertions failed:\n  - " + "\n  - ".join(issues))


def render_page(disp: pd.DataFrame, title: str, caption: str, out_stem: Path,
                page_tag: str = ""):
    """Render one landscape page of the table."""
    assert BODY_PT >= MIN_BODY_PT

    n = len(disp)
    # layout geometry in figure inches → normalised later via axes inches
    left_m, right_m = 0.40, 0.40
    top_m = 0.48
    bottom_notes = 1.35
    usable_w = PAGE_W - left_m - right_m
    usable_h = PAGE_H - top_m - bottom_notes

    # header block height (two rows)
    hdr_h = 0.48
    body_h = usable_h - hdr_h
    row_h = body_h / max(n, 1)
    if row_h < 0.22:
        raise AssertionError(
            f"row height {row_h:.3f} in too small for {n} rows; split the table instead of shrinking font")

    fig = plt.figure(figsize=(PAGE_W, PAGE_H), dpi=300, facecolor="white")
    ax = fig.add_axes([left_m / PAGE_W, bottom_notes / PAGE_H,
                       usable_w / PAGE_W, usable_h / PAGE_H])
    ax.set_xlim(0, 1)
    ax.set_ylim(1, 0)
    ax.axis("off")

    # title (figure coordinates)
    fig.text(0.5, 1 - 0.28 / PAGE_H, title, ha="center", va="top",
             fontsize=TITLE_PT, fontweight="bold", color=TEXT,
             fontfamily="DejaVu Sans")

    # column x edges
    xb = np.concatenate([[0.0], np.cumsum(COL_FRAC)])
    assert abs(xb[-1] - 1.0) < 1e-9

    # ---- header ----
    # top header band
    y0 = 0.0
    y_mid = hdr_h / usable_h * 0.48   # first header row depth (as fraction of axes)
    y_hdr_end = hdr_h / usable_h

    ax.add_patch(Rectangle((0, y0), 1, y_hdr_end, fc=HEADER_BG, ec="none", zorder=0))

    # Rank, Candidate, CRISPRi, Literature — span both header rows
    top_labels = [
        (0, 1, "Rank"),
        (1, 2, "Candidate"),
        (2, 3, "CRISPRi"),
        (3, 4, "Literature"),
    ]
    for i0, i1, lab in top_labels:
        cx = (xb[i0] + xb[i1]) / 2
        ax.text(cx, y_hdr_end / 2, lab, ha="center", va="center",
                fontsize=HDR_PT, fontweight="bold", color=TEXT, fontfamily="DejaVu Sans")

    # Grouped header: BTLA experimental evidence spanning cols 4–5
    ax.text((xb[4] + xb[6]) / 2, y_mid / 2, "BTLA experimental evidence",
            ha="center", va="center", fontsize=HDR_PT, fontweight="bold",
            color=TEXT, fontfamily="DejaVu Sans")
    # light rule under grouped title
    ax.plot([xb[4], xb[6]], [y_mid, y_mid], color="#888", lw=0.6, solid_capstyle="butt")
    # second-level headers
    for i0, i1, lab in [(4, 5, "mRNA context"), (5, 6, "Protein / phosphosite")]:
        cx = (xb[i0] + xb[i1]) / 2
        ax.text(cx, (y_mid + y_hdr_end) / 2, lab, ha="center", va="center",
                fontsize=HDR_PT - 0.3, fontweight="bold", color=TEXT,
                fontfamily="DejaVu Sans")

    # rules
    ax.plot([0, 1], [0, 0], color=RULE, lw=1.5, solid_capstyle="butt")
    ax.plot([0, 1], [y_hdr_end, y_hdr_end], color=RULE, lw=1.2, solid_capstyle="butt")

    # ---- body ----
    body_top = y_hdr_end
    body_bot = 1.0
    rh = (body_bot - body_top) / n

    # precompute wrap widths in chars
    wrap_w = []
    for j, fr in enumerate(COL_FRAC):
        _, nc = estimate_wrap("x", fr, usable_w, BODY_PT)
        wrap_w.append(nc)

    max_lines_seen = 1
    for i in range(n):
        r = disp.iloc[i]
        y = body_top + i * rh
        if i % 2 == 1:
            ax.add_patch(Rectangle((0, y), 1, rh, fc=ZEBRA, ec="none", zorder=0))

        values = [
            str(r["rank"]),
            r["cand"],
            r["crispri"],
            r["lit"],
            r["mrna"],
            r["indep"],
        ]
        for j, (key, val) in enumerate(zip(COL_KEYS, values)):
            bg = cell_bg(key, val if key != "cand" else "")
            if key != "cand" and bg:
                ax.add_patch(Rectangle((xb[j], y), COL_FRAC[j], rh,
                                       fc=bg, ec="none", zorder=1))

            lines = wrap_lines(val, wrap_w[j])
            max_lines_seen = max(max_lines_seen, len(lines))
            if len(lines) > 2:
                raise AssertionError(
                    f"cell wrap >2 lines at rank {r['rank']} col {key}: {lines}")

            if key == "cand" and r["shared"]:
                text = lines[0] + SUP_S
                fw = "bold"
            else:
                text = "\n".join(lines)
                fw = "normal"

            if key == "rank":
                ax.text((xb[j] + xb[j + 1]) / 2, y + rh / 2, text,
                        ha="center", va="center", fontsize=BODY_PT,
                        color=TEXT, fontfamily="DejaVu Sans", fontweight=fw, zorder=2)
            else:
                ax.text(xb[j] + 0.008, y + rh / 2, text,
                        ha="left", va="center", fontsize=BODY_PT,
                        color=TEXT, fontfamily="DejaVu Sans", fontweight=fw, zorder=2)

    ax.plot([0, 1], [1, 1], color=RULE, lw=1.5, solid_capstyle="butt")

    # ---- notes (two columns) beneath the axes ----
    note_top = bottom_notes / PAGE_H - 0.015
    col_x = [left_m / PAGE_W, 0.515]
    line_step = (NOTE_PT + 2.0) / 72.0 / PAGE_H
    for xi, group in zip(col_x, (NOTES[:3], NOTES[3:])):
        yy = note_top
        for note in group:
            wrapped = textwrap.fill(note, 85)
            nlines = wrapped.count("\n") + 1
            fig.text(xi, yy, wrapped, ha="left", va="top", fontsize=NOTE_PT,
                     color="#333333", fontfamily="DejaVu Sans", linespacing=1.3)
            yy -= nlines * line_step + 0.008

    # caption (short) at very bottom
    fig.text(left_m / PAGE_W, 0.018, caption, ha="left", va="bottom",
             fontsize=NOTE_PT + 0.5, fontweight="bold", color=TEXT,
             fontfamily="DejaVu Sans")

    stem = str(out_stem) + (f"_{page_tag}" if page_tag else "")
    fig.savefig(stem + ".svg", format="svg", bbox_inches="tight",
                facecolor="white", pad_inches=0.15)
    fig.savefig(stem + ".png", format="png", dpi=300, bbox_inches="tight",
                facecolor="white", pad_inches=0.15)
    plt.close(fig)
    return max_lines_seen, BODY_PT


def scientific_assertions(disp: pd.DataFrame, model: str, order: list[str]):
    assert len(disp) == 25, f"{model}: expected 25 rows"
    assert list(disp["rank"]) == list(range(1, 26))
    assert list(disp["cand"]) == order, f"{model}: order changed after join"

    n_cr = int(disp["crispri"].str.match(r"^(Concordant|Anti-concordant|Mixed)").sum())
    n_sm = int(disp["lit"].isin(["Strong", "Moderate"]).sum())
    n_mrna = int(disp["mrna"].str.startswith("TCR").sum())
    n_indep = int(disp["indep"].str.match(r"^(Protein:|Phospho:)").sum())
    exp = EXP[model]
    assert n_cr == exp["crispri"], f"{model} crispri {n_cr}!={exp['crispri']}"
    assert n_sm == exp["lit_sm"], f"{model} lit {n_sm}!={exp['lit_sm']}"
    assert n_mrna == exp["mrna"], f"{model} mrna {n_mrna}!={exp['mrna']}"
    assert n_indep == exp["indep"], f"{model} indep {n_indep}!={exp['indep']}"

    blob = " ".join(disp["indep"].astype(str) + disp["mrna"].astype(str) + disp["crispri"].astype(str))
    for bad in ["Inferred regulon", "BIONIC", "coIP", "Co-IP", "kinase", "synapse",
                "Detected —", "Phosphosite association", "Absent from"]:
        assert bad not in blob, f"legacy/verbose '{bad}' remains in display"
    assert not re.search(r"[A-Za-z0-9]+_[STY]\d+", blob)
    assert "REDACTED" not in blob
    assert BODY_PT >= MIN_BODY_PT


def render_model(model: str, csv_path: Path, out_dir: Path, canonical_order: list[str]):
    src = pd.read_csv(csv_path, keep_default_na=False)
    disp = to_display(src)
    # order assertion against canonical
    scientific_assertions(disp, model, canonical_order)
    assert_wrap_limits(disp, BODY_PT, PAGE_W - 0.9)

    title = (f"Table 4A. Triangulated evidence across GREmLN’s top 25 candidate regulators"
             if model == "GREmLN" else
             f"Table 4B. Triangulated evidence across GENIE3’s top 25 candidate regulators")
    caption = title if title.endswith(".") else title + "."
    # Fix caption: title already has period intent
    caption = title.rstrip(".") + "."

    stem = out_dir / ("table4a_gremln_evidence" if model == "GREmLN"
                      else "table4b_genie3_evidence")

    # Try single page; if row height too small, split.
    usable_h = PAGE_H - 0.48 - 1.35
    hdr_h = 0.48
    rh = (usable_h - hdr_h) / 25
    if rh < 0.22:
        # split 1–13 / 14–25
        pages = 2
        page_info = []
        for tag, sl in [("p1", disp.iloc[:13]), ("p2", disp.iloc[13:])]:
            t = title + (f" (ranks {sl.iloc[0]['rank']}–{sl.iloc[-1]['rank']})")
            ml, pt = render_page(sl.reset_index(drop=True), t, caption, stem, page_tag=tag)
            page_info.append((Path(str(stem) + f"_{tag}"), ml))
        import shutil
        for ext in (".svg", ".png"):
            shutil.copy(str(stem) + "_p1" + ext, str(stem) + ext)
    else:
        max_lines, pt = render_page(disp, title, caption, stem)
        pages = 1
        page_info = [(stem, max_lines)]

    return {
        "model": model,
        "pages": pages,
        "body_pt": BODY_PT,
        "title_pt": TITLE_PT,
        "header_pt": HDR_PT,
        "note_pt": NOTE_PT,
        "max_body_lines": max(m for _, m in page_info),
        "page_stems": [str(p) for p, _ in page_info],
    }


def load_order(repo: Path, model: str) -> list[str]:
    stem = ("table4a_gremln_top25_evidence" if model == "GREmLN"
            else "table4b_genie3_top25_evidence")
    for p in [
        repo / "results/multiomics/audit_v2" / f"{stem}.previous_publication.csv",
        repo / "results/publication_assets/tables" / f"{stem}.previous_publication.csv",
        repo / "results/multiomics/audit_v2" / ("table4a_gremln_evidence.csv" if model == "GREmLN"
                                                 else "table4b_genie3_evidence.csv"),
    ]:
        if p.exists():
            d = pd.read_csv(p, keep_default_na=False)
            col = "Candidate regulator" if "Candidate regulator" in d.columns else None
            if col:
                return [strip_s(x) for x in d[col]]
            if "cand" in d.columns:
                return list(d["cand"])
    raise FileNotFoundError(f"canonical order not found for {model}")


def main() -> int:
    repo = repo_root()
    out = repo / "results/multiomics/audit_v2"
    reports = []

    # shared independent = 0
    shared_indep = 0
    for model, csv_name in [
        ("GREmLN", "table4a_gremln_evidence.csv"),
        ("GENIE3", "table4b_genie3_evidence.csv"),
    ]:
        order = load_order(repo, model)
        # Prefer audited CSV order (must match canonical)
        src = pd.read_csv(out / csv_name, keep_default_na=False)
        src_order = [strip_s(x) for x in src["Candidate regulator"]]
        assert src_order == order, f"{model}: CSV order != canonical publication ranking"
        info = render_model(model, out / csv_name, out, order)
        reports.append(info)
        disp = to_display(src)
        # shared indep check across both — accumulate GENIE3 after
        if model == "GENIE3":
            pass

    # cross-model shared independent
    a = to_display(pd.read_csv(out / "table4a_gremln_evidence.csv", keep_default_na=False))
    b = to_display(pd.read_csv(out / "table4b_genie3_evidence.csv", keep_default_na=False))
    shared = set(a.loc[a["shared"], "cand"]) & set(b.loc[b["shared"], "cand"])
    for tf in shared:
        ia = a.loc[a["cand"] == tf, "indep"].iloc[0]
        ib = b.loc[b["cand"] == tf, "indep"].iloc[0]
        assert ia == ib
        assert not str(ia).startswith(("Protein:", "Phospho:")), f"shared indep: {tf}"
        shared_indep += int(str(ia).startswith(("Protein:", "Phospho:")))
    assert shared_indep == 0

    # mirror to publication_assets
    mirror = repo / "results/publication_assets/tables"
    mirror.mkdir(parents=True, exist_ok=True)
    for info in reports:
        for stem in info["page_stems"]:
            for ext in (".svg", ".png"):
                src = Path(stem + ext)
                if src.exists():
                    (mirror / src.name).write_bytes(src.read_bytes())
        # primary unsplit names
        base = ("table4a_gremln_evidence" if info["model"] == "GREmLN"
                else "table4b_genie3_evidence")
        for ext in (".svg", ".png"):
            src = out / f"{base}{ext}"
            if src.exists():
                (mirror / src.name).write_bytes(src.read_bytes())

    # write redesign report
    lines = [
        "# Table 4A/4B publication redesign report",
        "",
        "Presentation-only redesign. Evidence values unchanged from audited CSVs.",
        "",
        f"- Page size: {PAGE_W} × {PAGE_H} in landscape",
        f"- Title: {TITLE_PT} pt bold",
        f"- Headers: {HDR_PT} pt bold",
        f"- Body: {BODY_PT} pt (minimum enforced {MIN_BODY_PT} pt)",
        f"- Notes: {NOTE_PT} pt",
        "",
    ]
    for info in reports:
        lines += [
            f"## {info['model']}",
            f"- pages: {info['pages']}",
            f"- max body lines in any cell: {info['max_body_lines']}",
            f"- stems: {', '.join(info['page_stems'])}",
            "",
        ]
    lines += [
        "## Reconciliation (preserved)",
        "- GREmLN: CRISPRi 3; lit strong/moderate 7; mRNA 2; independent 4",
        "- GENIE3: CRISPRi 6; lit strong/moderate 9; mRNA 5; independent 3",
        "- Shared independent: 0",
        "",
        "All typography, wrap and scientific assertions passed.",
    ]
    (out / "table4_redesign_report.md").write_text("\n".join(lines))

    for info in reports:
        print(f"[{info['model']}] pages={info['pages']} body={info['body_pt']}pt "
              f"max_lines={info['max_body_lines']} -> {info['page_stems']}")
    print("wrote", out / "table4_redesign_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
