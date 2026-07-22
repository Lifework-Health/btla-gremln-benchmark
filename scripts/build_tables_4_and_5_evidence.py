#!/usr/bin/env python3
"""Update Tables 4A/4B and create Table 5 from final audited sources.

Sources (read only; no model rerun, no composite score):
  - Table 4A/4B Position order (canonical seed-excluded ordinals)
  - crispri_all_screen_sensitivity.csv (reconciled primary CRISPRi)
  - paperclip_pipeline_trace.csv final_usable_tier (Paperclip v2)
  - table_S7 / mrna_contrast_layers (multiomics audit_v2)

NR4A2: audited provisional Moderate (entity-corrected; pending human
verification) — pipeline_trace still records final_usable_tier=missing
after rubric failure; the audited provisional display tier is Moderate.

Does not edit manuscript prose.
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

DASH = "\u2014"
SUP_S = "\u02e2"
PG = "#E3F1E1"       # pale green — detected / qualifying independent / Strong lit
PYEL = "#FFF8DC"     # pale yellow — contextual mRNA evidence
GREY = "#EEEEEE"     # not represented / not measured / failed QC / None
HEADER_BG = "#f2f2f2"
ZEBRA = "#fafafa"
TEXT = "#222"

# Audited provisional override (documented entity correction)
TIER_OVERRIDE = {
    "NR4A2": "moderate",  # entity corrected; pending human verification
}


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def strip_shared(name: str) -> str:
    return re.sub(r"\u02e2$", "", str(name)).strip()


def fmt_tps(raw: str) -> str:
    """Normalise timepoint strings to '1 h, 4 h' style."""
    raw = str(raw).replace("min", " min").replace("h", " h")
    # collapse double spaces from already-spaced forms
    parts = [p.strip() for p in re.split(r"[;,]", raw) if p.strip()]
    out = []
    for p in parts:
        p = p.replace("  ", " ")
        # already like '2 min' / '4 h'
        if re.match(r"^\d+\s*(min|h)$", p):
            out.append(p)
        elif p.endswith("mn"):
            out.append(p.replace("mn", " min"))
        else:
            out.append(p)
    # unique preserve order
    seen = set()
    uniq = []
    for p in out:
        if p not in seen:
            seen.add(p)
            uniq.append(p)
    return ", ".join(uniq)


def crispri_label(row) -> str:
    if not bool(row["in_screen_8hr"]):
        return "Not represented in screen"
    if row["ontarget_qc_pass"] != True:
        return "Failed on-target quality control"
    npred = int(row["n_predicted_native"]) if pd.notna(row["n_predicted_native"]) else 0
    nconf = int(row["n_confirmed_native_8hr"]) if pd.notna(row["n_confirmed_native_8hr"]) else 0
    nd = f" ({nconf}/{npred})"
    d = str(row["response_direction"])
    return {
        "BTLA_concordant": f"Detected {DASH} BTLA-concordant" + nd,
        "anti_concordant": f"Detected {DASH} BTLA-anti-concordant" + nd,
        "mixed": f"Detected {DASH} mixed" + nd,
    }.get(d, "No detected response" + nd)


def lit_label(tier: str) -> str:
    t = str(tier).lower()
    if t in TIER_OVERRIDE.values() or t in ("strong", "moderate", "weak", "none"):
        return {"strong": "Strong", "moderate": "Moderate", "weak": "Weak", "none": "None"}[t]
    if t == "missing":
        return "None"  # should be overridden before call
    return "None"


def mrna_label(s7_row, mrna_row) -> str:
    st = mrna_row["TCR_activation__status"]
    if st == "measured_qualifying":
        tps = fmt_tps(mrna_row["TCR_activation__qual_timepoints"].replace(";", ", "))
        direction = str(mrna_row["TCR_activation__directions"])
        arrow = "↑" if direction == "up" else "↓" if direction == "down" else ""
        return f"TCR activation {arrow}: {tps}"
    if st == "absent_from_assay_universe":
        return "Absent from transcriptomic assay universe"
    if st == "failed_qc":
        return "Failed quality control"
    return "No qualifying transcriptional response"


def indep_label(s7_row) -> str:
    prot = str(s7_row["Protein abundance"])
    phos = str(s7_row["Phosphosite abundance"])
    indep = str(s7_row["Independent molecular association"])
    if indep.startswith("Yes"):
        bits = []
        if prot.startswith("Qualifying"):
            tps = prot.split(":", 1)[1].strip()
            bits.append(f"Protein association: {tps}")
        if phos.startswith("Qualifying"):
            # drop site-count parenthetical
            body = phos.split(":", 1)[1].strip()
            body = re.sub(r"\s*\(\d+ qualifying sites?\)", "", body).strip()
            bits.append(f"Phosphosite association: {body}")
        return "; ".join(bits) if bits else "No qualifying association in measured layers"
    # measured in at least one layer?
    measured = (
        prot.startswith("Measured") or prot.startswith("Qualifying") or
        phos.startswith("Measured") or phos.startswith("Qualifying")
    )
    if measured:
        return "No qualifying association in measured layers"
    return "Not measured in protein or phosphoproteomic assays"


def load_canonical_order(repo: Path, model: str) -> list[str]:
    """Prefer the frozen previous-publication Table 4 order; fall back to current."""
    stem = ("table4a_gremln_top25_evidence" if model == "GREmLN"
            else "table4b_genie3_top25_evidence")
    prev = repo / "results/publication_assets/tables" / f"{stem}.previous_publication.csv"
    cur = repo / "results/publication_assets/tables" / f"{stem}.csv"
    path = prev if prev.exists() else cur
    t4 = pd.read_csv(path)
    assert list(t4["Position"]) == list(range(1, 26))
    return [strip_shared(x) for x in t4["Candidate regulator"]]


def build_model_table(model: str, order: list[str], shared: set,
                      cr: pd.DataFrame, tiers: dict, s7: pd.DataFrame, mrna: pd.DataFrame):
    s7i = s7.set_index("Candidate regulator")
    mi = mrna.set_index("TF")
    crx = cr[cr["model"] == model].set_index("TF")
    rows = []
    for i, tf in enumerate(order):
        cand = f"{tf}{SUP_S}" if tf in shared else tf
        if tf not in crx.index:
            cris = "Not represented in screen"
        else:
            cris = crispri_label(crx.loc[tf])
        tier = tiers.get(tf, "none")
        if tf in TIER_OVERRIDE:
            tier = TIER_OVERRIDE[tf]
        lit = lit_label(tier)
        mrna_cell = mrna_label(s7i.loc[tf], mi.loc[tf])
        indep = indep_label(s7i.loc[tf])
        rows.append({
            "Position": i + 1,
            "Candidate regulator": cand,
            "CRISPRi evidence": cris,
            "Audited literature evidence": lit,
            "Candidate mRNA context": mrna_cell,
            "Independent BTLA molecular evidence": indep,
            "_TF": tf,
            "_tier": tier,
        })
    return pd.DataFrame(rows)


def cell_styles(df: pd.DataFrame):
    """Return {(row, col): bg} for manuscript-style shading."""
    cols = list(df.columns)
    cbg = {}
    for i in range(len(df)):
        # CRISPRi
        j = cols.index("CRISPRi evidence")
        v = df.iloc[i, j]
        if str(v).startswith("Detected"):
            cbg[(i, j)] = PG
        elif v in ("Not represented in screen", "Failed on-target quality control"):
            cbg[(i, j)] = GREY
        # Literature
        j = cols.index("Audited literature evidence")
        v = df.iloc[i, j]
        if v == "Strong":
            cbg[(i, j)] = PG
        elif v == "None":
            cbg[(i, j)] = GREY
        # mRNA context — pale yellow for qualifying contextual
        j = cols.index("Candidate mRNA context")
        v = df.iloc[i, j]
        if str(v).startswith("TCR activation"):
            cbg[(i, j)] = PYEL
        elif v == "Absent from transcriptomic assay universe" or v == "Failed quality control":
            cbg[(i, j)] = GREY
        # Independent
        j = cols.index("Independent BTLA molecular evidence")
        v = df.iloc[i, j]
        if "association:" in str(v):
            cbg[(i, j)] = PG
        elif v == "Not measured in protein or phosphoproteomic assays":
            cbg[(i, j)] = GREY
    return cbg


def render_table(df: pd.DataFrame, title: str, footnotes: list[str],
                 path_stem: Path, col_chars: list[int], cbg: dict):
    """Booktabs-style publication table matching existing Table 4 look."""
    cols = [c for c in df.columns if not c.startswith("_")]
    body_df = df[cols]
    wrapw = col_chars
    tot = sum(wrapw)
    frac = [w / tot for w in wrapw]
    xb = np.concatenate([[0.0], np.cumsum(frac)])

    hdr = [textwrap.fill(c, wrapw[j]) for j, c in enumerate(cols)]
    body = [[textwrap.fill(str(body_df.iloc[i][cols[j]]), wrapw[j])
             for j in range(len(cols))] for i in range(len(body_df))]
    fnw = [textwrap.fill(fn, max(95, int(tot * 1.05))) for fn in footnotes]

    def nlines(row):
        return max(1, max(s.count("\n") + 1 for s in row))

    fs = 7.4
    hL = nlines(hdr)
    rL = [nlines(r) for r in body]
    fL = sum(f.count("\n") + 1 for f in fnw)
    top = 2.2
    y0 = top + hL
    y_end = y0 + sum(rL)
    total = y_end + 1.2 + fL * 0.95 + 0.5

    figw = min(18.5, max(11.5, tot * 0.125))
    figh = max(9.5, total * 0.27)
    fig, ax = plt.subplots(figsize=(figw, figh))
    ax.set_xlim(0, 1)
    ax.set_ylim(total, 0)
    ax.axis("off")

    ax.text(0.0, top * 0.40, title, fontsize=fs + 3, fontweight="bold", va="center")
    ax.add_patch(Rectangle((0, top), 1, hL, fc=HEADER_BG, ec="none"))
    for j in range(len(cols)):
        ax.text(xb[j] + 0.003, top + hL / 2, hdr[j], fontsize=fs,
                fontweight="bold", va="center")
    for yy, lw in ((top, 1.35), (y0, 0.95), (y_end, 1.35)):
        ax.plot([0, 1], [yy, yy], color="#222", lw=lw, solid_capstyle="butt")

    cy = y0
    for i in range(len(body)):
        h = rL[i]
        if i % 2 == 1:
            ax.add_patch(Rectangle((0, cy), 1, h, fc=ZEBRA, ec="none"))
        for j in range(len(cols)):
            if (i, j) in cbg:
                ax.add_patch(Rectangle((xb[j], cy), frac[j], h, fc=cbg[(i, j)], ec="none"))
        for j in range(len(cols)):
            cell_val = str(body_df.iloc[i][cols[j]])
            bold = False
            if cols[j] == "Candidate regulator" and SUP_S in cell_val:
                bold = True
            if (cell_val.startswith("Detected") or cell_val == "Strong" or
                    "association:" in cell_val or cell_val.startswith("TCR activation")):
                bold = True
            ax.text(xb[j] + 0.003, cy + h / 2, body[i][j], fontsize=fs - 0.3, va="center",
                    fontweight="bold" if bold else "normal")
        cy += h

    fy = y_end + 0.95
    for f in fnw:
        ax.text(0.0, fy, f, fontsize=fs - 1.6, va="top", color="#333")
        fy += (f.count("\n") + 1) * 0.95

    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(Path(str(path_stem) + ext), dpi=400 if ext == ".png" else 300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


FOOT_SHARED = [
    "S, shared between the GREmLN and GENIE3 top-25 lists.",
    "Candidate regulators were restricted to the pre-specified pySCENIC human transcription-factor/regulator list.",
    "CRISPRi evidence is model-specific because each model can nominate a different target set for the same candidate. "
    "The (confirmed/predicted) denominator uses the native predicted-target budget (not a matched top-N budget). "
    "Not represented in screen and failed on-target quality control are shown explicitly and are not converted to zero response.",
    "Audited literature evidence uses the Paperclip v2 TCR-inhibition audit (fixed query \"<TF> TCR inhibition\"). "
    "Strong and moderate literature tiers are audited provisional annotations pending documented human verification.",
    "Candidate mRNA context is contextual evidence only (TCR versus unstimulated control). "
    "No candidate showed a qualifying incremental BTLA+TCR versus TCR mRNA response; that contrast is therefore omitted from the column.",
    "Independent BTLA molecular evidence is restricted to directly measured qualifying protein abundance or phosphosite change "
    "(BTLA+TCR versus TCR) at the source-native thresholds. Inferred regulon activity, kinase activity, BIONIC, co-IP and "
    "synapse annotations are not shown. Molecular direction is not interpreted as supportive or opposing because both "
    "model rankings are unsigned. Incomplete assay coverage is reported as not measured and is not negative evidence.",
    "Pale green marks detected CRISPRi response, Strong literature, or qualifying independent molecular association. "
    "Pale yellow marks qualifying contextual mRNA. Light grey marks not represented, failed QC, not measured, or None.",
]


def compare_cells(old: pd.DataFrame, new: pd.DataFrame, model: str) -> list[str]:
    """Report every cell that changed from the previous publication Table 4."""
    changes = []
    old = old.copy()
    old["_TF"] = old["Candidate regulator"].map(strip_shared)
    # map old columns to new
    colmap = {
        "CRISPRi evidence": "CRISPRi evidence",
        "Literature evidence": "Audited literature evidence",
        "BTLA experimental evidence": "Independent BTLA molecular evidence",
    }
    # also mRNA is new
    for _, r in new.iterrows():
        tf = r["_TF"]
        o = old[old["_TF"] == tf]
        if o.empty:
            changes.append(f"{model} {tf}: NEW ROW")
            continue
        o = o.iloc[0]
        # CRISPRi — normalise old labels for comparison of substance
        old_c = str(o["CRISPRi evidence"])
        new_c = str(r["CRISPRi evidence"])
        if old_c != new_c:
            changes.append(f"{model} {tf} | CRISPRi: '{old_c}' → '{new_c}'")
        old_l = str(o["Literature evidence"])
        new_l = str(r["Audited literature evidence"])
        if old_l != new_l:
            changes.append(f"{model} {tf} | Literature: '{old_l}' → '{new_l}'")
        # experimental → split into mRNA + independent (always a structural change)
        old_e = str(o["BTLA experimental evidence"])
        new_m = str(r["Candidate mRNA context"])
        new_i = str(r["Independent BTLA molecular evidence"])
        changes.append(
            f"{model} {tf} | BTLA experimental evidence '{old_e}' → "
            f"mRNA context '{new_m}' | independent '{new_i}'"
        )
    return changes


def main() -> int:
    repo = repo_root()
    out = repo / "results/publication_assets/tables"
    out.mkdir(parents=True, exist_ok=True)

    cr = pd.read_csv(repo / "results/publication_data/crispri_all_screen_sensitivity.csv")
    trace = pd.read_csv(repo / "results/paperclip/v2_tcr_inhibition/paperclip_pipeline_trace.csv")
    tiers = dict(zip(trace["TF"], trace["final_usable_tier"].astype(str).str.lower()))
    s7 = pd.read_csv(repo / "results/multiomics/audit_v2/table_S7_audited_multiomics_43_candidates.csv")
    mrna = pd.read_csv(repo / "results/multiomics/audit_v2/mrna_contrast_layers_summary_v2.csv")
    union = pd.read_csv(repo / "results/publication_data/top25_union_primary.csv")
    shared = set(union.loc[union["in_gremln_top25"] & union["in_genie3_top25"], "TF"])

    old_a_path = out / "table4a_gremln_top25_evidence.previous_publication.csv"
    old_b_path = out / "table4b_genie3_top25_evidence.previous_publication.csv"
    old_a = pd.read_csv(old_a_path if old_a_path.exists() else out / "table4a_gremln_top25_evidence.csv",
                        keep_default_na=False)
    old_b = pd.read_csv(old_b_path if old_b_path.exists() else out / "table4b_genie3_top25_evidence.csv",
                        keep_default_na=False)
    # If current file already has the new schema and no previous snapshot, skip cell-diff
    if "Literature evidence" not in old_a.columns:
        raise SystemExit("Missing previous_publication snapshot for Table 4 cell-diff")

    order_a = load_canonical_order(repo, "GREmLN")
    order_b = load_canonical_order(repo, "GENIE3")
    assert order_a == [strip_shared(x) for x in old_a["Candidate regulator"]]
    assert order_b == [strip_shared(x) for x in old_b["Candidate regulator"]]

    t4a = build_model_table("GREmLN", order_a, shared, cr, tiers, s7, mrna)
    t4b = build_model_table("GENIE3", order_b, shared, cr, tiers, s7, mrna)

    # ---------- assertions ----------
    A = []

    def check(name, cond):
        A.append((name, bool(cond)))
        if not cond:
            print("FAIL", name)

    check("t4a_25", len(t4a) == 25)
    check("t4b_25", len(t4b) == 25)
    check("t4a_order", list(t4a["_TF"]) == order_a)
    check("t4b_order", list(t4b["_TF"]) == order_b)

    # shared annotations identical
    for tf in shared:
        ra = t4a[t4a["_TF"] == tf].iloc[0]
        rb = t4b[t4b["_TF"] == tf].iloc[0]
        check(f"shared_lit_{tf}", ra["Audited literature evidence"] == rb["Audited literature evidence"])
        check(f"shared_mrna_{tf}", ra["Candidate mRNA context"] == rb["Candidate mRNA context"])
        check(f"shared_indep_{tf}",
              ra["Independent BTLA molecular evidence"] == rb["Independent BTLA molecular evidence"])

    def n_detected(df, model):
        return int(df["CRISPRi evidence"].str.startswith("Detected").sum())

    check("crispri_3_vs_6", n_detected(t4a, "GREmLN") == 3 and n_detected(t4b, "GENIE3") == 6)

    def lit_counts(df):
        vc = df["Audited literature evidence"].value_counts()
        return {k: int(vc.get(k, 0)) for k in ("Strong", "Moderate", "Weak", "None")}

    la, lb = lit_counts(t4a), lit_counts(t4b)
    check("lit_gremln", la == {"Strong": 1, "Moderate": 6, "Weak": 3, "None": 15})
    check("lit_genie3", lb == {"Strong": 2, "Moderate": 7, "Weak": 6, "None": 10})
    check("lit_sm_7_9",
          la["Strong"] + la["Moderate"] == 7 and lb["Strong"] + lb["Moderate"] == 9)

    def n_mrna(df):
        return int(df["Candidate mRNA context"].str.startswith("TCR activation").sum())

    check("mrna_2_vs_5", n_mrna(t4a) == 2 and n_mrna(t4b) == 5)

    def n_indep(df):
        return int(df["Independent BTLA molecular evidence"].str.contains("association:").sum())

    check("indep_4_vs_3", n_indep(t4a) == 4 and n_indep(t4b) == 3)

    shared_indep = [
        tf for tf in shared
        if "association:" in str(t4a[t4a["_TF"] == tf]["Independent BTLA molecular evidence"].iloc[0])
    ]
    check("no_shared_indep", shared_indep == [])

    # no legacy annotations
    blob = "\n".join(
        t4a["Independent BTLA molecular evidence"].astype(str).tolist() +
        t4b["Independent BTLA molecular evidence"].astype(str).tolist() +
        t4a["Candidate mRNA context"].astype(str).tolist() +
        t4b["Candidate mRNA context"].astype(str).tolist()
    )
    for bad in ["Inferred regulon", "BIONIC", "Co-IP", "coIP", "kinase", "synapse", "regulon"]:
        check(f"no_legacy_{bad}", bad.lower() not in blob.lower())

    # no restricted values
    check("no_residues", not re.search(r"[A-Za-z0-9]+_[STY]\d+", blob))
    check("no_redacted", "REDACTED" not in blob)

    failed = [n for n, ok in A if not ok]
    if failed:
        raise SystemExit(f"ASSERTION FAILURE: {failed}")

    # ---------- Table 5 ----------
    # CRISPRi denominators
    def crispri_stats(df, model):
        tfs = set(df["_TF"])
        sub = cr[(cr["model"] == model) & (cr["TF"].isin(tfs))]
        represented = int(sub["in_screen_8hr"].sum())
        evaluable = int(((sub["in_screen_8hr"]) & (sub["ontarget_qc_pass"] == True)).sum())
        detected = n_detected(df, model)
        return detected, represented, evaluable

    ga_d, ga_r, ga_e = crispri_stats(t4a, "GREmLN")
    gb_d, gb_r, gb_e = crispri_stats(t4b, "GENIE3")

    t5 = pd.DataFrame([
        {
            "Evidence source": "CRISPRi perturbational evidence",
            "GREmLN": f"{ga_d} candidates with a detected response among {ga_r} represented",
            "GENIE3": f"{gb_d} candidates with a detected response among {gb_e} evaluable after QC",
            "Interpretation": "Broader perturbational evidence for GENIE3 in the related CD4 T cell screen",
        },
        {
            "Evidence source": "Audited literature evidence",
            "GREmLN": f"{la['Strong'] + la['Moderate']}/25 strong or moderate",
            "GENIE3": f"{lb['Strong'] + lb['Moderate']}/25 strong or moderate",
            "Interpretation": "Modest literature lean towards GENIE3; tiers remain provisional",
        },
        {
            "Evidence source": "TCR activation mRNA context",
            "GREmLN": f"{n_mrna(t4a)}/25",
            "GENIE3": f"{n_mrna(t4b)}/25",
            "Interpretation": ("Activation-associated context was more frequent for GENIE3; "
                               "no candidate showed a qualifying incremental BTLA-specific mRNA response"),
        },
        {
            "Evidence source": "Independent BTLA molecular association",
            "GREmLN": f"{n_indep(t4a)}/25",
            "GENIE3": f"{n_indep(t4b)}/25",
            "Interpretation": ("Slight numerical lean towards GREmLN, but protein and phosphosite "
                               "coverage was sparse and unequal"),
        },
        {
            "Evidence source": "Overall",
            "GREmLN": "Distinct candidates and molecular hypotheses",
            "GENIE3": "Broader perturbational and literature corroboration",
            "Interpretation": "No clear or consistent superiority across evidence sources",
        },
    ])

    # Table 5 must reconcile with 4A/4B
    assert f"{ga_d} candidates" in t5.iloc[0]["GREmLN"]
    assert f"{gb_d} candidates" in t5.iloc[0]["GENIE3"]
    assert t5.iloc[1]["GREmLN"].startswith("7/25")
    assert t5.iloc[1]["GENIE3"].startswith("9/25")
    assert t5.iloc[2]["GREmLN"] == "2/25" and t5.iloc[2]["GENIE3"] == "5/25"
    assert t5.iloc[3]["GREmLN"] == "4/25" and t5.iloc[3]["GENIE3"] == "3/25"

    # ---------- write CSVs ----------
    pub_cols = [c for c in t4a.columns if not c.startswith("_")]
    t4a[pub_cols].to_csv(out / "table4a_gremln_evidence.csv", index=False)
    t4b[pub_cols].to_csv(out / "table4b_genie3_evidence.csv", index=False)
    t5.to_csv(out / "table5_triangulated_evidence_summary.csv", index=False)
    # also refresh legacy filenames used by the notebook for continuity
    t4a[pub_cols].to_csv(out / "table4a_gremln_top25_evidence.csv", index=False)
    t4b[pub_cols].to_csv(out / "table4b_genie3_top25_evidence.csv", index=False)

    # ---------- render ----------
    COLW = [6, 12, 28, 12, 22, 28]
    render_table(t4a, "Table 4A. Evidence across GREmLN's top-25 candidate regulators",
                 FOOT_SHARED, out / "table4a_gremln_evidence", COLW, cell_styles(t4a))
    render_table(t4b, "Table 4B. Evidence across GENIE3's top-25 candidate regulators",
                 FOOT_SHARED, out / "table4b_genie3_evidence", COLW, cell_styles(t4b))

    # Table 5 — simpler wider text columns
    t5_cbg = {}
    render_table(t5, "Table 5. Triangulated evidence across the two top-25 candidate lists",
                 [
                     "Descriptive comparison of evidence coverage across each model's top-25; "
                     "no weighted score, Spearman correlation, bootstrap interval or superiority margin is implied.",
                     "CRISPRi: GREmLN denominator is candidates represented in the screen; "
                     "GENIE3 denominator is candidates evaluable after on-target QC.",
                     "Literature tiers are audited provisional Paperclip v2 annotations pending documented human verification.",
                     "Candidate mRNA is contextual; independent molecular association is restricted to qualifying "
                     "protein or phosphosite measurements. Assay coverage is incomplete and unequal across candidates.",
                 ],
                 out / "table5_triangulated_evidence_summary",
                 [22, 28, 30, 40], t5_cbg)

    # also write svg/png under the legacy stem names for notebook continuity
    for src, dst in [
        ("table4a_gremln_evidence", "table4a_gremln_top25_evidence"),
        ("table4b_genie3_evidence", "table4b_genie3_top25_evidence"),
    ]:
        for ext in (".png", ".svg", ".pdf"):
            s, d = out / f"{src}{ext}", out / f"{dst}{ext}"
            if s.exists():
                d.write_bytes(s.read_bytes())

    # ---------- reconciliation report ----------
    changes = compare_cells(old_a, t4a, "4A") + compare_cells(old_b, t4b, "4B")
    lines = [
        "# Table evidence reconciliation",
        "",
        "Updated Tables 4A/4B and new Table 5 from final audited sources "
        "(CRISPRi primary, Paperclip v2, multiomics audit_v2). "
        "Manuscript prose was not edited.",
        "",
        "## Headline counts",
        "",
        f"- CRISPRi detected: GREmLN {ga_d}/{ga_r} represented; GENIE3 {gb_d}/{gb_e} evaluable after QC",
        f"- Literature strong/moderate: GREmLN {la['Strong']+la['Moderate']}/25 "
        f"({la}); GENIE3 {lb['Strong']+lb['Moderate']}/25 ({lb})",
        f"- TCR activation mRNA: GREmLN {n_mrna(t4a)}/25; GENIE3 {n_mrna(t4b)}/25",
        f"- Independent molecular: GREmLN {n_indep(t4a)}/25; GENIE3 {n_indep(t4b)}/25; shared 0",
        f"- NR4A2 literature display tier: Moderate "
        f"(audited provisional; pipeline_trace final_usable_tier was missing; "
        f"entity corrected, pending human verification)",
        "",
        f"## Cell changes from previous publication Table 4 ({len(changes)} entries)",
        "",
    ]
    for ch in changes:
        lines.append(f"- {ch}")
    lines += [
        "",
        "## Assertions",
        "",
        f"All {len(A)} executable checks passed.",
        "",
        "## Outputs",
        "",
        "- `table4a_gremln_evidence.csv/.png/.svg`",
        "- `table4b_genie3_evidence.csv/.png/.svg`",
        "- `table5_triangulated_evidence_summary.csv/.png/.svg`",
        "- `table_evidence_reconciliation.md`",
    ]
    (out / "table_evidence_reconciliation.md").write_text("\n".join(lines))

    print(f"[Tables 4A/4B/5] OK  lit={la}/{lb} mrna={n_mrna(t4a)}/{n_mrna(t4b)} "
          f"indep={n_indep(t4a)}/{n_indep(t4b)} crispri={ga_d}/{gb_d}")
    print(f"cell changes logged: {len(changes)}")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
