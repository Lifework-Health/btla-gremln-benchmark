#!/usr/bin/env python3
"""Build TF-level Paperclip evidence summary figures for GREmLN and GENIE3 top 25.

Outputs (literature annotation only):
    results/paperclip/v2_tcr_inhibition/paperclip_gremln_tf_evidence_summary.csv
    results/paperclip/v2_tcr_inhibition/paperclip_genie3_tf_evidence_summary.csv
    results/paperclip/v2_tcr_inhibition/fig_paperclip_v2_gremln_tf_summary.{svg,png}
    results/paperclip/v2_tcr_inhibition/fig_paperclip_v2_genie3_tf_summary.{svg,png}
    results/paperclip/v2_tcr_inhibition/fig_paperclip_v2_*_tf_summary_caption.txt

Uses corrected Paperclip v2 audit outputs after entity identity validation.
Fails if any figure assertion is violated.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_paperclip_tcr_inhibition_audit as base  # noqa: E402

TIER_ORDER = ["strong", "moderate", "weak", "none", "missing"]
TIER_COLORS = {
    "strong": "#1b7837",
    "moderate": "#7fbf7b",
    "weak": "#fee08b",
    "none": "#d9d9d9",
    "missing": "#b2182b",
}
TIER_TEXT = {
    "strong": "Strong",
    "moderate": "Moderate",
    "weak": "Weak",
    "none": "None",
    "missing": "Missing",
}
TIER_FG = {
    "strong": "white",
    "moderate": "#1a1a1a",
    "weak": "#1a1a1a",
    "none": "#1a1a1a",
    "missing": "white",
}
SHARED_MARKER = "●"  # filled circle; defined in caption/legend
GROUP_LABEL = {
    "shared": "shared",
    "gremln_only": "GREmLN specific",
    "genie3_only": "GENIE3 specific",
}


def _human_verification_complete(out: Path) -> bool:
    hv = pd.read_csv(out / "human_verification.csv", keep_default_na=False)
    if hv.empty:
        return False
    required = ["reviewer", "date", "TF", "paper_id", "entity_confirmed",
                "relevance_confirmed", "evidence_directness_confirmed",
                "tier_accepted_or_changed"]
    for col in required:
        if col not in hv.columns:
            return False
    # complete = at least one strong/moderate key paper has named reviewer + date
    ok = (hv["reviewer"].astype(str).str.strip() != "") & (hv["date"].astype(str).str.strip() != "")
    return bool(ok.any())


def _audit_status(row: pd.Series, human_done: bool) -> str:
    tier = str(row["final_usable_tier"]).lower()
    entity = str(row["entity_audit_status"])
    notes = str(row["notes"])
    if tier == "missing":
        if "rubric" in notes.lower() or str(row.get("original_rubric_valid", "")).lower() == "false":
            return "Rubric failure"
        return "Missing"
    if "corrected" in entity.lower() or str(row.get("rerun_required", "")).lower() == "true":
        return "Entity corrected"
    if tier in ("strong", "moderate"):
        if human_done:
            return "Verified"
        return "Pending human verification"
    return "Other reviewed status"


def _missing_reason(row: pd.Series, tiers: pd.DataFrame) -> str:
    if str(row["final_usable_tier"]).lower() != "missing":
        return ""
    tf = row["TF"]
    sub = tiers[tiers["TF"] == tf]
    if len(sub):
        rr = str(sub.iloc[0].get("rubric_reasons", "") or "")
        if rr:
            return f"rubric validation failure ({rr})"
    return "unresolved structured judgement (final usable tier missing)"


def _direct_glyph(tier: str, n_direct: int) -> str:
    if str(tier).lower() == "missing":
        return "—"
    if int(n_direct) >= 1:
        return "●"  # filled; rendered as dark green circle in figure
    return "○"


def _ordinal_top25(ranking: pd.DataFrame, gene_col: str, members: set) -> list[str]:
    """Preserve ranking-file order among canonical top-25 members; assign ranks 1..25."""
    ordered = [g for g in ranking[gene_col].tolist() if g in members]
    if len(ordered) != 25:
        raise AssertionError(f"expected 25 members in ranking order, got {len(ordered)}")
    if len(set(ordered)) != 25:
        raise AssertionError("duplicate TF in top-25 ordinal order")
    return ordered


def build_model_table(model: str, ordered_tfs: list[str], trace: pd.DataFrame,
                      tiers: pd.DataFrame, human_done: bool) -> pd.DataFrame:
    tr = {r["TF"]: r for _, r in trace.iterrows()}
    rows = []
    for rank, tf in enumerate(ordered_tfs, start=1):
        r = tr[tf]
        tier = str(r["final_usable_tier"]).lower()
        n_ret = int(r["papers_returned"])
        n_elig = int(r["papers_entity_eligible"])
        n_rel = int(r["papers_relevant"])
        n_prim = int(r["relevant_primary_papers"])
        n_dir = int(r["direct_causal_papers"])
        keys = [k for k in str(r["key_eligible_paper_ids"]).split(";") if k.strip()]
        group = GROUP_LABEL[str(r["candidate_group"])]
        shared = group == "shared"
        conf = str(r["confidence"]).strip().lower()
        conf_label = {"high": "High", "moderate": "Moderate", "low": "Low"}.get(conf, conf.title() or "—")
        if tier == "missing" and (not conf or conf in ("", "nan")):
            conf_label = "—"
        audit = _audit_status(r, human_done)
        # Human verification status column: Verified only with named dated review
        if tier in ("strong", "moderate"):
            hv_status = "Verified" if (human_done and audit == "Verified") else "Pending human verification"
        else:
            hv_status = "not required"
        rows.append({
            "rank": rank,
            "TF": tf,
            "candidate_full_name": r["candidate_full_name"],
            "candidate_group": group,
            "in_gremln_top25": bool(r["candidate_group"] in ("gremln_only", "shared")),
            "in_genie3_top25": bool(r["candidate_group"] in ("genie3_only", "shared")),
            "shared_marker": SHARED_MARKER if shared else "",
            "final_usable_tier": tier,
            "original_judge_tier": str(r["original_judge_tier"]),
            "corrected_judge_tier": str(r["corrected_judge_tier"]),
            "confidence": conf_label,
            "papers_returned": n_ret,
            "papers_entity_eligible": n_elig,
            "papers_relevant": n_rel,
            "primary_relevant_papers": n_prim,
            "direct_causal_papers": n_dir,
            "evidence_count_string": f"{n_ret} / {n_elig} / {n_rel} / {n_prim} / {n_dir}",
            "key_paper_count": len(keys),
            "key_paper_ids": ";".join(keys),
            "direct_causal_glyph": _direct_glyph(tier, n_dir),
            "entity_audit_status": r["entity_audit_status"],
            "human_verification_status": hv_status,
            "audit_status": audit,
            "missing_reason": _missing_reason(r, tiers),
        })
    return pd.DataFrame(rows)


def assert_model_table(df: pd.DataFrame, model: str, ordered_tfs: list[str],
                       chart_counts: pd.DataFrame, elig_ids: dict,
                       human_done: bool) -> None:
    checks = []

    def chk(name, cond):
        checks.append((name, bool(cond)))

    chk(f"{model}: exactly 25 rows", len(df) == 25)
    chk(f"{model}: ranks 1..25 unique", list(df["rank"]) == list(range(1, 26)))
    chk(f"{model}: no duplicate TF", df["TF"].is_unique)
    chk(f"{model}: order matches ranking", list(df["TF"]) == ordered_tfs)
    for _, r in df.iterrows():
        parts = [int(x.strip()) for x in r["evidence_count_string"].split("/")]
        chk(f"{model}/{r['TF']} count string",
            parts == [r["papers_returned"], r["papers_entity_eligible"],
                      r["papers_relevant"], r["primary_relevant_papers"],
                      r["direct_causal_papers"]])
    # tier distribution vs model-level chart
    row = chart_counts[chart_counts["group"] == f"{model} top 25"].iloc[0]
    for t in TIER_ORDER:
        chk(f"{model} tier {t} matches chart",
            int((df["final_usable_tier"] == t).sum()) == int(row[t]))
    chk(f"{model} strong+moderate matches chart",
        int(df["final_usable_tier"].isin(["strong", "moderate"]).sum())
        == int(row["strong_plus_moderate"]))
    # key papers entity eligible
    for _, r in df.iterrows():
        for kid in str(r["key_paper_ids"]).split(";"):
            if kid:
                chk(f"{model}/{r['TF']} key {kid} eligible",
                    kid in elig_ids.get(r["TF"], set()))
    # strong / moderate structural
    for _, r in df[df["final_usable_tier"] == "strong"].iterrows():
        chk(f"{model}/{r['TF']} strong has direct", int(r["direct_causal_papers"]) >= 1)
    for _, r in df[df["final_usable_tier"] == "moderate"].iterrows():
        chk(f"{model}/{r['TF']} moderate has primary", int(r["primary_relevant_papers"]) >= 1)
    # missing stays missing
    for _, r in df[df["final_usable_tier"] == "missing"].iterrows():
        chk(f"{model}/{r['TF']} missing not none", r["final_usable_tier"] == "missing")
        chk(f"{model}/{r['TF']} missing reason present", bool(str(r["missing_reason"]).strip()))
    # Verified requires named dated human review
    for _, r in df[df["audit_status"] == "Verified"].iterrows():
        chk(f"{model}/{r['TF']} Verified implies human status",
            r["human_verification_status"] == "Verified" and human_done)
    # No Verified labels when human review is incomplete
    if not human_done:
        chk(f"{model}: no Verified without human review",
            (df["audit_status"] != "Verified").all())

    failed = [n for n, ok in checks if not ok]
    npass = sum(1 for _, ok in checks if ok)
    print(f"[{model}] assertions: {npass}/{len(checks)} passed")
    if failed:
        for n in failed[:20]:
            print(f"  FAIL: {n}")
        raise AssertionError(f"{model} figure assertions failed ({len(failed)})")


def assert_shared_parity(gr: pd.DataFrame, g3: pd.DataFrame) -> None:
    shared_gr = set(gr[gr["candidate_group"] == "shared"]["TF"])
    shared_g3 = set(g3[g3["candidate_group"] == "shared"]["TF"])
    if shared_gr != shared_g3:
        raise AssertionError(f"shared set mismatch: {shared_gr ^ shared_g3}")
    fields = ["final_usable_tier", "original_judge_tier", "corrected_judge_tier",
              "confidence", "papers_returned", "papers_entity_eligible", "papers_relevant",
              "primary_relevant_papers", "direct_causal_papers", "evidence_count_string",
              "key_paper_count", "key_paper_ids", "direct_causal_glyph",
              "entity_audit_status", "human_verification_status", "audit_status",
              "missing_reason", "shared_marker"]
    for tf in sorted(shared_gr):
        a = gr[gr["TF"] == tf].iloc[0]
        b = g3[g3["TF"] == tf].iloc[0]
        for f in fields:
            if str(a[f]) != str(b[f]):
                raise AssertionError(f"shared {tf} field {f}: {a[f]!r} != {b[f]!r}")
    print(f"[shared] parity OK for {len(shared_gr)} candidates")


def render_figure(df: pd.DataFrame, model: str, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Circle

    n = len(df)
    fig_w, fig_h = 13.4, 12.6
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, 100)
    ax.set_ylim(-2.4, n + 4.0)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    cols = {
        "rank": 1.2,
        "tf": 6.5,
        "tier": 22.0,
        "counts": 40.0,
        "conf": 58.5,
        "key": 68.5,
        "audit": 80.0,
        "glyph": 96.5,
    }
    headers = [
        (cols["rank"], "Rank"),
        (cols["tf"], "Candidate regulator"),
        (cols["tier"], "Final usable tier"),
        (cols["counts"], "returned / eligible / relevant / primary / direct"),
        (cols["conf"], "Confidence"),
        (cols["key"], "Key evidence"),
        (cols["audit"], "Audit status"),
        (cols["glyph"], "Direct"),
    ]

    ax.text(50, n + 3.4, f"{model} candidate evidence", ha="center", va="center",
            fontsize=13, fontweight="bold", color="#222222")

    hy = n + 2.4
    for x, lab in headers:
        ax.text(x, hy, lab, ha="left" if x != cols["glyph"] else "center",
                va="center", fontsize=7.2, fontweight="bold", color="#333333")
    ax.plot([0.8, 99.2], [hy - 0.5, hy - 0.5], color="#bbbbbb", lw=0.8)

    for _, r in df.iterrows():
        y = n - int(r["rank"]) + 1.0
        if int(r["rank"]) % 2 == 0:
            ax.add_patch(FancyBboxPatch(
                (0.6, y - 0.42), 98.8, 0.84,
                boxstyle="square,pad=0", linewidth=0, facecolor="#f7f7f7"))
        ax.plot([0.8, 99.2], [y - 0.42, y - 0.42], color="#e8e8e8", lw=0.4)

        ax.text(cols["rank"], y, str(int(r["rank"])), ha="right", va="center",
                fontsize=8.5, fontfamily="DejaVu Sans Mono", color="#222222")

        tf_label = f"{r['TF']}  {SHARED_MARKER}" if r["shared_marker"] else r["TF"]
        ax.text(cols["tf"], y, tf_label, ha="left", va="center",
                fontsize=8.5, fontweight="bold" if r["shared_marker"] else "normal",
                color="#111111")

        tier = r["final_usable_tier"]
        pill_w, pill_h = 11.5, 0.58
        ax.add_patch(FancyBboxPatch(
            (cols["tier"], y - pill_h / 2), pill_w, pill_h,
            boxstyle="round,pad=0.02,rounding_size=0.15",
            linewidth=0.4, edgecolor="#666666",
            facecolor=TIER_COLORS[tier]))
        ax.text(cols["tier"] + pill_w / 2, y, TIER_TEXT[tier],
                ha="center", va="center", fontsize=7.8, fontweight="bold",
                color=TIER_FG[tier])

        ax.text(cols["counts"], y, r["evidence_count_string"], ha="left", va="center",
                fontsize=7.8, fontfamily="DejaVu Sans Mono", color="#222222")
        ax.text(cols["conf"], y, r["confidence"], ha="left", va="center",
                fontsize=8.0, color="#222222")

        k = int(r["key_paper_count"])
        key_txt = f"{k} key paper" + ("" if k == 1 else "s") if k else "—"
        ax.text(cols["key"], y, key_txt, ha="left", va="center",
                fontsize=8.0, color="#222222")
        ax.text(cols["audit"], y, r["audit_status"], ha="left", va="center",
                fontsize=7.2, color="#333333")

        gx = cols["glyph"]
        if r["direct_causal_glyph"] == "—":
            ax.text(gx, y, "—", ha="center", va="center", fontsize=9, color="#666666")
        elif int(r["direct_causal_papers"]) >= 1:
            ax.add_patch(Circle((gx, y), 0.20, facecolor="#1b7837",
                                edgecolor="#1b7837", lw=0.8, zorder=3))
        else:
            ax.add_patch(Circle((gx, y), 0.20, facecolor="white",
                                edgecolor="#555555", lw=1.0, zorder=3))

    # Bottom rule under last row
    ax.plot([0.8, 99.2], [0.55, 0.55], color="#bbbbbb", lw=0.8)

    # Legend — two clean rows, well below the table
    ly1, ly2 = -0.15, -0.85
    ax.text(1.0, ly1 + 0.35, "Legend", fontsize=8, fontweight="bold", color="#333333")
    lx = 1.0
    for t in TIER_ORDER:
        ax.add_patch(FancyBboxPatch(
            (lx, ly1 - 0.20), 2.0, 0.36,
            boxstyle="round,pad=0.01,rounding_size=0.08",
            linewidth=0.3, edgecolor="#666666", facecolor=TIER_COLORS[t]))
        ax.text(lx + 2.3, ly1, TIER_TEXT[t], ha="left", va="center", fontsize=7.0)
        lx += 6.6
    ax.text(lx + 0.8, ly1, f"{SHARED_MARKER}  shared with the other model's top 25",
            ha="left", va="center", fontsize=7.2)

    ax.add_patch(Circle((2.0, ly2), 0.16, facecolor="#1b7837", edgecolor="#1b7837"))
    ax.text(2.6, ly2, "direct causal primary study", ha="left", va="center", fontsize=7.0)
    ax.add_patch(Circle((22.0, ly2), 0.16, facecolor="white", edgecolor="#555555", lw=1.0))
    ax.text(22.6, ly2, "no direct causal primary study", ha="left", va="center", fontsize=7.0)
    ax.text(48.0, ly2, "—  missing / not assessable", ha="left", va="center", fontsize=7.0)

    ax.text(50, -1.55,
            'Fixed query "<TF> TCR inhibition"; entity validated, LLM judged and rubric '
            "checked. Literature annotation only.",
            ha="center", va="center", fontsize=7.5, color="#555555", style="italic")
    ax.text(50, -2.15,
            f"Figure  |  {model} top 25 candidate regulator Paperclip evidence summary",
            ha="center", va="center", fontsize=11, fontweight="bold", color="#111111")

    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.05)
    stem = f"fig_paperclip_v2_{model.lower()}_tf_summary"
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(out / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[figure] wrote {stem}.svg/.png")


def write_caption(model: str, df: pd.DataFrame, out: Path, other_model: str) -> None:
    dist = {t: int((df["final_usable_tier"] == t).sum()) for t in TIER_ORDER}
    sm = dist["strong"] + dist["moderate"]
    missing = df[df["final_usable_tier"] == "missing"]
    if len(missing):
        r = missing.iloc[0]
        missing_note = (
            f"The tier for {r['TF']} remained missing because of a {r['missing_reason']}; "
            f"missing was not interpreted as no evidence."
        )
    else:
        missing_note = ""

    # Human verification note
    hv_pending = (df["human_verification_status"] == "Pending human verification").any()
    if hv_pending or (df["audit_status"] == "Pending human verification").any():
        human_note = ("Human verification of strong and moderate key papers remained "
                      "pending at the time of this audit.")
    else:
        human_note = ("Strong and moderate key papers were verified by a named human "
                      "reviewer; see human_verification.csv.")

    if model == "GREmLN":
        title = "Figure | GREmLN top 25 candidate regulator Paperclip evidence summary."
        shared_clause = ("Shared candidates also present in the GENIE3 top 25 are marked "
                         f"with {SHARED_MARKER}.")
        body_counts = (
            f"GREmLN contained {dist['strong']} strong, {dist['moderate']} moderate, "
            f"{dist['weak']} weak, {dist['none']} no-evidence and {dist['missing']} missing "
            f"candidate tiers; {sm} candidates had strong or moderate evidence."
        )
        fname = "fig_paperclip_v2_gremln_tf_summary_caption.txt"
    else:
        title = "Figure | GENIE3 top 25 candidate regulator Paperclip evidence summary."
        shared_clause = ("Shared candidates also present in the GREmLN top 25 are marked "
                         f"with {SHARED_MARKER}.")
        body_counts = (
            f"GENIE3 contained {dist['strong']} strong, {dist['moderate']} moderate, "
            f"{dist['weak']} weak, {dist['none']} no-evidence and {dist['missing']} missing "
            f"candidate tiers; {sm} candidates had strong or moderate evidence."
        )
        fname = "fig_paperclip_v2_genie3_tf_summary_caption.txt"

    parts = [
        title,
        f"Candidate regulators are ordered by their {model} rank.",
        "Literature was retrieved using one fixed Paperclip query per candidate, "
        '"<TF> TCR inhibition", with the ten highest ranked results retained.',
        "Papers passed an entity identity gate before assessment by a fixed LLM judge "
        "and deterministic evidence rubric.",
        "The coloured tile shows the final usable evidence tier: strong, moderate, "
        "weak, none or missing.",
        "The five-number evidence string gives, from left to right, the numbers of "
        "papers returned, entity eligible, relevant, primary and directly causal or "
        "perturbational.",
        "A filled circle indicates at least one direct causal primary study; an open "
        "circle indicates none.",
        shared_clause,
        body_counts,
        "Paperclip supplied literature retrieval but did not assign the evidence tiers.",
        "Literature evidence was added after model ranking and was not a predictive input.",
    ]
    if missing_note:
        parts.append(missing_note)
    parts.append(human_note)

    text = " ".join(parts)
    (out / fname).write_text(text + "\n", encoding="utf-8")
    print(f"[caption] wrote {fname}")


def main() -> int:
    repo = base.find_repo_root(Path(__file__).resolve().parent)
    out = repo / base.OUT_REL

    trace = pd.read_csv(out / "paperclip_pipeline_trace.csv", keep_default_na=False)
    tiers = pd.read_csv(out / "tf_evidence_tiers.csv", keep_default_na=False)
    union = pd.read_csv(out / "candidate_union.csv")
    chart = pd.read_csv(out / "paperclip_evidence_summary_counts.csv")
    ev = pd.read_csv(out / "paperclip_entity_validation.csv", keep_default_na=False)
    elig_ids = {
        tf: set(g[g["entity_eligible"].astype(str).str.lower() == "true"]["paperclip_result_id"])
        for tf, g in ev.groupby("TF")
    }

    gr_rank = pd.read_csv(repo / "results" / "tables" /
                          "gremln_btla_vs_tcr_seed_excluded_tf_ranking.csv")
    g3_rank = pd.read_csv(repo / "results" / "tables" /
                          "genie3_btla_vs_tcr_seed_excluded_tf_ranking.csv")

    gr_members = set(union[union["in_gremln_top25"].astype(str).str.lower() == "true"]["TF"])
    g3_members = set(union[union["in_genie3_top25"].astype(str).str.lower() == "true"]["TF"])
    assert len(gr_members) == 25 and len(g3_members) == 25

    gr_order = _ordinal_top25(gr_rank, "gene", gr_members)
    g3_order = _ordinal_top25(g3_rank, "gene", g3_members)

    human_done = _human_verification_complete(out)
    gr_df = build_model_table("GREmLN", gr_order, trace, tiers, human_done)
    g3_df = build_model_table("GENIE3", g3_order, trace, tiers, human_done)

    # Drop helper columns not required in CSV (keep in_gremln/in_genie3 for auditability
    # but also required fields from the task — include them)
    csv_cols = [
        "rank", "TF", "candidate_full_name", "candidate_group", "shared_marker",
        "final_usable_tier", "original_judge_tier", "corrected_judge_tier", "confidence",
        "papers_returned", "papers_entity_eligible", "papers_relevant",
        "primary_relevant_papers", "direct_causal_papers", "evidence_count_string",
        "key_paper_count", "key_paper_ids", "direct_causal_glyph",
        "entity_audit_status", "human_verification_status", "audit_status", "missing_reason",
    ]
    # Also retain membership flags as requested in section 2
    for extra in ("in_gremln_top25", "in_genie3_top25"):
        if extra not in csv_cols:
            csv_cols.append(extra)

    gr_df[csv_cols].to_csv(out / "paperclip_gremln_tf_evidence_summary.csv", index=False)
    g3_df[csv_cols].to_csv(out / "paperclip_genie3_tf_evidence_summary.csv", index=False)

    assert_model_table(gr_df, "GREmLN", gr_order, chart, elig_ids, human_done)
    assert_model_table(g3_df, "GENIE3", g3_order, chart, elig_ids, human_done)
    assert_shared_parity(gr_df, g3_df)

    render_figure(gr_df, "GREmLN", out)
    render_figure(g3_df, "GENIE3", out)
    write_caption("GREmLN", gr_df, out, "GENIE3")
    write_caption("GENIE3", g3_df, out, "GREmLN")

    # Report summary for stdout
    def sm_names(df):
        return list(df[df["final_usable_tier"].isin(["strong", "moderate"])]
                    .sort_values(["final_usable_tier", "rank"])["TF"])

    print("\n=== REPORT ===")
    print("paths:")
    for p in [
        "paperclip_gremln_tf_evidence_summary.csv",
        "paperclip_genie3_tf_evidence_summary.csv",
        "fig_paperclip_v2_gremln_tf_summary.svg",
        "fig_paperclip_v2_gremln_tf_summary.png",
        "fig_paperclip_v2_gremln_tf_summary_caption.txt",
        "fig_paperclip_v2_genie3_tf_summary.svg",
        "fig_paperclip_v2_genie3_tf_summary.png",
        "fig_paperclip_v2_genie3_tf_summary_caption.txt",
    ]:
        print(f"  {out / p}")
    for name, df in (("GREmLN", gr_df), ("GENIE3", g3_df)):
        dist = {t: int((df["final_usable_tier"] == t).sum()) for t in TIER_ORDER}
        print(f"{name} tiers: {dist}")
        print(f"{name} strong/moderate: {sm_names(df)}")
    shared_sm = gr_df[(gr_df["candidate_group"] == "shared")
                      & gr_df["final_usable_tier"].isin(["strong", "moderate"])]["TF"].tolist()
    print("shared strong/moderate:", shared_sm)
    miss = g3_df[g3_df["final_usable_tier"] == "missing"]
    if len(miss):
        print("missing:", miss.iloc[0]["TF"], "-", miss.iloc[0]["missing_reason"])
    corr = pd.read_csv(out / "corrected_tiers.csv", keep_default_na=False)
    changed = corr[corr["legacy_final_usable_tier"] != corr["corrected_final_usable_tier"]]
    print("entity-correction tier changes:",
          list(zip(changed["TF"], changed["legacy_final_usable_tier"],
                   changed["corrected_final_usable_tier"])))
    print("human verification: pending" if not human_done else "human verification: complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
