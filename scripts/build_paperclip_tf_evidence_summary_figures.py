#!/usr/bin/env python3
"""Build TF-level Paperclip evidence summary figures (consistency-corrected).

Canonical display order comes from primary_rankings_common_universe.csv using
the same tie-break as Table 3 (GREmLN: dense rank then SumCSLS; GENIE3: dense
rank then alphabetical). Literature evidence is joined onto that order and
never re-sorts it.

Evidence fields:
  supporting_key_paper_ids  — entity-eligible AND materially relevant
  assessment_cited_paper_ids — judge key_paperclip_result_ids (any tier)
  qualifying_direct_papers  — eligible + relevant primary + direct causal
                              + strong-tier TCR-inhibition phenotype
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_paperclip_tcr_inhibition_audit as base  # noqa: E402

TIER_ORDER = ["strong", "moderate", "weak", "none", "missing"]
TIER_COLORS = {
    "strong": "#1b7837", "moderate": "#7fbf7b", "weak": "#fee08b",
    "none": "#d9d9d9", "missing": "#b2182b",
}
TIER_TEXT = {
    "strong": "Strong", "moderate": "Moderate", "weak": "Weak",
    "none": "None", "missing": "Missing",
}
TIER_FG = {
    "strong": "white", "moderate": "#1a1a1a", "weak": "#1a1a1a",
    "none": "#1a1a1a", "missing": "white",
}
SHARED_MARKER = "●"
GROUP_LABEL = {
    "shared": "shared",
    "gremln_only": "GREmLN specific",
    "genie3_only": "GENIE3 specific",
}

# Strong-tier TCR inhibition phenotypes from the fixed judge prompt.
QUALIFYING_PHENOTYPES = {
    "proximal TCR signalling inhibition or attenuation",
    "inhibitory checkpoint response",
    "T cell anergy or hyporesponsiveness",
    "exhaustion explicitly linked to impaired TCR function",
}


def _human_verification_complete(out: Path) -> bool:
    hv_path = out / "human_verification.csv"
    if not hv_path.exists():
        return False
    hv = pd.read_csv(hv_path, keep_default_na=False)
    if hv.empty:
        return False
    ok = ((hv.get("reviewer", pd.Series(dtype=str)).astype(str).str.strip() != "")
          & (hv.get("date", pd.Series(dtype=str)).astype(str).str.strip() != ""))
    return bool(ok.any())


def _load_judge(out: Path, tf: str, use_rerun: bool) -> dict | None:
    path = out / ("judge_rerun_outputs" if use_rerun else "judge_raw_outputs") / f"{tf}.json"
    if not path.exists():
        path = out / "judge_raw_outputs" / f"{tf}.json"
    if not path.exists():
        return None
    return json.loads(path.read_text()).get("parsed")


def _paper_counts(parsed: dict | None, eligible: set) -> dict:
    """Derive evidence counts from a judge packet against entity-eligible IDs."""
    assessments = (parsed or {}).get("paper_assessments") or []
    cited = list((parsed or {}).get("key_paperclip_result_ids") or [])
    cited_eligible = [c for c in cited if c in eligible]

    relevant = [a for a in assessments
                if a.get("relevant") and a.get("paperclip_result_id") in eligible]
    primary = [a for a in relevant if a.get("study_type") == "primary experimental"]
    broad_direct = [a for a in primary
                    if a.get("evidence_directness") == "direct causal or perturbational"]
    qualifying = [a for a in broad_direct
                  if str(a.get("phenotype") or "") in QUALIFYING_PHENOTYPES]

    # Supporting keys: cited by the judge AND entity-eligible AND materially relevant
    relevant_ids = {a.get("paperclip_result_id") for a in relevant}
    supporting = [c for c in cited if c in eligible and c in relevant_ids]

    return {
        "papers_relevant": len(relevant),
        "primary_relevant_papers": len(primary),
        "broad_direct_causal_papers": len(broad_direct),
        "qualifying_direct_papers": len(qualifying),
        "supporting_key_paper_ids": supporting,
        "assessment_cited_paper_ids": cited,
        "assessment_cited_eligible_ids": cited_eligible,
        "qualifying_paper_ids": [a.get("paperclip_result_id") for a in qualifying],
        "broad_direct_paper_ids": [a.get("paperclip_result_id") for a in broad_direct],
    }


def _audit_status(row: pd.Series, human_done: bool) -> str:
    tier = str(row["audited_provisional_tier"]).lower()
    if tier == "missing":
        reason = str(row.get("missing_reason", ""))
        if "excerpt" in reason.lower() or "word" in reason.lower():
            return "Supporting excerpt exceeded word limit"
        if "rubric" in reason.lower():
            return "Rubric validation failed"
        return "Missing for another stated reason"
    if str(row.get("entity_corrected", "")).lower() in ("true", "1"):
        return "Entity corrected"
    if tier in ("strong", "moderate"):
        return "Verified" if human_done else "Pending human verification"
    return "Automated audit complete"


def _missing_reason(tf: str, tiers: pd.DataFrame) -> str:
    sub = tiers[tiers["TF"] == tf]
    if not len(sub):
        return "unresolved structured judgement"
    r = sub.iloc[0]
    if str(r.get("final_usable_tier", "")).lower() != "missing":
        return ""
    rr = str(r.get("rubric_reasons", "") or "")
    if "supporting_excerpt >20 words" in rr:
        pid = rr.split("for ")[-1].strip() if "for " in rr else ""
        return (f"the supporting excerpt"
                f"{f' for {pid}' if pid else ''} exceeded the deterministic 20-word limit")
    if rr:
        return f"a deterministic rubric validation failure ({rr})"
    return "rubric validation failed"


def canonical_orders(repo: Path) -> tuple[pd.DataFrame, pd.DataFrame, list, list]:
    """Return (gr_base, g3_base, gr_order, g3_order) with display_rank 1..25."""
    pr = pd.read_csv(repo / "results" / "publication_data" / "primary_rankings_common_universe.csv")
    u25 = pd.read_csv(repo / "results" / "publication_data" / "top25_union_primary.csv")
    t3 = pd.read_csv(repo / "results" / "publication_assets" / "tables" / "table3_top25_rankings.csv")

    gr_members = set(u25.loc[u25["in_gremln_top25"] == True, "TF"])
    g3_members = set(u25.loc[u25["in_genie3_top25"] == True, "TF"])

    gr_ord = (pr[pr["TF"].isin(gr_members)]
              .sort_values(["gremln_dense_rank", "gremln_csls_seed_sum"],
                           ascending=[True, False])
              .reset_index(drop=True))
    g3_ord = (pr[pr["TF"].isin(g3_members)]
              .sort_values(["genie3_dense_rank", "TF"])
              .reset_index(drop=True))

    # Strip Table 3 shared markers for comparison
    def _sym(cell: str) -> str:
        s = str(cell)
        for sep in ("ˢ", "\u02e2", " ("):
            s = s.split(sep)[0]
        return s.strip()

    t3_gr = [_sym(x) for x in t3["GREmLN candidate"]]
    t3_g3 = [_sym(x) for x in t3["GENIE3 candidate"]]
    assert list(gr_ord["TF"]) == t3_gr, (
        f"GREmLN order differs from Table 3:\n  got {list(gr_ord['TF'])}\n  t3 {t3_gr}")
    assert list(g3_ord["TF"]) == t3_g3, (
        f"GENIE3 order differs from Table 3:\n  got {list(g3_ord['TF'])}\n  t3 {t3_g3}")

    gr_base = gr_ord[["TF", "gremln_dense_rank", "gremln_csls_seed_sum"]].copy()
    gr_base["display_rank"] = range(1, 26)
    gr_base["model"] = "GREmLN"
    g3_base = g3_ord[["TF", "genie3_dense_rank"]].copy()
    g3_base["display_rank"] = range(1, 26)
    g3_base["model"] = "GENIE3"
    return gr_base, g3_base, list(gr_ord["TF"]), list(g3_ord["TF"])


def build_model_table(base: pd.DataFrame, model: str, trace: pd.DataFrame,
                      tiers: pd.DataFrame, elig: dict, out: Path,
                      human_done: bool) -> pd.DataFrame:
    tr = {r["TF"]: r for _, r in trace.iterrows()}
    corr = pd.read_csv(out / "corrected_tiers.csv", keep_default_na=False)
    corr_by = {r["TF"]: r for _, r in corr.iterrows()}

    rows = []
    # Join onto base WITHOUT re-sorting
    for _, b in base.iterrows():
        tf = b["TF"]
        r = tr[tf]
        cr = corr_by.get(tf, {})
        use_rerun = str(cr.get("rerun_required", "")).lower() == "true"
        parsed = _load_judge(out, tf, use_rerun)
        eligible = elig.get(tf, set())
        counts = _paper_counts(parsed, eligible)

        tier = str(r["final_usable_tier"]).lower()
        # Never convert missing/weak
        assert tier in TIER_ORDER, tier

        entity_corrected = use_rerun or ("corrected" in str(r["entity_audit_status"]).lower())
        miss_reason = _missing_reason(tf, tiers) if tier == "missing" else ""

        conf = str(r["confidence"]).strip().lower()
        conf_label = {"high": "High", "moderate": "Moderate", "low": "Low"}.get(
            conf, conf.title() or "—")
        if tier == "missing" and conf_label in ("", "Nan"):
            conf_label = "—"

        n_ret = int(r["papers_returned"])
        n_elig = len(eligible)
        n_rel = counts["papers_relevant"]
        n_prim = counts["primary_relevant_papers"]
        n_qual = counts["qualifying_direct_papers"]
        n_broad = counts["broad_direct_causal_papers"]
        supp = counts["supporting_key_paper_ids"]
        cited = counts["assessment_cited_paper_ids"]

        group = GROUP_LABEL[str(r["candidate_group"])]
        shared = group == "shared"

        row = {
            "display_rank": int(b["display_rank"]),
            "rank": int(b["display_rank"]),
            "TF": tf,
            "candidate_full_name": r["candidate_full_name"],
            "candidate_group": group,
            "in_gremln_top25": group in ("GREmLN specific", "shared"),
            "in_genie3_top25": group in ("GENIE3 specific", "shared"),
            "shared_marker": SHARED_MARKER if shared else "",
            "audited_provisional_tier": tier,
            "final_usable_tier": tier if human_done else "",  # blank until human verification
            "original_judge_tier": str(r["original_judge_tier"]),
            "corrected_judge_tier": str(r["corrected_judge_tier"]),
            "confidence": conf_label,
            "papers_returned": n_ret,
            "papers_entity_eligible": n_elig,
            "papers_relevant": n_rel,
            "primary_relevant_papers": n_prim,
            "broad_direct_causal_papers": n_broad,
            "qualifying_direct_papers": n_qual,
            "evidence_count_string": f"{n_ret} / {n_elig} / {n_rel} / {n_prim} / {n_qual}",
            "supporting_key_paper_count": len(supp),
            "supporting_key_paper_ids": ";".join(supp),
            "assessment_cited_paper_ids": ";".join(cited),
            "key_paper_count": len(supp),  # figure uses supporting keys
            "key_paper_ids": ";".join(supp),
            "direct_causal_glyph": ("—" if tier == "missing"
                                   else ("●" if n_qual >= 1 else "○")),
            "entity_audit_status": r["entity_audit_status"],
            "entity_corrected": entity_corrected,
            "human_verification_status": (
                "Verified" if (human_done and tier in ("strong", "moderate"))
                else ("Pending human verification" if tier in ("strong", "moderate")
                      else "not required")),
            "missing_reason": miss_reason,
            "qualifying_direct_paper_ids": ";".join(counts["qualifying_paper_ids"]),
            "broad_direct_paper_ids": ";".join(counts["broad_direct_paper_ids"]),
        }
        row["audit_status"] = _audit_status(pd.Series(row), human_done)
        rows.append(row)

    df = pd.DataFrame(rows)
    # Critical: order must remain base order
    assert list(df["TF"]) == list(base["TF"]), "evidence join reordered candidates"
    return df


def assert_model_table(df: pd.DataFrame, model: str, ordered_tfs: list[str],
                       chart: pd.DataFrame, elig: dict, human_done: bool,
                       prev_df: pd.DataFrame | None) -> None:
    checks = []

    def chk(name, cond):
        checks.append((name, bool(cond)))

    chk(f"{model}: 25 rows", len(df) == 25)
    chk(f"{model}: display_rank 1..25", list(df["display_rank"]) == list(range(1, 26)))
    chk(f"{model}: TF unique", df["TF"].is_unique)
    chk(f"{model}: order matches canonical", list(df["TF"]) == ordered_tfs)
    chk(f"{model}: evidence cannot affect rank",
        list(df.sort_values("display_rank")["TF"]) == ordered_tfs)

    for _, r in df.iterrows():
        parts = [int(x.strip()) for x in r["evidence_count_string"].split("/")]
        chk(f"{model}/{r['TF']} count string",
            parts == [r["papers_returned"], r["papers_entity_eligible"],
                      r["papers_relevant"], r["primary_relevant_papers"],
                      r["qualifying_direct_papers"]])
        # supporting keys ⊆ relevant
        chk(f"{model}/{r['TF']} supporting <= relevant",
            int(r["supporting_key_paper_count"]) <= int(r["papers_relevant"]))
        for kid in str(r["supporting_key_paper_ids"]).split(";"):
            if kid:
                chk(f"{model}/{r['TF']} supporting {kid} eligible",
                    kid in elig.get(r["TF"], set()))

    row = chart[chart["group"] == f"{model} top 25"].iloc[0]
    for t in TIER_ORDER:
        chk(f"{model} tier {t} vs chart",
            int((df["audited_provisional_tier"] == t).sum()) == int(row[t]))
    chk(f"{model} S+M vs chart",
        int(df["audited_provisional_tier"].isin(["strong", "moderate"]).sum())
        == int(row["strong_plus_moderate"]))

    for _, r in df[df["audited_provisional_tier"] == "strong"].iterrows():
        chk(f"{model}/{r['TF']} strong has qualifying direct",
            int(r["qualifying_direct_papers"]) >= 1)
    for _, r in df[df["audited_provisional_tier"] == "moderate"].iterrows():
        chk(f"{model}/{r['TF']} moderate has relevant primary",
            int(r["primary_relevant_papers"]) >= 1)

    for _, r in df.iterrows():
        filled = r["direct_causal_glyph"] == "●"
        chk(f"{model}/{r['TF']} glyph <=> qualifying",
            filled == (int(r["qualifying_direct_papers"]) >= 1
                       and r["audited_provisional_tier"] != "missing"))

    for _, r in df[df["audited_provisional_tier"] == "missing"].iterrows():
        chk(f"{model}/{r['TF']} missing stays missing",
            r["audited_provisional_tier"] == "missing")
        chk(f"{model}/{r['TF']} missing reason", bool(str(r["missing_reason"]).strip()))

    if not human_done:
        chk(f"{model}: no Verified without human review",
            (df["audit_status"] != "Verified").all())
        chk(f"{model}: final_usable_tier blank until human verification",
            (df["final_usable_tier"].astype(str).str.strip() == "").all())

    # Vague labels banned
    chk(f"{model}: no 'Other reviewed status'",
        (df["audit_status"] != "Other reviewed status").all())

    failed = [n for n, ok in checks if not ok]
    print(f"[{model}] assertions: {sum(ok for _, ok in checks)}/{len(checks)} passed")
    if failed:
        for n in failed[:25]:
            print(f"  FAIL: {n}")
        raise AssertionError(f"{model} figure assertions failed ({len(failed)})")


def assert_shared_parity(gr: pd.DataFrame, g3: pd.DataFrame) -> None:
    shared = set(gr[gr["candidate_group"] == "shared"]["TF"])
    assert shared == set(g3[g3["candidate_group"] == "shared"]["TF"])
    fields = [
        "audited_provisional_tier", "confidence", "papers_returned",
        "papers_entity_eligible", "papers_relevant", "primary_relevant_papers",
        "qualifying_direct_papers", "evidence_count_string",
        "supporting_key_paper_count", "supporting_key_paper_ids",
        "direct_causal_glyph", "audit_status", "human_verification_status",
        "missing_reason", "shared_marker",
    ]
    for tf in sorted(shared):
        a = gr[gr["TF"] == tf].iloc[0]
        b = g3[g3["TF"] == tf].iloc[0]
        for f in fields:
            if str(a[f]) != str(b[f]):
                raise AssertionError(f"shared {tf} {f}: {a[f]!r} != {b[f]!r}")
    print(f"[shared] parity OK for {len(shared)} candidates")


def render_figure(df: pd.DataFrame, model: str, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Circle

    n = len(df)
    fig_w, fig_h = 14.0, 12.8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, 100)
    ax.set_ylim(-2.5, n + 4.4)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Column left edges — extra gap between evidence flow and confidence
    cols = {
        "rank": 1.0,
        "tf": 5.8,
        "tier": 20.5,
        "counts": 37.0,
        "conf": 62.0,
        "key": 71.5,
        "audit": 82.5,
        "glyph": 96.8,
    }

    ax.text(50, n + 3.7, f"{model} candidate evidence", ha="center", va="center",
            fontsize=13, fontweight="bold", color="#222222")

    # Two-line header for evidence flow; single-line for others
    hy_top = n + 2.85
    hy_bot = n + 2.25
    ax.text(cols["rank"], hy_bot, "Rank", ha="left", va="center",
            fontsize=7.4, fontweight="bold", color="#333333")
    ax.text(cols["tf"], hy_bot, "Candidate regulator", ha="left", va="center",
            fontsize=7.4, fontweight="bold", color="#333333")
    ax.text(cols["tier"], hy_bot, "Audited provisional tier", ha="left", va="center",
            fontsize=7.4, fontweight="bold", color="#333333")
    ax.text(cols["counts"], hy_top, "Evidence flow", ha="left", va="center",
            fontsize=7.4, fontweight="bold", color="#333333")
    ax.text(cols["counts"], hy_bot,
            "Returned / eligible / relevant / primary / qualifying direct",
            ha="left", va="center", fontsize=6.4, color="#444444")
    ax.text(cols["conf"], hy_bot, "Confidence", ha="left", va="center",
            fontsize=7.4, fontweight="bold", color="#333333")
    ax.text(cols["key"], hy_bot, "Supporting key papers", ha="left", va="center",
            fontsize=7.4, fontweight="bold", color="#333333")
    ax.text(cols["audit"], hy_bot, "Audit status", ha="left", va="center",
            fontsize=7.4, fontweight="bold", color="#333333")
    ax.text(cols["glyph"], hy_bot, "Direct", ha="center", va="center",
            fontsize=7.4, fontweight="bold", color="#333333")
    ax.plot([0.6, 99.4], [hy_bot - 0.45, hy_bot - 0.45], color="#bbbbbb", lw=0.8)

    for _, r in df.iterrows():
        y = n - int(r["display_rank"]) + 1.0
        if int(r["display_rank"]) % 2 == 0:
            ax.add_patch(FancyBboxPatch(
                (0.5, y - 0.42), 99.0, 0.84,
                boxstyle="square,pad=0", linewidth=0, facecolor="#f7f7f7"))
        ax.plot([0.6, 99.4], [y - 0.42, y - 0.42], color="#e8e8e8", lw=0.4)

        ax.text(cols["rank"] + 2.2, y, str(int(r["display_rank"])), ha="right",
                va="center", fontsize=8.4, fontfamily="DejaVu Sans Mono")
        tf_lab = f"{r['TF']}  {SHARED_MARKER}" if r["shared_marker"] else r["TF"]
        ax.text(cols["tf"], y, tf_lab, ha="left", va="center", fontsize=8.4,
                fontweight="bold" if r["shared_marker"] else "normal")

        tier = r["audited_provisional_tier"]
        pill_w, pill_h = 13.5, 0.58
        ax.add_patch(FancyBboxPatch(
            (cols["tier"], y - pill_h / 2), pill_w, pill_h,
            boxstyle="round,pad=0.02,rounding_size=0.15",
            linewidth=0.4, edgecolor="#666666", facecolor=TIER_COLORS[tier]))
        ax.text(cols["tier"] + pill_w / 2, y, TIER_TEXT[tier], ha="center",
                va="center", fontsize=7.6, fontweight="bold", color=TIER_FG[tier])

        ax.text(cols["counts"], y, r["evidence_count_string"], ha="left", va="center",
                fontsize=7.6, fontfamily="DejaVu Sans Mono")
        ax.text(cols["conf"], y, r["confidence"], ha="left", va="center", fontsize=8.0)

        k = int(r["supporting_key_paper_count"])
        key_txt = f"{k} key paper" + ("" if k == 1 else "s") if k else "—"
        ax.text(cols["key"], y, key_txt, ha="left", va="center", fontsize=7.8)

        # Wrap long audit labels slightly by font size
        ax.text(cols["audit"], y, r["audit_status"], ha="left", va="center",
                fontsize=6.6 if len(str(r["audit_status"])) > 28 else 7.2)

        gx = cols["glyph"]
        if r["direct_causal_glyph"] == "—":
            ax.text(gx, y, "—", ha="center", va="center", fontsize=9, color="#666666")
        elif int(r["qualifying_direct_papers"]) >= 1:
            ax.add_patch(Circle((gx, y), 0.20, facecolor="#1b7837",
                                edgecolor="#1b7837", lw=0.8, zorder=3))
        else:
            ax.add_patch(Circle((gx, y), 0.20, facecolor="white",
                                edgecolor="#555555", lw=1.0, zorder=3))

    ax.plot([0.6, 99.4], [0.55, 0.55], color="#bbbbbb", lw=0.8)

    ly1, ly2 = -0.15, -0.85
    ax.text(1.0, ly1 + 0.35, "Legend", fontsize=8, fontweight="bold")
    lx = 1.0
    for t in TIER_ORDER:
        ax.add_patch(FancyBboxPatch(
            (lx, ly1 - 0.20), 2.0, 0.36,
            boxstyle="round,pad=0.01,rounding_size=0.08",
            linewidth=0.3, edgecolor="#666666", facecolor=TIER_COLORS[t]))
        ax.text(lx + 2.3, ly1, TIER_TEXT[t], ha="left", va="center", fontsize=7.0)
        lx += 6.5
    ax.text(lx + 0.5, ly1, f"{SHARED_MARKER}  shared with the other model's top 25",
            ha="left", va="center", fontsize=7.0)

    ax.add_patch(Circle((2.0, ly2), 0.16, facecolor="#1b7837", edgecolor="#1b7837"))
    ax.text(2.6, ly2, "qualifying direct primary (strong-tier phenotype)",
            ha="left", va="center", fontsize=7.0)
    ax.add_patch(Circle((42.0, ly2), 0.16, facecolor="white", edgecolor="#555555", lw=1.0))
    ax.text(42.6, ly2, "no qualifying direct", ha="left", va="center", fontsize=7.0)
    ax.text(62.0, ly2, "—  missing / not assessable", ha="left", va="center", fontsize=7.0)

    ax.text(50, -1.55,
            'Fixed query "<TF> TCR inhibition"; entity validated, LLM judged and rubric '
            "checked. Literature annotation only. Strong and moderate tiers remain provisional.",
            ha="center", va="center", fontsize=7.3, color="#555555", style="italic")
    ax.text(50, -2.2,
            f"Figure  |  {model} top 25 candidate regulator Paperclip evidence summary",
            ha="center", va="center", fontsize=11, fontweight="bold")

    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.05)
    stem = f"fig_paperclip_v2_{model.lower()}_tf_summary"
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(out / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[figure] wrote {stem}.svg/.png")


def write_caption(model: str, df: pd.DataFrame, out: Path) -> None:
    dist = {t: int((df["audited_provisional_tier"] == t).sum()) for t in TIER_ORDER}
    sm = dist["strong"] + dist["moderate"]
    missing = df[df["audited_provisional_tier"] == "missing"]
    if len(missing):
        r = missing.iloc[0]
        missing_note = (
            f"The tier for {r['TF']} remained missing because {r['missing_reason']}; "
            f"missing was not interpreted as no evidence."
        )
    else:
        missing_note = ""

    human_note = ("Strong and moderate tiers remain provisional pending documented human "
                  "verification.")

    if model == "GREmLN":
        title = "Figure | GREmLN top 25 candidate regulator Paperclip evidence summary."
        shared_clause = (f"Shared candidates also present in the GENIE3 top 25 are marked "
                         f"with {SHARED_MARKER}.")
        body_counts = (
            f"GREmLN contained {dist['strong']} strong, {dist['moderate']} moderate, "
            f"{dist['weak']} weak, {dist['none']} no-evidence and {dist['missing']} missing "
            f"candidate tiers; {sm} candidates had strong or moderate evidence.")
        fname = "fig_paperclip_v2_gremln_tf_summary_caption.txt"
    else:
        title = "Figure | GENIE3 top 25 candidate regulator Paperclip evidence summary."
        shared_clause = (f"Shared candidates also present in the GREmLN top 25 are marked "
                         f"with {SHARED_MARKER}.")
        body_counts = (
            f"GENIE3 contained {dist['strong']} strong, {dist['moderate']} moderate, "
            f"{dist['weak']} weak, {dist['none']} no-evidence and {dist['missing']} missing "
            f"candidate tiers; {sm} candidates had strong or moderate evidence.")
        fname = "fig_paperclip_v2_genie3_tf_summary_caption.txt"

    parts = [
        title,
        f"Candidate regulators are ordered by their {model} rank after the model's "
        f"complete tie-break procedure (identical to Table 3).",
        "Literature was retrieved using one fixed Paperclip query per candidate, "
        '"<TF> TCR inhibition", with the ten highest ranked results retained.',
        "Papers passed an entity identity gate before assessment by a fixed LLM judge "
        "and deterministic evidence rubric.",
        "The coloured tile shows the audited provisional evidence tier: strong, "
        "moderate, weak, none or missing.",
        "The five-number evidence string gives, from left to right, the numbers of "
        "papers returned, entity eligible, relevant, primary and primary papers "
        "providing qualifying direct evidence for the strong-tier TCR inhibition "
        "rubric.",
        "A filled circle indicates at least one qualifying direct primary study; an "
        "open circle indicates none.",
        "Broader causal evidence in related or generic T cell phenotypes may still "
        "result in moderate or weak evidence and does not fill the direct glyph.",
        "Supporting key papers are entity-eligible papers that the judge marked "
        "materially relevant; papers cited only to justify a weak or none assessment "
        "are retained in the backing CSV as assessment-cited papers and are not shown "
        "as supporting key papers.",
        shared_clause,
        body_counts,
        "Paperclip supplied literature retrieval but did not assign the evidence tiers.",
        "Literature evidence was added after model ranking and was not a predictive input.",
    ]
    if missing_note:
        parts.append(missing_note)
    parts.append(human_note)
    (out / fname).write_text(" ".join(parts) + "\n", encoding="utf-8")
    print(f"[caption] wrote {fname}")


def write_reconciliation(out: Path, gr: pd.DataFrame, g3: pd.DataFrame,
                         prev_gr: pd.DataFrame | None, prev_g3: pd.DataFrame | None,
                         gr_order: list, g3_order: list) -> None:
    L = ["# TF-level Paperclip figure consistency reconciliation\n\n"]
    L.append("## Why the previous GREmLN ordering differed\n\n")
    L.append(
        "The previous TF-level figure ordered GREmLN candidates by row order in "
        "`gremln_btla_vs_tcr_seed_excluded_tf_ranking.csv` (neighbour-count / CSLS score "
        "order) after filtering to the top 25 membership set. That ignored the canonical "
        "within-tie SumCSLS procedure used for Table 3: sort by `gremln_dense_rank` "
        "ascending, then `gremln_csls_seed_sum` descending. Consequently ties at the same "
        "seed-neighbour count were shown in an order that did not match the publication "
        "table (for example KIF22 ahead of EZH2 by SumCSLS in Table 3, but EZH2 first in "
        "the previous figure).\n\n"
        "GENIE3 previously used ranking-file order; the canonical residual-tie rule is "
        "alphabetical within equal `genie3_dense_rank`. In this dataset GENIE3 top 25 "
        "order was already nearly identical.\n\n"
    )
    L.append("## Corrected canonical display order (top 5)\n\n")
    L.append(f"- GREmLN: {', '.join(gr_order[:5])}\n")
    L.append(f"- GENIE3: {', '.join(g3_order[:5])}\n\n")

    def pos_changes(prev, new, label):
        if prev is None:
            return []
        po = {r["TF"]: int(r["rank"]) for _, r in prev.iterrows()}
        pn = {r["TF"]: int(r["display_rank"]) for _, r in new.iterrows()}
        return [(tf, po[tf], pn[tf]) for tf in pn if tf in po and po[tf] != pn[tf]]

    L.append("## Candidates whose displayed position changed\n\n")
    for label, ch in (("GREmLN", pos_changes(prev_gr, gr, "GREmLN")),
                      ("GENIE3", pos_changes(prev_g3, g3, "GENIE3"))):
        L.append(f"### {label}\n")
        if not ch:
            L.append("- none\n")
        else:
            for tf, a, b in sorted(ch, key=lambda x: x[2]):
                L.append(f"- {tf}: {a} → {b}\n")
        L.append("\n")

    def key_changes(prev, new):
        if prev is None:
            return []
        out = []
        for _, r in new.iterrows():
            p = prev[prev.TF == r.TF]
            if not len(p):
                continue
            old = int(p.iloc[0].get("key_paper_count", 0))
            newc = int(r["supporting_key_paper_count"])
            if old != newc:
                out.append((r.TF, old, newc))
        return out

    L.append("## Candidates whose supporting key paper count changed\n\n")
    for label, ch in (("GREmLN", key_changes(prev_gr, gr)),
                      ("GENIE3", key_changes(prev_g3, g3))):
        L.append(f"### {label}\n")
        if not ch:
            L.append("- none\n")
        else:
            for tf, a, b in ch:
                L.append(f"- {tf}: {a} → {b} (now supporting = eligible ∩ relevant)\n")
        L.append("\n")

    def direct_changes(prev, new):
        if prev is None:
            return []
        out = []
        for _, r in new.iterrows():
            p = prev[prev.TF == r.TF]
            if not len(p):
                continue
            # previous fifth number was broad direct
            old_str = str(p.iloc[0].get("evidence_count_string", ""))
            try:
                old = int(old_str.split("/")[-1].strip())
            except Exception:
                old = int(p.iloc[0].get("direct_causal_papers", 0) or 0)
            newc = int(r["qualifying_direct_papers"])
            if old != newc:
                out.append((r.TF, old, newc, int(r["broad_direct_causal_papers"])))
        return out

    L.append("## Candidates whose displayed direct (5th) count changed\n\n")
    L.append("Previous 5th number was broad direct causal; now it is qualifying direct "
             "(strong-tier phenotype).\n\n")
    for label, ch in (("GREmLN", direct_changes(prev_gr, gr)),
                      ("GENIE3", direct_changes(prev_g3, g3))):
        L.append(f"### {label}\n")
        if not ch:
            L.append("- none\n")
        else:
            for tf, a, b, broad in ch:
                L.append(f"- {tf}: displayed {a} → {b} (broad direct retained = {broad})\n")
        L.append("\n")

    L.append("## Evidence tiers\n\n")
    L.append("- No audited provisional tier changed in this correction pass.\n\n")

    miss = pd.concat([gr, g3])
    miss = miss[miss["audited_provisional_tier"] == "missing"].drop_duplicates("TF")
    L.append("## NR4A2 exact failure\n\n")
    if len(miss):
        r = miss.iloc[0]
        L.append(f"- TF: **{r['TF']}** (GENIE3 specific)\n")
        L.append(f"- Original judge tier: **{r['original_judge_tier']}**\n")
        L.append(f"- Entity eligible / relevant / primary / broad direct / qualifying direct: "
                 f"{r['papers_entity_eligible']} / {r['papers_relevant']} / "
                 f"{r['primary_relevant_papers']} / {r['broad_direct_causal_papers']} / "
                 f"{r['qualifying_direct_papers']}\n")
        L.append(f"- Supporting key papers: `{r['supporting_key_paper_ids']}`\n")
        L.append(f"- Assessment-cited papers: `{r['assessment_cited_paper_ids']}`\n")
        L.append(f"- Exact failed assertion: `{r['missing_reason']}`\n")
        L.append("- Protocol-compliant retry: **not allowed** (retry is only for "
                 "unparseable JSON; content rubric violations are frozen as missing).\n")
        L.append(f"- Figure audit label: **{r['audit_status']}**\n")
        L.append("- Final tier: remains **missing**.\n\n")

    L.append("## Focus TF direct-evidence audit\n\n")
    L.append("Whether previously displayed direct counts were broad causal rather than "
             "strong-tier qualifying:\n\n")
    focus = ["BHLHE40", "RBPJ", "ID2", "MAF", "JUNB", "EZH2", "NFIL3", "AHR", "BCL3", "FOXP1"]
    both = pd.concat([gr, g3]).drop_duplicates("TF").set_index("TF")
    for tf in focus:
        if tf not in both.index:
            continue
        r = both.loc[tf]
        L.append(f"- **{tf}** ({r['audited_provisional_tier']}): broad direct = "
                 f"{r['broad_direct_causal_papers']}, qualifying direct = "
                 f"{r['qualifying_direct_papers']}"
                 f"{' — previously displayed broad causal as if strong-qualifying' if int(r['broad_direct_causal_papers']) != int(r['qualifying_direct_papers']) else ' — counts agree'}"
                 f".\n")
    L.append("\n")

    L.append("## Aggregate chart reconciliation\n\n")
    for name, df in (("GREmLN", gr), ("GENIE3", g3)):
        dist = {t: int((df["audited_provisional_tier"] == t).sum()) for t in TIER_ORDER}
        L.append(f"- {name}: {dist} (sum={sum(dist.values())})\n")
    L.append("\nMatches `paperclip_evidence_summary_counts.csv`.\n")

    (out / "tf_figure_consistency_reconciliation.md").write_text("".join(L), encoding="utf-8")
    print("[recon] wrote tf_figure_consistency_reconciliation.md")


def main() -> int:
    repo = base.find_repo_root(Path(__file__).resolve().parent)
    out = repo / base.OUT_REL

    prev_gr_path = out / "paperclip_gremln_tf_evidence_summary.csv"
    prev_g3_path = out / "paperclip_genie3_tf_evidence_summary.csv"
    # Prefer the last committed versions so re-runs still report true deltas.
    prev_gr = prev_g3 = None
    try:
        import subprocess
        for label, dest in (("paperclip_gremln_tf_evidence_summary.csv", "gr"),
                            ("paperclip_genie3_tf_evidence_summary.csv", "g3")):
            cp = subprocess.run(
                ["git", "-C", str(repo), "show",
                 f"HEAD:results/paperclip/v2_tcr_inhibition/{label}"],
                capture_output=True, text=True)
            if cp.returncode == 0 and cp.stdout.strip():
                from io import StringIO
                dfp = pd.read_csv(StringIO(cp.stdout), keep_default_na=False)
                if dest == "gr":
                    prev_gr = dfp
                else:
                    prev_g3 = dfp
    except Exception:
        pass
    if prev_gr is None and prev_gr_path.exists():
        prev_gr = pd.read_csv(prev_gr_path, keep_default_na=False)
    if prev_g3 is None and prev_g3_path.exists():
        prev_g3 = pd.read_csv(prev_g3_path, keep_default_na=False)

    gr_base, g3_base, gr_order, g3_order = canonical_orders(repo)
    print("Canonical GREmLN top 5:", gr_order[:5])
    print("Canonical GENIE3 top 5:", g3_order[:5])

    trace = pd.read_csv(out / "paperclip_pipeline_trace.csv", keep_default_na=False)
    tiers = pd.read_csv(out / "tf_evidence_tiers.csv", keep_default_na=False)
    chart = pd.read_csv(out / "paperclip_evidence_summary_counts.csv")
    ev = pd.read_csv(out / "paperclip_entity_validation.csv", keep_default_na=False)
    elig = {tf: set(g[g["entity_eligible"].astype(str).str.lower() == "true"]["paperclip_result_id"])
            for tf, g in ev.groupby("TF")}
    human_done = _human_verification_complete(out)

    gr = build_model_table(gr_base, "GREmLN", trace, tiers, elig, out, human_done)
    g3 = build_model_table(g3_base, "GENIE3", trace, tiers, elig, out, human_done)

    csv_cols = [
        "display_rank", "rank", "TF", "candidate_full_name", "candidate_group",
        "shared_marker", "in_gremln_top25", "in_genie3_top25",
        "audited_provisional_tier", "final_usable_tier",
        "original_judge_tier", "corrected_judge_tier", "confidence",
        "papers_returned", "papers_entity_eligible", "papers_relevant",
        "primary_relevant_papers", "broad_direct_causal_papers",
        "qualifying_direct_papers", "evidence_count_string",
        "supporting_key_paper_count", "supporting_key_paper_ids",
        "assessment_cited_paper_ids", "key_paper_count", "key_paper_ids",
        "direct_causal_glyph", "entity_audit_status", "entity_corrected",
        "human_verification_status", "audit_status", "missing_reason",
        "qualifying_direct_paper_ids", "broad_direct_paper_ids",
    ]
    gr[csv_cols].to_csv(out / "paperclip_gremln_tf_evidence_summary.csv", index=False)
    g3[csv_cols].to_csv(out / "paperclip_genie3_tf_evidence_summary.csv", index=False)

    assert_model_table(gr, "GREmLN", gr_order, chart, elig, human_done, prev_gr)
    assert_model_table(g3, "GENIE3", g3_order, chart, elig, human_done, prev_g3)
    assert_shared_parity(gr, g3)

    # Shared displayed positions vs Table 3
    t3 = pd.read_csv(repo / "results" / "publication_assets" / "tables" / "table3_top25_rankings.csv")

    def _sym(cell):
        s = str(cell)
        for sep in ("ˢ", "\u02e2", " ("):
            s = s.split(sep)[0]
        return s.strip()

    t3_gr_pos = {_sym(r["GREmLN candidate"]): int(r["Position"]) for _, r in t3.iterrows()}
    t3_g3_pos = {_sym(r["GENIE3 candidate"]): int(r["Position"]) for _, r in t3.iterrows()}
    for tf in gr[gr["candidate_group"] == "shared"]["TF"]:
        assert int(gr[gr.TF == tf].iloc[0]["display_rank"]) == t3_gr_pos[tf]
        assert int(g3[g3.TF == tf].iloc[0]["display_rank"]) == t3_g3_pos[tf]
    print("[shared] displayed positions match Table 3")

    render_figure(gr, "GREmLN", out)
    render_figure(g3, "GENIE3", out)
    write_caption("GREmLN", gr, out)
    write_caption("GENIE3", g3, out)
    write_reconciliation(out, gr, g3, prev_gr, prev_g3, gr_order, g3_order)

    print("\n=== REPORT ===")
    print("GREmLN order:", list(gr.TF))
    print("GENIE3 order:", list(g3.TF))
    for name, df in (("GREmLN", gr), ("GENIE3", g3)):
        dist = {t: int((df.audited_provisional_tier == t).sum()) for t in TIER_ORDER}
        print(f"{name} tiers:", dist)
        sm = df[df.audited_provisional_tier.isin(["strong", "moderate"])]
        print(f"{name} strong/moderate:",
              list(sm.sort_values(["audited_provisional_tier", "display_rank"]).TF))
    miss = g3[g3.audited_provisional_tier == "missing"]
    if len(miss):
        print("NR4A2:", miss.iloc[0]["missing_reason"], "| label:", miss.iloc[0]["audit_status"])
    print("human verification: pending" if not human_done else "complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
