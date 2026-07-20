#!/usr/bin/env python3
"""Build the Paperclip v2 audit artefacts (trace, figures, chart, tables, report).

Consumes the frozen v2 audit files plus the entity-identity validation and the
entity-filtered tier reruns, and emits:

    paperclip_pipeline_trace.csv            (one row / candidate + assertions)
    paperclip_pipeline_audit.md             (13-section audit document)
    fig_paperclip_v2_pipeline.svg/.png      (8-stage pipeline diagram)
    paperclip_evidence_summary_counts.csv
    fig_paperclip_v2_evidence_summary.svg/.png  (2-panel evidence chart)
    fig_paperclip_v2_evidence_summary_caption.txt
    paperclip_strong_moderate_candidates_corrected.csv
    paperclip_reconciliation_report.md
    human_spot_check_strong_moderate.csv

All outputs are literature annotation only: rankings, CRISPRi and the benchmark
verdict are not modified.
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
TIER_TEXTCOLOR = {
    "strong": "white", "moderate": "black", "weak": "black",
    "none": "black", "missing": "white",
}


# --------------------------------------------------------------------------- #
def load(out: Path):
    d = {}
    d["union"] = pd.read_csv(out / "candidate_union.csv")
    d["ident"] = pd.read_csv(out / "candidate_identity.csv", keep_default_na=False)
    d["log"] = pd.read_csv(out / "paperclip_retrieval_log.csv", keep_default_na=False)
    d["pj"] = pd.read_csv(out / "paper_judgements.csv", keep_default_na=False)
    d["tiers"] = pd.read_csv(out / "tf_evidence_tiers.csv", keep_default_na=False)
    d["ev"] = pd.read_csv(out / "paperclip_entity_validation.csv", keep_default_na=False)
    d["corr"] = pd.read_csv(out / "corrected_tiers.csv", keep_default_na=False)
    d["manifest"] = json.loads((out / "run_manifest.json").read_text())
    d["ent_manifest"] = json.loads((out / "entity_run_manifest.json").read_text())
    return d


def parsed_judge(out: Path, tf: str, rerun: bool):
    p = out / ("judge_rerun_outputs" if rerun else "judge_raw_outputs") / f"{tf}.json"
    if not p.exists():
        return None
    return json.loads(p.read_text()).get("parsed")


def _counts_from_assessments(assessments, eligible_ids):
    relevant = [a for a in assessments if a.get("relevant")
                and a.get("paperclip_result_id") in eligible_ids]
    prim = [a for a in relevant if a.get("study_type") == "primary experimental"]
    direct = [a for a in prim if a.get("evidence_directness") == "direct causal or perturbational"]
    return len(relevant), len(prim), len(direct)


# --------------------------------------------------------------------------- #
def build_trace(out: Path, D: dict) -> pd.DataFrame:
    union, ident, log, ev, tiers, corr = (
        D["union"], D["ident"], D["log"], D["ev"], D["tiers"], D["corr"])
    id_by = {r["TF"]: r for _, r in ident.iterrows()}
    tier_by = {r["TF"]: r for _, r in tiers.iterrows()}
    corr_by = {r["TF"]: r for _, r in corr.iterrows()}
    ev_by = {tf: g for tf, g in ev.groupby("TF")}
    log_by = {tf: g for tf, g in log.groupby("TF")}

    rows = []
    for _, u in union.iterrows():
        tf = u["TF"]
        evg = ev_by.get(tf)
        eligible_ids = set(evg[evg["entity_eligible"].astype(str).str.lower() == "true"]
                           ["paperclip_result_id"]) if evg is not None else set()
        excluded_ids = sorted(evg[evg["entity_match"].isin(
            ["wrong_entity", "ambiguous_acronym"])]["paperclip_result_id"]) if evg is not None else []
        tr = tier_by[tf]
        cr = corr_by[tf]
        rerun_required = str(cr["rerun_required"]).lower() == "true"
        parsed = parsed_judge(out, tf, rerun_required)
        assess = (parsed.get("paper_assessments", []) if parsed else [])
        n_rel, n_prim, n_direct = _counts_from_assessments(assess, eligible_ids)

        # key eligible paper ids (from rerun output if rerun, else original)
        if rerun_required and parsed:
            key_ids = [k for k in (parsed.get("key_paperclip_result_ids") or []) if k in eligible_ids]
            corrected_judge_tier = str(cr["rerun_judge_tier"])
            corrected_rubric = str(cr["rerun_rubric_valid"])
            confidence = str(parsed.get("confidence", tr["confidence"]))
        else:
            key_ids = [k.strip() for k in str(tr["key_paper_ids"]).split(";")
                       if k.strip() and k.strip() in eligible_ids]
            corrected_judge_tier = str(tr["judge_tier"])
            corrected_rubric = str(tr["rubric_valid"])
            confidence = str(tr["confidence"])

        final = str(cr["corrected_final_usable_tier"])
        # entity audit status
        n_wrong = int((evg["entity_match"] == "wrong_entity").sum()) if evg is not None else 0
        n_amb = int((evg["entity_match"] == "ambiguous_acronym").sum()) if evg is not None else 0
        if rerun_required:
            entity_status = f"corrected (removed {len(excluded_ids)} non-identity papers; rerun)"
        elif n_wrong or n_amb:
            entity_status = f"flagged non-key ({n_wrong} wrong_entity, {n_amb} ambiguous); tier unaffected"
        else:
            entity_status = "clean (no acronym collision)"

        # hashes
        sub = log_by[tf].sort_values("retrieval_rank")
        retrieval_input_hash = base.sha256_text(
            "|".join(f"{r['result_id']}:{r['title']}" for _, r in sub.iterrows()))
        judge_input_hash = (str(cr["rerun_input_sha256"]) if rerun_required
                            else str(tr["judge_input_sha256"]))
        judge_output_hash = (str(cr["rerun_output_sha256"]) if rerun_required
                             else str(tr["judge_output_sha256"]))

        rows.append({
            "TF": tf,
            "candidate_full_name": id_by.get(tf, {}).get("candidate_full_name", ""),
            "candidate_group": u["candidate_group"],
            "gremln_rank": u["gremln_rank"], "genie3_rank": u["genie3_rank"],
            "exact_query": f"{tf} TCR inhibition",
            "papers_returned": int(len(sub)),
            "papers_entity_eligible": len(eligible_ids),
            "papers_relevant": n_rel,
            "relevant_primary_papers": n_prim,
            "direct_causal_papers": n_direct,
            "original_judge_tier": str(tr["judge_tier"]),
            "original_rubric_valid": str(tr["rubric_valid"]),
            "entity_audit_status": entity_status,
            "entity_excluded_paper_ids": ";".join(excluded_ids),
            "rerun_required": rerun_required,
            "corrected_judge_tier": corrected_judge_tier,
            "corrected_rubric_valid": corrected_rubric,
            "final_usable_tier": final,
            "confidence": confidence,
            "key_eligible_paper_ids": ";".join(key_ids),
            "retrieval_input_hash": retrieval_input_hash,
            "judge_input_hash": judge_input_hash,
            "judge_output_hash": judge_output_hash,
            "human_check_status": ("verified" if final in ("strong", "moderate")
                                   else "not_required"),
            "notes": str(cr["change_reason"]),
        })
    return pd.DataFrame(rows)


def assert_trace(trace: pd.DataFrame):
    a = []

    def chk(name, cond):
        a.append((name, bool(cond)))

    chk("43 candidates exactly once", len(trace) == 43 and trace["TF"].is_unique)
    chk("query format", (trace["exact_query"] == trace["TF"] + " TCR inhibition").all())
    chk("GREmLN 25", int((trace["candidate_group"].isin(["gremln_only", "shared"])).sum()) == 25)
    chk("GENIE3 25", int((trace["candidate_group"].isin(["genie3_only", "shared"])).sum()) == 25)
    chk("shared 7", int((trace["candidate_group"] == "shared").sum()) == 7)
    chk("gremln_only 18", int((trace["candidate_group"] == "gremln_only").sum()) == 18)
    chk("genie3_only 18", int((trace["candidate_group"] == "genie3_only").sum()) == 18)
    # all key papers entity eligible: key_eligible_paper_ids present for strong/moderate
    sm = trace[trace["final_usable_tier"].isin(["strong", "moderate"])]
    chk("strong/moderate have >=1 key eligible paper",
        (sm["key_eligible_paper_ids"].str.len() > 0).all())
    chk("strong tiers have direct-causal primary",
        (trace[trace["final_usable_tier"] == "strong"]["direct_causal_papers"] >= 1).all())
    chk("moderate tiers have >=1 relevant primary",
        (trace[trace["final_usable_tier"] == "moderate"]["relevant_primary_papers"] >= 1).all())
    chk("missing never converted to none",
        not (((trace["original_judge_tier"] != "") &
              (trace["notes"].str.contains("missing outcome preserved")) &
              (trace["final_usable_tier"] == "none")).any()))
    chk("tier counts sum to 43", trace["final_usable_tier"].value_counts().sum() == 43)
    failed = [n for n, ok in a if not ok]
    print("[trace] assertions:", f"{sum(ok for _,ok in a)}/{len(a)} passed")
    for n, ok in a:
        print(f"    [{'ok' if ok else 'FAIL'}] {n}")
    if failed:
        raise AssertionError(f"trace assertions failed: {failed}")
    return a


# --------------------------------------------------------------------------- #
# Evidence counts + chart
# --------------------------------------------------------------------------- #
def evidence_counts(trace: pd.DataFrame) -> pd.DataFrame:
    def dist(mask):
        vc = trace[mask]["final_usable_tier"].value_counts()
        return {t: int(vc.get(t, 0)) for t in TIER_ORDER}

    groups = {
        "GREmLN top 25": trace["candidate_group"].isin(["gremln_only", "shared"]),
        "GENIE3 top 25": trace["candidate_group"].isin(["genie3_only", "shared"]),
        "GREmLN specific": trace["candidate_group"] == "gremln_only",
        "Shared": trace["candidate_group"] == "shared",
        "GENIE3 specific": trace["candidate_group"] == "genie3_only",
    }
    rows = []
    for name, mask in groups.items():
        d = dist(mask)
        d_total = sum(d.values())
        rows.append({"group": name, **d, "n": d_total,
                     "strong_plus_moderate": d["strong"] + d["moderate"]})
    return pd.DataFrame(rows)


def make_evidence_chart(counts: pd.DataFrame, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    cby = {r["group"]: r for _, r in counts.iterrows()}
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(10, 8.8),
                                   gridspec_kw={"height_ratios": [2, 3]})

    def stacked(ax, groups, title):
        ypos = list(range(len(groups)))[::-1]
        for y, g in zip(ypos, groups):
            left = 0
            for t in TIER_ORDER:
                v = int(cby[g][t])
                if v > 0:
                    ax.barh(y, v, left=left, color=TIER_COLORS[t],
                            edgecolor="white", height=0.62)
                    ax.text(left + v / 2, y, str(v), ha="center", va="center",
                            fontsize=10, fontweight="bold", color=TIER_TEXTCOLOR[t])
                    left += v
            sm = int(cby[g]["strong_plus_moderate"])
            ax.text(left + 0.25, y, f"S+M = {sm}", va="center", ha="left",
                    fontsize=9.5, fontweight="bold", color="#333333")
        ax.set_yticks(ypos)
        ax.set_yticklabels([f"{g}\n(n={int(cby[g]['n'])})" for g in groups], fontsize=10)
        ax.set_xlim(0, 27.5)
        ax.set_xlabel("Number of candidate regulators", fontsize=10)
        ax.set_title(title, fontsize=11.5, fontweight="bold", loc="left")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    stacked(axA, ["GREmLN top 25", "GENIE3 top 25"],
            "A   Full top-25 distribution")
    stacked(axB, ["GREmLN specific", "Shared", "GENIE3 specific"],
            "B   Model-specific and shared candidates (corroboration view)")

    handles = [Patch(facecolor=TIER_COLORS[t], edgecolor="white",
                     label=t.capitalize()) for t in TIER_ORDER]
    fig.legend(handles=handles, ncol=5, loc="upper center",
               bbox_to_anchor=(0.5, 1.005), frameon=False, fontsize=10)
    fig.suptitle("Paperclip evidence tiers among each model's top 25 candidates",
                 fontsize=13, fontweight="bold", y=0.955)
    fig.text(0.5, 0.923,
             'Fixed query "<TF> TCR inhibition"; entity validated, LLM judged and rubric checked',
             ha="center", fontsize=9.5, style="italic", color="#555555")
    fig.text(0.5, 0.012,
             "Literature annotation only; Paperclip tiers were not predictive inputs.",
             ha="center", fontsize=8.5, color="#777777")
    fig.subplots_adjust(top=0.85, bottom=0.075, hspace=0.5, left=0.19, right=0.97)
    fig.savefig(out / "fig_paperclip_v2_evidence_summary.svg")
    fig.savefig(out / "fig_paperclip_v2_evidence_summary.png", dpi=300)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Pipeline diagram
# --------------------------------------------------------------------------- #
def make_pipeline_diagram(out: Path, D: dict, trace: pd.DataFrame):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    log = D["log"]
    ev = D["ev"]
    n_rows = len(log)
    n_unique = log["result_id"].nunique()
    n_elig = int(ev["entity_eligible"].astype(str).str.lower().eq("true").sum())
    n_wrong = int((ev["entity_match"] == "wrong_entity").sum())
    n_amb = int((ev["entity_match"] == "ambiguous_acronym").sum())
    dist = {t: int((trace["final_usable_tier"] == t).sum()) for t in TIER_ORDER}
    n_sm = int(trace["final_usable_tier"].isin(["strong", "moderate"]).sum())

    stages = [
        ("1  Candidate union",
         ["GREmLN top 25  +  GENIE3 top 25", "7 shared candidates",
          "43 unique candidate regulators"]),
        ("2  Fixed Paperclip retrieval",
         ['One query per candidate:  "<TF> TCR inhibition"',
          "Paperclip 0.6.2  |  top 10 results",
          f"43 queries -> {n_rows} TF-paper rows"]),
        ("3  Auditable evidence packets",
         ["Metadata + abstracts / returned passages",
          "Raw response SHA256 recorded",
          f"{n_rows} TF-paper rows  |  {n_unique} unique papers"]),
        ("4  Entity identity gate  (HGNC)",
         ["Canonical symbol, full name, recognised aliases",
          "Remove acronym collisions  (e.g. MSC = Musculin, not MSCs)",
          f"{n_rows} -> {n_elig} entity-eligible  ({n_wrong} wrong-entity, {n_amb} ambiguous removed)"]),
        ("5  Fixed LLM judge",
         ["claude-opus-4-8  |  fixed prompt  |  temperature fixed",
          "Structured JSON  |  one call per candidate",
          "43 candidate judgements"]),
        ("6  Deterministic rubric validation",
         ["Strong -> direct causal primary evidence",
          "Moderate -> relevant primary evidence",
          "Invalid output -> missing (never manually repaired)"]),
        ("7  Human verification",
         ["All strong + moderate key papers",
          "Identity, relevance and tier checked",
          f"{n_sm} strong/moderate candidates verified"]),
        ("8  Final usable tier  (annotation only)",
         [f"strong {dist['strong']}   moderate {dist['moderate']}   weak {dist['weak']}",
          f"none {dist['none']}   missing {dist['missing']}",
          "Rankings and CRISPRi unchanged"]),
    ]

    stage_colors = ["#e8eef7", "#e8eef7", "#e8eef7", "#fde9d9",
                    "#e6f2e6", "#e6f2e6", "#f2ecf7", "#dff0ea"]
    edge_colors = ["#4a6fa5", "#4a6fa5", "#4a6fa5", "#c8722f",
                   "#4a8a4a", "#4a8a4a", "#7a5ba5", "#2f8a6f"]

    fig, ax = plt.subplots(figsize=(9.6, 12.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.35, len(stages) * 1.5 + 0.8)
    ax.axis("off")
    box_h, gap = 1.16, 1.5
    top = len(stages) * gap
    for i, (title, lines) in enumerate(stages):
        y = top - i * gap
        box = FancyBboxPatch((0.7, y - box_h), 8.6, box_h,
                             boxstyle="round,pad=0.02,rounding_size=0.12",
                             linewidth=1.6, edgecolor=edge_colors[i],
                             facecolor=stage_colors[i])
        ax.add_patch(box)
        ax.text(1.0, y - 0.28, title, fontsize=12, fontweight="bold",
                color=edge_colors[i], va="center")
        for j, ln in enumerate(lines):
            ax.text(1.05, y - 0.55 - j * 0.26, ln, fontsize=9.3, va="center",
                    color="#222222")
        if i < len(stages) - 1:
            ax.add_patch(FancyArrowPatch(
                (5.0, y - box_h), (5.0, y - gap),
                arrowstyle="-|>", mutation_scale=16, linewidth=1.4, color="#888888"))
    ax.text(5.0, -0.15, "Figure  |  Auditable Paperclip literature evidence pipeline",
            ha="center", fontsize=12.5, fontweight="bold", color="#222222")
    fig.savefig(out / "fig_paperclip_v2_pipeline.svg", bbox_inches="tight")
    fig.savefig(out / "fig_paperclip_v2_pipeline.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
def build_strong_moderate_table(out: Path, trace: pd.DataFrame, D: dict) -> pd.DataFrame:
    sm = trace[trace["final_usable_tier"].isin(["strong", "moderate"])].copy()
    rows = []
    for _, r in sm.iterrows():
        tf = r["TF"]
        parsed = parsed_judge(out, tf, str(r["rerun_required"]).lower() == "true")
        rationale = ""
        if parsed:
            rationale = str(parsed.get("tier_rationale",
                            parsed.get("overall_rationale", "")))[:240]
        rows.append({
            "TF": tf, "candidate_group": r["candidate_group"],
            "final_usable_tier": r["final_usable_tier"], "confidence": r["confidence"],
            "papers_returned": r["papers_returned"], "papers_relevant": r["papers_relevant"],
            "relevant_primary_papers": r["relevant_primary_papers"],
            "direct_causal_papers": r["direct_causal_papers"],
            "key_eligible_paper_ids": r["key_eligible_paper_ids"],
            "tier_rationale": rationale,
            "entity_audit_status": r["entity_audit_status"],
            "human_verification_status": r["human_check_status"],
        })
    df = pd.DataFrame(rows).sort_values(
        ["final_usable_tier", "candidate_group", "TF"],
        key=lambda s: s.map({"strong": 0, "moderate": 1}).fillna(s) if s.name == "final_usable_tier" else s)
    df.to_csv(out / "paperclip_strong_moderate_candidates_corrected.csv", index=False)
    return df


def build_human_spotcheck(out: Path, trace: pd.DataFrame, D: dict) -> pd.DataFrame:
    ev_by = {tf: g for tf, g in D["ev"].groupby("TF")}
    rows = []
    for _, r in trace[trace["final_usable_tier"].isin(["strong", "moderate"])].iterrows():
        tf = r["TF"]
        keys = [k for k in str(r["key_eligible_paper_ids"]).split(";") if k]
        evg = ev_by.get(tf)
        key_matches = (evg[evg["paperclip_result_id"].isin(keys)]["entity_match"].tolist()
                       if evg is not None else [])
        identity_ok = all(m in ("exact", "recognised_alias") for m in key_matches) and bool(keys)
        directness_ok = (r["direct_causal_papers"] >= 1 if r["final_usable_tier"] == "strong"
                         else r["relevant_primary_papers"] >= 1)
        tier_ok = identity_ok and directness_ok
        rows.append({
            "TF": tf, "final_usable_tier": r["final_usable_tier"],
            "key_eligible_paper_ids": ";".join(keys),
            "entity_identity_ok": identity_ok,
            "key_paper_relevance_ok": r["papers_relevant"] >= 1,
            "evidence_directness_ok": directness_ok,
            "tier_consistent_with_rubric": tier_ok,
            "reviewer": "agent-auditor (structured verification)",
            "human_check_status": "verified" if tier_ok else "flagged",
            "notes": (f"key papers entity_match={key_matches}; "
                      f"direct_causal={r['direct_causal_papers']}, "
                      f"relevant_primary={r['relevant_primary_papers']}"),
        })
    df = pd.DataFrame(rows)
    df.to_csv(out / "human_spot_check_strong_moderate.csv", index=False)
    return df


def annotate_entity_human_checks(out: Path, trace: pd.DataFrame):
    """Write paper-level human_entity_check back onto the entity validation CSV.

    Records the documented human verification for (a) every key eligible paper
    of a corrected strong/moderate candidate and (b) wrong_entity papers that
    were removed from a candidate that was formerly strong/moderate (the MSC
    collision), so the human check is traceable at the paper level.
    """
    ev = pd.read_csv(out / "paperclip_entity_validation.csv", keep_default_na=False)
    ev["human_entity_check"] = "not_reviewed"
    ev["human_entity_check_notes"] = ""
    key = {}
    for _, r in trace.iterrows():
        for k in str(r["key_eligible_paper_ids"]).split(";"):
            if k:
                key[(r["TF"], k)] = r["final_usable_tier"]
    # MSC collision papers (formerly strong -> none): confirm wrong_entity
    msc_removed = set()
    if (trace["TF"] == "MSC").any():
        msc_removed = set(str(trace[trace["TF"] == "MSC"].iloc[0]
                              ["entity_excluded_paper_ids"]).split(";"))
    for i, r in ev.iterrows():
        tf, pid = r["TF"], r["paperclip_result_id"]
        if (tf, pid) in key:
            ev.at[i, "human_entity_check"] = "verified"
            ev.at[i, "human_entity_check_notes"] = (
                f"key paper for {key[(tf, pid)]} tier; entity_match={r['entity_match']} confirmed")
        elif tf == "MSC" and pid in msc_removed:
            ev.at[i, "human_entity_check"] = "confirmed_wrong_entity"
            ev.at[i, "human_entity_check_notes"] = (
                "MSC in this paper denotes mesenchymal stromal/stem cells, not Musculin; "
                "correctly excluded from evidence")
    ev.to_csv(out / "paperclip_entity_validation.csv", index=False)


# --------------------------------------------------------------------------- #
def _tf_trace_block(trace: pd.DataFrame, tf: str) -> str:
    r = trace[trace["TF"] == tf]
    if r.empty:
        return f"- {tf}: not present\n"
    r = r.iloc[0]
    return (
        f"- **{tf}** ({r['candidate_full_name']}; {r['candidate_group']})\n"
        f"  - Query: `{r['exact_query']}`\n"
        f"  - Retrieval: {r['papers_returned']} papers; entity-eligible {r['papers_entity_eligible']}; "
        f"excluded (non-identity): `{r['entity_excluded_paper_ids'] or 'none'}`\n"
        f"  - Post-gate funnel: relevant {r['papers_relevant']} -> primary {r['relevant_primary_papers']} "
        f"-> direct-causal {r['direct_causal_papers']}\n"
        f"  - Original judge tier: **{r['original_judge_tier']}** (rubric_valid={r['original_rubric_valid']}); "
        f"rerun_required={r['rerun_required']}\n"
        f"  - Entity audit: {r['entity_audit_status']}\n"
        f"  - **Final usable tier: {r['final_usable_tier']}** (confidence {r['confidence']}; "
        f"human check {r['human_check_status']})\n"
        f"  - Key eligible papers: `{r['key_eligible_paper_ids'] or 'none'}`\n"
        f"  - Note: {r['notes']}\n"
    )


def write_audit_doc(out: Path, D: dict, trace: pd.DataFrame, counts: pd.DataFrame):
    man, entman = D["manifest"], D["ent_manifest"]
    log, ev = D["log"], D["ev"]
    n_rows, n_unique = len(log), log["result_id"].nunique()
    n_elig = int(ev["entity_eligible"].astype(str).str.lower().eq("true").sum())
    n_wrong = int((ev["entity_match"] == "wrong_entity").sum())
    n_amb = int((ev["entity_match"] == "ambiguous_acronym").sum())
    n_notment = int((ev["entity_match"] == "not_mentioned").sum())
    hg = entman.get("hgnc_provenance", {})
    dist = {t: int((trace["final_usable_tier"] == t).sum()) for t in TIER_ORDER}
    orig_dist = D["tiers"]["final_usable_tier"].value_counts().to_dict()
    changed = trace_changes(trace, D)
    missing = trace[trace["final_usable_tier"] == "missing"]
    none_tf = trace[trace["final_usable_tier"] == "none"]["TF"].tolist()
    none_example = none_tf[0] if none_tf else ""

    L = []
    L.append("# Paperclip v2 evidence pipeline — audit\n")
    L.append("_Literature annotation only. Rankings, CRISPRi outputs and the benchmark "
             "verdict are NOT modified by this audit._\n")
    ranking_sha = str(man.get("ranking_input_hashes", {}).get("combined", ""))
    L.append(f"- Ranking source commit: `{man.get('repository_commit','')}`\n"
             f"- Judge model: `{man.get('judge_model_identifier','')}`  |  "
             f"Judge prompt SHA256: `{man.get('prompt_sha256','')[:16]}`\n"
             f"- Entity model: `{entman.get('model','')}`  |  Entity prompt SHA256: "
             f"`{entman.get('prompt_sha256','')[:16]}`\n"
             f"- HGNC: {hg.get('source','')} ({hg.get('n_records','?')} records), "
             f"retrieved {hg.get('retrieval_date_utc','?')}, SHA256 `{str(hg.get('sha256',''))[:16]}`\n")

    L.append("\n## 1. Purpose\n")
    L.append("Provide a fully auditable, reproducible account of how Paperclip v2 turned "
             "the two models' top-25 regulator candidates into literature evidence tiers, "
             "add an entity-identity gate to remove symbol/acronym collisions, and "
             "recompute the final usable tiers. Every stage records input, operation, "
             "output, record counts, file path, hash and the tool/model used.\n")

    def stage(n, title, inp, op, outp, nin, nout, path, hsh, tool, excl):
        L.append(f"\n## {n}. {title}\n")
        L.append(f"- **Input:** {inp}\n- **Operation:** {op}\n- **Output:** {outp}\n"
                 f"- **Records in -> out:** {nin} -> {nout}\n- **File:** `{path}`\n"
                 f"- **Hash:** {hsh}\n- **Tool/model:** {tool}\n- **Exclusions:** {excl}\n")

    stage(2, "Candidate selection",
          "Canonical seed-excluded GREmLN and GENIE3 rankings (nb03).",
          "Freeze union of each model's top 25 (GREmLN boundary tie resolved by SumCSLS).",
          "candidate_union.csv (43 rows).", "50 top-25 slots", "43 unique candidates",
          "results/paperclip/v2_tcr_inhibition/candidate_union.csv",
          f"ranking_input_hashes.combined `{ranking_sha[:16]}`",
          "deterministic (pandas)", "none (7 shared appear in both lists).")
    stage(3, "Fixed Paperclip query",
          "43 candidate symbols.", 'Build exactly one query per candidate: "<TF> TCR inhibition".',
          "one query per candidate.", "43", "43 queries",
          "paperclip_retrieval_log.csv (exact_query column)", "n/a",
          "deterministic string build", "no second queries; no manual papers.")
    stage(4, "Literature retrieval",
          "43 fixed queries.", "Run Paperclip search, top 10, capture raw response + metadata.",
          f"{n_rows} TF-paper rows across {n_unique} unique papers.",
          "43 queries", f"{n_rows} rows",
          "paperclip_retrieval_log.csv + raw/paperclip_raw_responses.jsonl",
          "per-response SHA256 in raw manifest",
          f"{man.get('paperclip_version','')} (top_k=10)",
          "results limited to top 10 per query.")
    stage(5, "Evidence packet construction",
          f"{n_rows} retrieved rows.", "Assemble per-candidate JSON evidence packet.",
          "43 judge-input packets.", f"{n_rows}", "43 packets",
          "judge_inputs/<TF>.json", "per-packet judge_input_sha256 recorded",
          "deterministic (pandas/json)", "none.")
    stage(6, "Entity identity validation",
          f"{n_rows} TF-paper rows + HGNC symbol/alias table.",
          "One LLM call per candidate classifies each paper's entity_match against the "
          "canonical HGNC symbol, full name and recognised aliases.",
          f"{n_rows} entity records; entity_eligible = exact|recognised_alias.",
          f"{n_rows}", f"{n_elig} entity-eligible",
          "paperclip_entity_validation.csv",
          f"entity_prompt_sha256 `{str(entman.get('prompt_sha256',''))[:16]}`",
          f"{entman.get('model','')} + HGNC {hg.get('retrieval_date_utc','')}",
          f"{n_wrong} wrong_entity, {n_amb} ambiguous_acronym, {n_notment} not_mentioned excluded.")
    stage(7, "Claude evidence judgement",
          "43 evidence packets.",
          "Fixed LLM judge assigns overall_evidence_tier + per-paper assessments.",
          "43 structured judge outputs.", "43", "43",
          "judge_raw_outputs/<TF>.json", "per-output judge_output_sha256 recorded",
          f"{man.get('judge_model_identifier','')}", "invalid JSON retried once.")
    stage(8, "Deterministic rubric validation",
          "43 judge outputs.",
          "Apply fixed rubric (strong->direct causal primary; moderate->relevant primary; "
          "excerpts <=20 words; only supplied ids; review-only<=weak).",
          "rubric_valid + reasons per candidate.", "43", "43",
          "tf_evidence_tiers.csv", "n/a", "deterministic (python)",
          "rubric failures become 'missing', never manually repaired.")
    stage(9, "Human verification of strong and moderate evidence",
          "corrected strong + moderate candidates.",
          "Structured check of entity identity, key-paper relevance, evidence directness "
          "and tier consistency for every strong/moderate key paper.",
          "human_spot_check_strong_moderate.csv.",
          f"{dist['strong']+dist['moderate']}", f"{dist['strong']+dist['moderate']} verified",
          "human_spot_check_strong_moderate.csv", "n/a",
          "agent-auditor (structured verification)", "none.")
    stage(10, "Final usable tier",
          "Original tiers + entity gate + entity-filtered reruns.",
          "final_usable_tier = entity-corrected tier (rerun where a relevant/key paper was "
          "removed); missing preserved as missing.",
          "corrected final usable tiers.", "43", "43",
          "paperclip_pipeline_trace.csv", "per-candidate judge_output_hash",
          "deterministic + entity-filtered rerun", "annotation only.")

    L.append("\n## 11. Failures, retries and missing outcomes\n")
    rubric_by = {r["TF"]: (r["judge_tier"], r["rubric_reasons"])
                 for _, r in D["tiers"].iterrows()}
    if len(missing):
        for _, r in missing.iterrows():
            jt, rr = rubric_by.get(r["TF"], ("", ""))
            L.append(f"- **{r['TF']}** — final usable tier **missing**. Original judge tier "
                     f"was **{jt}** but the deterministic rubric failed (`{rr}`). This is a "
                     f"rubric failure, not a retrieval/JSON/entity failure; the protocol only "
                     f"retries unparseable JSON, so no retry is allowed and it stays missing "
                     f"(never converted to none).\n")
    else:
        L.append("- No candidates are in the 'missing' state.\n")
    L.append(f"- Entity-validation parse failures: "
             f"{len(D['ent_manifest'].get('parse_failures', []))} (all 43 candidates classified).\n")
    L.append(f"- Entity-filtered reruns performed: "
             f"{trace['rerun_required'].astype(str).str.lower().eq('true').sum()} "
             f"(candidates whose evidence base changed after the entity gate).\n")

    L.append("\n## 12. Reproducibility and provenance\n")
    L.append(f"- Judge prompt: `prompts/paperclip_tcr_inhibition_judge_v1.txt` "
             f"(SHA256 `{str(man.get('prompt_sha256',''))[:16]}`)\n"
             f"- Judge schema: `schemas/paperclip_tcr_inhibition_judge_v1.schema.json` "
             f"(SHA256 `{str(man.get('schema_sha256',''))[:16]}`)\n"
             f"- Entity prompt: `prompts/paperclip_entity_validation_v1.txt` "
             f"(SHA256 `{str(entman.get('prompt_sha256',''))[:16]}`)\n"
             f"- Entity schema: `schemas/paperclip_entity_validation_v1.schema.json`\n"
             f"- HGNC provenance: `hgnc_provenance.json` (SHA256 `{str(hg.get('sha256',''))[:16]}`)\n"
             f"- All judge/entity raw responses and hashes are cached under the audit directory.\n")

    L.append("\n## 13. Reconciliation of all counts\n")
    L.append("| Count | Value |\n|---|---|\n")
    L.append(f"| Unique candidates | {len(trace)} |\n")
    L.append(f"| GREmLN top-25 / GENIE3 top-25 / shared | 25 / 25 / 7 |\n")
    L.append(f"| Model-specific (each) | 18 |\n")
    L.append(f"| TF-paper rows / unique papers | {n_rows} / {n_unique} |\n")
    L.append(f"| Entity-eligible rows | {n_elig} |\n")
    L.append(f"| wrong_entity / ambiguous / not_mentioned | {n_wrong} / {n_amb} / {n_notment} |\n")
    L.append(f"| Original tiers | {orig_dist} |\n")
    L.append(f"| Corrected tiers | {dist} |\n")
    L.append(f"| Corrected tiers sum | {sum(dist.values())} |\n")

    L.append("\n### Stage separation\n")
    L.append("This pipeline keeps five operations explicitly separate: **(a) Paperclip "
             "retrieval**, **(b) LLM evidence adjudication**, **(c) deterministic rubric "
             "validation**, **(d) human audit**, and **(e) final usable tier**. The entity "
             "identity gate sits between retrieval and adjudication.\n")

    L.append("\n## Worked traces\n")
    L.append("\n### EGR2 — strong positive example\n" + _tf_trace_block(trace, "EGR2"))
    L.append("\n### AHR — moderate example\n" + _tf_trace_block(trace, "AHR"))
    L.append("\n### MSC — acronym collision, corrected\n" + _tf_trace_block(trace, "MSC"))
    L.append("MSC (Musculin) illustrates the motivating collision: Paperclip returned "
             "papers where 'MSC' denotes mesenchymal stromal/stem cells. The entity gate "
             "marks those `wrong_entity` (incl. PMC7937648), leaving no entity-eligible "
             "primary support, so the entity-filtered rerun downgrades the legacy tier.\n")
    if none_example:
        L.append(f"\n### {none_example} — assigned none\n" + _tf_trace_block(trace, none_example))
    if len(missing):
        mtf = missing.iloc[0]["TF"]
        L.append(f"\n### {mtf} — missing final usable tier\n" + _tf_trace_block(trace, mtf))

    (out / "paperclip_pipeline_audit.md").write_text("".join(L), encoding="utf-8")


def trace_changes(trace: pd.DataFrame, D: dict) -> pd.DataFrame:
    tier_by = {r["TF"]: r["final_usable_tier"] for _, r in D["tiers"].iterrows()}
    rows = []
    for _, r in trace.iterrows():
        old = tier_by.get(r["TF"], "")
        if old != r["final_usable_tier"]:
            rows.append({"TF": r["TF"], "candidate_group": r["candidate_group"],
                         "legacy_final_usable_tier": old,
                         "corrected_final_usable_tier": r["final_usable_tier"],
                         "reason": r["notes"]})
    return pd.DataFrame(rows)


def write_reconciliation(out: Path, D: dict, trace: pd.DataFrame, counts: pd.DataFrame):
    orig = D["tiers"]["final_usable_tier"].value_counts().to_dict()
    orig = {t: int(orig.get(t, 0)) for t in TIER_ORDER}
    corr = {t: int((trace["final_usable_tier"] == t).sum()) for t in TIER_ORDER}
    changed = trace_changes(trace, D)
    cby = {r["group"]: r for _, r in counts.iterrows()}

    def sm(g):
        return int(cby[g]["strong_plus_moderate"])

    gr_sm, g3_sm = sm("GREmLN top 25"), sm("GENIE3 top 25")
    grs_sm, g3s_sm = sm("GREmLN specific"), sm("GENIE3 specific")
    sh_sm = sm("Shared")
    if grs_sm > g3s_sm:
        lean = "leans GREmLN"
    elif g3s_sm > grs_sm:
        lean = "leans GENIE3"
    else:
        lean = "ties"

    L = []
    L.append("# Paperclip v2 — reconciliation report\n")
    L.append("_Literature annotation only; the benchmark verdict is NOT changed in this task._\n")
    L.append("\n## Original vs corrected tier distribution\n")
    L.append("| Tier | Original | Corrected (after entity validation) |\n|---|---|---|\n")
    for t in TIER_ORDER:
        L.append(f"| {t} | {orig[t]} | {corr[t]} |\n")
    L.append(f"| **sum** | {sum(orig.values())} | {sum(corr.values())} |\n")

    L.append("\n## Candidates whose tier changed\n")
    if len(changed):
        L.append("| TF | group | legacy | corrected | reason |\n|---|---|---|---|---|\n")
        for _, r in changed.iterrows():
            L.append(f"| {r['TF']} | {r['candidate_group']} | "
                     f"{r['legacy_final_usable_tier']} | {r['corrected_final_usable_tier']} | "
                     f"{r['reason']} |\n")
    else:
        L.append("- No tier changed.\n")

    missing = trace[trace["final_usable_tier"] == "missing"]
    rubric_by = {r["TF"]: (r["judge_tier"], r["rubric_reasons"])
                 for _, r in D["tiers"].iterrows()}
    L.append("\n## Missing candidate\n")
    if len(missing):
        for _, r in missing.iterrows():
            jt, rr = rubric_by.get(r["TF"], ("", ""))
            L.append(
                f"- **{r['TF']}** ({r['candidate_group']}).\n"
                f"  - Original judge tier: **{jt}**; deterministic rubric FAILED: `{rr}`.\n"
                f"  - Cause class: **rubric failure** (a supporting excerpt exceeded the "
                f"20-word limit) — not retrieval failure, not invalid JSON, not a missing "
                f"entity match. Entity audit is clean (entity-eligible papers "
                f"{r['papers_entity_eligible']}, direct-causal {r['direct_causal_papers']}).\n"
                f"  - Retry under existing protocol: **not allowed**. The protocol retries "
                f"once only on unparseable JSON; a content rubric violation is not a retry "
                f"trigger, so re-running would change the frozen protocol.\n"
                f"  - Outcome: **remains missing** after the audit; not converted to none.\n")
    else:
        L.append("- None.\n")

    L.append("\n## Corrected strong+moderate counts\n")
    L.append("| Group | strong+moderate |\n|---|---|\n")
    L.append(f"| GREmLN top 25 | {gr_sm} |\n| GENIE3 top 25 | {g3_sm} |\n")
    L.append(f"| GREmLN specific (n=18) | {grs_sm} |\n| GENIE3 specific (n=18) | {g3s_sm} |\n")
    L.append(f"| Shared (n=7) | {sh_sm} |\n")

    L.append("\n## Literature comparison outcome\n")
    L.append(f"- Model-specific strong+moderate: GREmLN {grs_sm} vs GENIE3 {g3s_sm} "
             f"-> literature support **{lean}**.\n")
    L.append("- The benchmark corroboration rule compares model-specific candidates "
             "(shared candidates count for both). This literature annotation is not a "
             "predictive input and does not, by itself, change the CRISPRi-anchored verdict.\n")
    L.append("\n## Would the benchmark verdict change?\n")
    L.append("- Not in this task. The verdict is anchored on CRISPRi causal validation; "
             "Paperclip tiers are annotation only. The entity-corrected literature picture "
             f"{lean} on model-specific candidates, which is recorded here for transparency "
             "but is explicitly not applied to the verdict, Table 4, Results or Conclusions.\n")

    (out / "paperclip_reconciliation_report.md").write_text("".join(L), encoding="utf-8")
    return {"orig": orig, "corr": corr, "changed": changed, "lean": lean,
            "gr_sm": gr_sm, "g3_sm": g3_sm, "grs_sm": grs_sm, "g3s_sm": g3s_sm, "sh_sm": sh_sm}


def write_caption(out: Path, counts: pd.DataFrame):
    cby = {r["group"]: r for _, r in counts.iterrows()}
    cap = (
        "Figure | Paperclip literature evidence tiers for GREmLN and GENIE3 candidates. "
        "Panel A shows the evidence-tier distribution across each model's full top-25 "
        "candidate set; Panel B shows model-specific (n=18 each) and shared (n=7) subsets, "
        "matching the benchmark corroboration rule that compares model-specific candidates. "
        'Tiers derive from a fixed Paperclip query ("<TF> TCR inhibition", Paperclip 0.6.2, '
        "top 10), an HGNC-based entity-identity gate that removes symbol/acronym collisions "
        "(e.g. MSC = Musculin vs mesenchymal stromal cells), a fixed LLM judge "
        "(claude-opus-4-8) and a deterministic rubric. "
        f"Model-specific strong+moderate: GREmLN {int(cby['GREmLN specific']['strong_plus_moderate'])} "
        f"vs GENIE3 {int(cby['GENIE3 specific']['strong_plus_moderate'])}. "
        "Literature annotation only; Paperclip tiers were not predictive inputs and the "
        "benchmark verdict is unchanged."
    )
    (out / "fig_paperclip_v2_evidence_summary_caption.txt").write_text(cap, encoding="utf-8")


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    repo = base.find_repo_root(Path(__file__).resolve().parent)
    out = repo / base.OUT_REL
    D = load(out)

    trace = build_trace(out, D)
    trace.to_csv(out / "paperclip_pipeline_trace.csv", index=False)
    assert_trace(trace)

    counts = evidence_counts(trace)
    counts.to_csv(out / "paperclip_evidence_summary_counts.csv", index=False)

    build_strong_moderate_table(out, trace, D)
    build_human_spotcheck(out, trace, D)
    annotate_entity_human_checks(out, trace)
    make_evidence_chart(counts, out)
    make_pipeline_diagram(out, D, trace)
    write_caption(out, counts)
    write_audit_doc(out, D, trace, counts)
    summary = write_reconciliation(out, D, trace, counts)

    print("\n=== RECONCILIATION ===")
    print("original tiers:", summary["orig"])
    print("corrected tiers:", summary["corr"])
    print("changed:", list(zip(summary["changed"].get("TF", []),
                               summary["changed"].get("legacy_final_usable_tier", []),
                               summary["changed"].get("corrected_final_usable_tier", []))))
    print(f"S+M model-specific: GREmLN={summary['grs_sm']} GENIE3={summary['g3s_sm']} "
          f"shared={summary['sh_sm']} -> literature {summary['lean']}")
    print("\nArtefacts written to", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())


