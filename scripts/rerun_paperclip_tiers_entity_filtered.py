#!/usr/bin/env python3
"""Recompute Paperclip v2 evidence tiers after the entity-identity gate.

For every candidate we determine which retrieved papers are entity_eligible
(entity_match in {exact, recognised_alias}). A candidate's evidence base only
changes when the entity gate removes a paper that the original judge had marked
relevant=true (or cited as a key paper). For those candidates we RE-RUN the
same fixed evidence judge (claude-opus-4-8, identical prompt) on a packet that
contains only the entity-eligible papers, then re-apply the deterministic
rubric. All other candidates keep their original tier unchanged.

The special 'missing' outcome (rubric failure under the original protocol) is
never rerun here and is never converted to 'none'.

Outputs:
    results/paperclip/v2_tcr_inhibition/judge_rerun_outputs/<TF>.json
    results/paperclip/v2_tcr_inhibition/corrected_tiers.csv
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_paperclip_tcr_inhibition_audit as base  # noqa: E402

ELIGIBLE = {"exact", "recognised_alias"}
RERUNNABLE = {"strong", "moderate", "weak"}  # 'missing'/'none' are not rerun


def load_parsed_judge(rec_path: Path):
    if not rec_path.exists():
        return None
    rec = json.loads(rec_path.read_text())
    return rec.get("parsed")


def main(argv=None) -> int:
    import pandas as pd

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo-root", default=None)
    ap.add_argument("--model", default=base.DEFAULT_JUDGE_MODEL)
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--plan-only", action="store_true",
                    help="only print which candidates would be rerun")
    args = ap.parse_args(argv)

    repo = (Path(args.repo_root).resolve() if args.repo_root
            else base.find_repo_root(Path(__file__).resolve().parent))
    out = repo / base.OUT_REL
    rerun_dir = out / "judge_rerun_outputs"
    rerun_dir.mkdir(exist_ok=True)

    union = pd.read_csv(out / "candidate_union.csv")
    log = pd.read_csv(out / "paperclip_retrieval_log.csv", keep_default_na=False)
    tiers = pd.read_csv(out / "tf_evidence_tiers.csv", keep_default_na=False)
    ev = pd.read_csv(out / "paperclip_entity_validation.csv", keep_default_na=False)

    prompt_template = (repo / base.PROMPT_REL).read_text(encoding="utf-8")
    schema = json.loads((repo / base.SCHEMA_REL).read_text(encoding="utf-8"))

    elig_ids = {tf: set(g[g["entity_eligible"].astype(str).str.lower() == "true"]["paperclip_result_id"])
                for tf, g in ev.groupby("TF")}
    all_ids = {tf: set(g["paperclip_result_id"]) for tf, g in ev.groupby("TF")}
    tier_by = {r["TF"]: r for _, r in tiers.iterrows()}
    jraw = out / "judge_raw_outputs"

    # ---- decide rerun set ----
    plan = []
    for tf in union["TF"]:
        row = tier_by[tf]
        legacy_final = str(row["final_usable_tier"])
        legacy_judge = str(row["judge_tier"])
        parsed = load_parsed_judge(jraw / f"{tf}.json")
        eligible = elig_ids.get(tf, set())
        ineligible = all_ids.get(tf, set()) - eligible
        # which ineligible papers did the judge treat as relevant / key?
        removed_relevant, removed_key = [], []
        if parsed:
            for a in parsed.get("paper_assessments", []) or []:
                rid = a.get("paperclip_result_id")
                if rid in ineligible and a.get("relevant"):
                    removed_relevant.append(rid)
            for rid in parsed.get("key_paperclip_result_ids", []) or []:
                if rid in ineligible:
                    removed_key.append(rid)
        changed = bool(removed_relevant or removed_key)
        rerun = changed and legacy_final in RERUNNABLE
        plan.append({
            "TF": tf, "legacy_judge": legacy_judge, "legacy_final": legacy_final,
            "eligible": eligible, "removed_relevant": removed_relevant,
            "removed_key": removed_key, "rerun": rerun,
        })

    rerun_tfs = [p["TF"] for p in plan if p["rerun"]]
    print(f"[rerun] candidates needing rerun ({len(rerun_tfs)}): {rerun_tfs}")
    for p in plan:
        if p["removed_relevant"] or p["removed_key"]:
            print(f"   {p['TF']}: legacy={p['legacy_final']} removed_relevant={p['removed_relevant']} "
                  f"removed_key={p['removed_key']} -> rerun={p['rerun']}")
    if args.plan_only:
        return 0

    api_key = os.environ.get("CURSOR_API_KEY")
    if rerun_tfs and not api_key:
        print("[rerun] CURSOR_API_KEY not set -> cannot rerun judge.")
        return 2

    log_by = {tf: g.sort_values("retrieval_rank") for tf, g in log.groupby("TF")}
    results = []
    for p in plan:
        tf = p["TF"]
        base_row = {
            "TF": tf,
            "candidate_group": tier_by[tf]["candidate_group"],
            "legacy_judge_tier": p["legacy_judge"],
            "legacy_final_usable_tier": p["legacy_final"],
            "entity_ineligible_removed_relevant": ";".join(p["removed_relevant"]),
            "entity_ineligible_removed_key": ";".join(p["removed_key"]),
            "rerun_required": p["rerun"],
        }
        if not p["rerun"]:
            base_row.update({
                "rerun_judge_tier": "", "rerun_rubric_valid": "",
                "corrected_final_usable_tier": p["legacy_final"],
                "rerun_input_sha256": "", "rerun_output_sha256": "",
                "change_reason": ("unchanged: missing outcome preserved (not rerun)"
                                  if p["legacy_final"] == "missing"
                                  else "unchanged: entity gate removed no relevant/key paper"),
            })
            results.append(base_row)
            continue

        # build entity-filtered packet (mirrors stage_inputs, eligible papers only)
        sub = log_by[tf]
        sub = sub[sub["result_id"].isin(p["eligible"])]
        papers = [{
            "paperclip_result_id": r["result_id"], "rank": int(r["retrieval_rank"]),
            "title": r["title"], "year": r["year"],
            "publication_type": r["publication_type"], "abstract": r["abstract"],
            "relevant_passages": r["snippet"], "doi_or_pmcid": r["doi"] or r["pmcid"],
        } for _, r in sub.iterrows()]
        packet = {"tf": tf, "query": base.build_query(tf), "top_k": base.TOP_K,
                  "papers_retrieved": len(papers), "papers": papers,
                  "entity_filtered": True,
                  "excluded_result_ids": sorted(all_ids[tf] - p["eligible"])}
        packet_text = json.dumps(packet, indent=2, ensure_ascii=False)
        in_sha = base.sha256_text(packet_text)
        prompt = (base._render_prompt(prompt_template, tf, len(papers))
                  + "\n\nRETRIEVED PAPERS (JSON):\n" + packet_text)

        cache_file = rerun_dir / f"{tf}.json"
        if (not args.no_cache) and cache_file.exists():
            record = json.loads(cache_file.read_text())
            print(f"[rerun] {tf}: cached")
        else:
            parsed, attempts = None, []
            for attempt in (1, 2):
                extra = "" if attempt == 1 else (
                    "\n\nYour previous response was not valid JSON. Return ONLY the JSON object.")
                try:
                    text, run_id, agent_id, _, capture = base._call_cursor_judge(
                        prompt + extra, args.model, api_key)
                except Exception as e:
                    attempts.append({"attempt": attempt, "startup_error": str(e)})
                    break
                obj = base.extract_json(text)
                attempts.append({"attempt": attempt, "run_id": run_id, "agent_id": agent_id,
                                 "capture_method": capture, "raw": text, "parsed_ok": obj is not None})
                if obj is not None:
                    parsed = obj
                    break
            record = {"tf": tf, "model": args.model, "timestamp_utc": base.utcnow(),
                      "entity_filtered": True, "judge_input_sha256": in_sha,
                      "excluded_result_ids": packet["excluded_result_ids"],
                      "attempts": attempts, "parsed": parsed}
            cache_file.write_text(json.dumps(record, indent=2, ensure_ascii=False))
            print(f"[rerun] {tf}: {'ok' if parsed else 'PARSE FAIL'} "
                  f"(eligible papers={len(papers)})")

        parsed = record.get("parsed")
        out_sha = base.sha256_text(json.dumps(record, sort_keys=True))
        eligible_ids_list = list(p["eligible"])
        rubric_valid, reasons = base.run_rubric_checks(
            tf, parsed, eligible_ids_list, base.build_query(tf), len(papers))
        rerun_tier = str(parsed.get("overall_evidence_tier", "")).lower() if parsed else ""
        if parsed is None:
            corrected = "missing"
            reason = "rerun produced unparseable JSON -> missing (protocol)"
        elif rubric_valid:
            corrected = rerun_tier
            reason = f"entity-filtered rerun tier={rerun_tier}; rubric passed"
        else:
            corrected = "missing"
            reason = f"entity-filtered rerun failed rubric: {'; '.join(reasons)[:160]}"
        base_row.update({
            "rerun_judge_tier": rerun_tier, "rerun_rubric_valid": rubric_valid,
            "corrected_final_usable_tier": corrected,
            "rerun_input_sha256": in_sha, "rerun_output_sha256": out_sha,
            "change_reason": reason,
        })
        results.append(base_row)

    cdf = pd.DataFrame(results)
    cdf.to_csv(out / "corrected_tiers.csv", index=False)
    changed = cdf[cdf["legacy_final_usable_tier"] != cdf["corrected_final_usable_tier"]]
    print(f"\n[rerun] wrote corrected_tiers.csv ({len(cdf)} rows); "
          f"tiers changed for {len(changed)}: "
          f"{list(zip(changed['TF'], changed['legacy_final_usable_tier'], changed['corrected_final_usable_tier']))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
