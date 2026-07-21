#!/usr/bin/env python3
"""Entity-identity validation gate for the Paperclip v2 evidence pipeline.

For every candidate regulator, one fixed LLM call (claude-opus-4-8, via the
Cursor SDK) classifies each of its retrieved papers against the regulator's
canonical HGNC symbol, full name and recognised aliases, detecting symbol /
acronym collisions (e.g. MSC = Musculin vs mesenchymal stromal cells).

A paper is entity_eligible only when entity_match in {exact, recognised_alias}.
Ambiguous / wrong-entity / not-mentioned papers are excluded before evidence
tier adjudication.

Outputs:
    results/paperclip/v2_tcr_inhibition/paperclip_entity_validation.csv
    results/paperclip/v2_tcr_inhibition/entity_raw_outputs/<TF>.json
    results/paperclip/v2_tcr_inhibition/entity_run_manifest.json

Reruns reuse cached entity_raw_outputs/<TF>.json unless --no-cache is given.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_paperclip_tcr_inhibition_audit as base  # noqa: E402

PROMPT_REL = "prompts/paperclip_entity_validation_v1.txt"
SCHEMA_REL = "schemas/paperclip_entity_validation_v1.schema.json"
OUT_REL = "results/paperclip/v2_tcr_inhibition"
DEFAULT_MODEL = os.environ.get("PAPERCLIP_JUDGE_MODEL", "claude-opus-4-8")
DRY_TFS = ["MSC", "EGR2", "MAF"]


def render_prompt(template: str, tf: str, full_name: str, aliases: str, n: int) -> str:
    return (template.replace("{TF}", tf)
                    .replace("{FULL_NAME}", full_name or "(not in HGNC)")
                    .replace("{ALIASES}", aliases or "(none)")
                    .replace("{N}", str(n)))


def main(argv=None) -> int:
    import pandas as pd

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo-root", default=None)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--all", action="store_true")
    g.add_argument("--dry", action="store_true")
    g.add_argument("--tfs", default=None)
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--no-cache", action="store_true")
    args = ap.parse_args(argv)

    repo = (Path(args.repo_root).resolve() if args.repo_root
            else base.find_repo_root(Path(__file__).resolve().parent))
    out = repo / OUT_REL
    raw_dir = out / "entity_raw_outputs"
    raw_dir.mkdir(parents=True, exist_ok=True)
    use_cache = not args.no_cache

    union = pd.read_csv(out / "candidate_union.csv")
    ident = pd.read_csv(out / "candidate_identity.csv", keep_default_na=False)
    log = pd.read_csv(out / "paperclip_retrieval_log.csv", keep_default_na=False)
    prompt_template = (repo / PROMPT_REL).read_text(encoding="utf-8")
    schema = json.loads((repo / SCHEMA_REL).read_text(encoding="utf-8"))
    prompt_sha = base.sha256_file(repo / PROMPT_REL)

    id_by = {r["TF"]: r for _, r in ident.iterrows()}
    all_tfs = list(union["TF"])
    if args.tfs:
        tfs = [t.strip() for t in args.tfs.split(",") if t.strip()]
    elif args.dry:
        tfs = [t for t in DRY_TFS if t in all_tfs]
    else:
        tfs = all_tfs

    api_key = os.environ.get("CURSOR_API_KEY")
    if not api_key:
        print("[entity] CURSOR_API_KEY not set -> cannot run entity validation.")
        return 2

    ELIGIBLE = {"exact", "recognised_alias"}
    rows = []
    stats = {"attempted": 0, "completed": 0, "retries": 0, "parse_failures": [],
             "collisions": []}

    for tf in tfs:
        idr = id_by.get(tf, {})
        full_name = str(idr.get("candidate_full_name", "") or "")
        aliases = str(idr.get("recognised_gene_aliases", "") or "")
        sub = log[log["TF"] == tf].sort_values("retrieval_rank")
        papers = [{"paperclip_result_id": r["result_id"], "retrieval_rank": int(r["retrieval_rank"]),
                   "title": r["title"], "abstract": r["abstract"], "returned_passage": r["snippet"]}
                  for _, r in sub.iterrows()]
        packet = {"tf": tf, "candidate_symbol": tf, "candidate_full_name": full_name,
                  "recognised_gene_aliases": aliases, "papers": papers}
        prompt = (render_prompt(prompt_template, tf, full_name, aliases, len(papers))
                  + "\n\nRETRIEVED PAPERS (JSON):\n"
                  + json.dumps(packet, ensure_ascii=False, indent=2))

        cache_file = raw_dir / f"{tf}.json"
        if use_cache and cache_file.exists():
            record = json.loads(cache_file.read_text())
            print(f"[entity] {tf}: cached")
        else:
            stats["attempted"] += 1
            attempts, parsed = [], None
            for attempt in (1, 2):
                extra = "" if attempt == 1 else (
                    "\n\nYour previous response was not valid JSON. Return ONLY the JSON object.")
                try:
                    text, run_id, agent_id, _, capture = base._call_cursor_judge(
                        prompt + extra, args.model, api_key,
                        file_markers=("paper_entities",))
                except Exception as e:
                    attempts.append({"attempt": attempt, "startup_error": str(e)})
                    break
                obj = base.extract_json(text)
                attempts.append({"attempt": attempt, "run_id": run_id, "agent_id": agent_id,
                                 "capture_method": capture, "raw": text, "parsed_ok": obj is not None})
                if obj is not None:
                    parsed = obj
                    break
                stats["retries"] += 1
            record = {"tf": tf, "model": args.model, "timestamp_utc": base.utcnow(),
                      "attempts": attempts, "parsed": parsed}
            cache_file.write_text(json.dumps(record, indent=2, ensure_ascii=False))
            if parsed is None:
                stats["parse_failures"].append(tf)
            print(f"[entity] {tf}: {'ok' if parsed else 'PARSE FAIL'} ({len(papers)} papers)")

        parsed = record.get("parsed")
        if parsed is not None:
            stats["completed"] += 1
            schema_errs = base.validate_judge_json(parsed, schema)
            ent_by_id = {p.get("paperclip_result_id"): p
                         for p in parsed.get("paper_entities", [])}
            # Fallback: the model occasionally copies a wrong result_id; match by
            # retrieval_rank so every retrieved paper still receives a class.
            ent_by_rank = {int(p.get("retrieval_rank", -1)): p
                           for p in parsed.get("paper_entities", [])}
        else:
            schema_errs = ["no entity output"]
            ent_by_id, ent_by_rank = {}, {}

        for _, r in sub.iterrows():
            e = ent_by_id.get(r["result_id"]) or ent_by_rank.get(int(r["retrieval_rank"]), {})
            match = e.get("entity_match", "unresolved")
            eligible = match in ELIGIBLE
            if match == "wrong_entity":
                stats["collisions"].append((tf, r["result_id"]))
            rows.append({
                "TF": tf, "paperclip_result_id": r["result_id"],
                "retrieval_rank": int(r["retrieval_rank"]),
                "candidate_symbol": tf, "candidate_full_name": full_name,
                "recognised_gene_aliases": aliases,
                "entity_mentioned": e.get("entity_mentioned", ""),
                "entity_match": match,
                "entity_match_rationale": e.get("entity_match_rationale", ""),
                "entity_validation_model": args.model,
                "entity_validation_prompt_sha256": prompt_sha,
                "entity_validation_confidence": e.get("entity_validation_confidence", ""),
                "entity_eligible": eligible,
                "schema_valid": not schema_errs,
                "human_entity_check": "not_reviewed",
                "human_entity_check_notes": "",
            })

    ev = pd.DataFrame(rows)
    ev.to_csv(out / "paperclip_entity_validation.csv", index=False)
    manifest = {
        "model": args.model, "prompt_sha256": prompt_sha,
        "schema_sha256": base.sha256_file(repo / SCHEMA_REL),
        "hgnc_provenance": json.loads((out / "hgnc_provenance.json").read_text())
        if (out / "hgnc_provenance.json").exists() else {},
        "candidates": len(tfs), "rows": len(ev),
        "attempted": stats["attempted"], "completed": stats["completed"],
        "retries": stats["retries"], "parse_failures": stats["parse_failures"],
        "wrong_entity_collisions": stats["collisions"],
        "completion_timestamp_utc": base.utcnow(),
    }
    (out / "entity_run_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\n[entity] rows={len(ev)} completed={stats['completed']}/{len(tfs)} "
          f"collisions(wrong_entity)={stats['collisions']} parse_failures={stats['parse_failures']}")
    if len(ev):
        print("[entity] entity_match distribution:")
        print(ev["entity_match"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
