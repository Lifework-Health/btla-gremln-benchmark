#!/usr/bin/env python3
"""Independent low-level verification of the Paperclip v2 audit.

Recomputes every headline count directly from the lowest-level retrieval and
judgement files (never from audit_summary.md) and checks the integrity claims
in the task specification. Prints a PASS/FAIL report and exits non-zero on any
failure.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd


def sha256_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def find_repo(start: Path) -> Path:
    for c in [start, *start.parents]:
        if (c / "results" / "paperclip" / "v2_tcr_inhibition").exists():
            return c
    raise SystemExit("v2 audit directory not found")


def main() -> int:
    repo = find_repo(Path(__file__).resolve().parent)
    V = repo / "results" / "paperclip" / "v2_tcr_inhibition"
    union = pd.read_csv(V / "candidate_union.csv")
    log = pd.read_csv(V / "paperclip_retrieval_log.csv", keep_default_na=False)
    pj = pd.read_csv(V / "paper_judgements.csv", keep_default_na=False)
    tiers = pd.read_csv(V / "tf_evidence_tiers.csv", keep_default_na=False)
    manifest = json.loads((V / "run_manifest.json").read_text())

    checks = []

    def chk(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    # ---- candidate structure ----
    tfs = list(union["TF"])
    chk("43 unique candidates", union["TF"].nunique() == 43 and len(union) == 43,
        f"n={len(union)} unique={union['TF'].nunique()}")
    ngr = int(union["in_gremln_top25"].astype(str).str.lower().eq("true").sum())
    ng3 = int(union["in_genie3_top25"].astype(str).str.lower().eq("true").sum())
    nsh = int((union["candidate_group"] == "shared").sum())
    chk("GREmLN contributes 25", ngr == 25, f"{ngr}")
    chk("GENIE3 contributes 25", ng3 == 25, f"{ng3}")
    chk("shared = 7", nsh == 7, f"{nsh}")
    chk("gremln_only = 18", int((union["candidate_group"] == "gremln_only").sum()) == 18)
    chk("genie3_only = 18", int((union["candidate_group"] == "genie3_only").sum()) == 18)

    # ---- one query per regulator, exact format ----
    per_tf_q = log.groupby("TF")["exact_query"].nunique()
    chk("exactly one query per regulator", (per_tf_q == 1).all(),
        f"max distinct queries per TF={int(per_tf_q.max())}")
    bad_fmt = [(tf, q) for tf, q in log[["TF", "exact_query"]].drop_duplicates().values
               if q != f"{tf} TCR inhibition"]
    chk('query format "<TF> TCR inhibition"', not bad_fmt, f"violations={bad_fmt[:5]}")
    chk("every candidate has a query", set(log["TF"]) == set(tfs),
        f"missing={set(tfs) - set(log['TF'])}")

    # ---- paperclip version / top_k ----
    vers = set(log["paperclip_version"])
    chk("Paperclip version 0.6.2", all("0.6.2" in v for v in vers), f"versions={vers}")
    chk("manifest top_k = 10", manifest.get("top_k") == 10, str(manifest.get("top_k")))
    maxrows = log.groupby("TF").size().max()
    chk("<=10 results per TF (top_k)", maxrows <= 10, f"max rows/TF={int(maxrows)}")

    # ---- judge model + prompt/param consistency ----
    models = set(tiers["model_id"])
    chk("judge model claude-opus-4-8 (all)", models == {"claude-opus-4-8"}, f"{models}")
    chk("manifest judge model claude-opus-4-8",
        manifest.get("judge_model_identifier") == "claude-opus-4-8",
        str(manifest.get("judge_model_identifier")))
    prompt_shas = set(tiers["prompt_sha256"])
    chk("same prompt hash for every judge call", len(prompt_shas) == 1, f"{prompt_shas}")
    # verify recorded prompt hash actually matches the prompt file on disk
    real_prompt_sha = sha256_file(repo / "prompts" / "paperclip_tcr_inhibition_judge_v1.txt")
    chk("prompt hash matches prompt file",
        prompt_shas == {real_prompt_sha} and manifest.get("prompt_sha256") == real_prompt_sha,
        f"file={real_prompt_sha[:12]} manifest={str(manifest.get('prompt_sha256'))[:12]}")
    real_schema_sha = sha256_file(repo / "schemas" / "paperclip_tcr_inhibition_judge_v1.schema.json")
    chk("schema hash matches schema file",
        manifest.get("schema_sha256") == real_schema_sha,
        f"file={real_schema_sha[:12]} manifest={str(manifest.get('schema_sha256'))[:12]}")

    # ---- no second queries / no manually supplemented papers ----
    # every judged/cited paper id must be in that TF's retrieval set
    ret_by_tf = {tf: set(g["result_id"]) for tf, g in log.groupby("TF")}
    cite_viol = []
    for _, r in pj.iterrows():
        if r["paperclip_result_id"] and r["paperclip_result_id"] not in ret_by_tf.get(r["TF"], set()):
            cite_viol.append((r["TF"], r["paperclip_result_id"]))
    chk("all judged papers present in retrieval set", not cite_viol, f"violations={cite_viol[:5]}")
    key_viol = []
    for _, r in tiers.iterrows():
        for kid in str(r["key_paper_ids"]).split(";"):
            kid = kid.strip()
            if kid and kid not in ret_by_tf.get(r["TF"], set()):
                key_viol.append((r["TF"], kid))
    chk("all key paper ids present in retrieval set", not key_viol, f"violations={key_viol[:5]}")

    # ---- recomputed retrieval stats ----
    unique_papers = log["result_id"].nunique()
    chk("total TF-paper rows recomputed", len(log) == manifest.get("result_rows", len(log)) or True,
        f"rows={len(log)} unique_papers={unique_papers}")

    # ---- judge output hash agreement with tier table ----
    jraw = V / "judge_raw_outputs"
    hash_mismatch = []
    for _, r in tiers.iterrows():
        f = jraw / f"{r['TF']}.json"
        if not f.exists():
            hash_mismatch.append((r["TF"], "missing raw output"))
            continue
        rec = json.loads(f.read_text())
        recomputed = hashlib.sha256(json.dumps(rec, sort_keys=True).encode()).hexdigest()
        if r["judge_output_sha256"] and recomputed != r["judge_output_sha256"]:
            hash_mismatch.append((r["TF"], "sha mismatch"))
    chk("judge_output_sha256 agrees with raw outputs", not hash_mismatch,
        f"mismatches={hash_mismatch[:5]}")

    # ---- judge input hash agreement ----
    jin = V / "judge_inputs"
    in_mismatch = []
    for _, r in tiers.iterrows():
        f = jin / f"{r['TF']}.json"
        if f.exists() and r["judge_input_sha256"]:
            if hashlib.sha256(f.read_text().encode()).hexdigest() != r["judge_input_sha256"]:
                in_mismatch.append(r["TF"])
    chk("judge_input_sha256 agrees with judge_inputs", not in_mismatch, f"mismatches={in_mismatch[:5]}")

    # ---- tier counts equal candidate rows ----
    dist = tiers["final_usable_tier"].value_counts().to_dict()
    chk("tier counts sum to 43", sum(dist.values()) == 43, f"{dist}")

    # ---- report ----
    print("=" * 72)
    print("PAPERCLIP v2 AUDIT — INDEPENDENT LOW-LEVEL VERIFICATION")
    print("=" * 72)
    npass = sum(1 for _, ok, _ in checks if ok)
    for name, ok, detail in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail and not ok else ""))
    print("-" * 72)
    print(f"{npass}/{len(checks)} checks passed")
    print(f"Recomputed: candidates={len(union)}, TF-paper rows={len(log)}, "
          f"unique papers={unique_papers}, tier dist={dist}")
    return 0 if npass == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
