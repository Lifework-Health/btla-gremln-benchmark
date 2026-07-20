#!/usr/bin/env python3
"""Auditable Paperclip TCR-inhibition literature-evidence pipeline (v2).

Single fixed query per candidate regulator:

    <TF> TCR inhibition

Stages
------
freeze     Freeze the candidate union from the canonical seed-excluded rankings.
retrieve   Run one Paperclip query per TF, top-10, with a full audit trail.
inputs     Build one judge-input evidence packet per TF (+ SHA256).
judge      Run one fixed LLM judge call per TF via the Cursor SDK.
rubric     Deterministic post-judge rubric validation.
summary    Emit tf_evidence_tiers.csv, paper_judgements.csv, run_manifest.json,
           audit_summary.md and the legacy-vs-audited diagnostic.
all        freeze -> retrieve -> inputs -> judge -> rubric -> summary.

Design notes
------------
* The Paperclip CLI is the only retrieval backend. We never simulate it with a
  web search and we fail explicitly if it is unavailable.
* Raw Paperclip responses (which may contain copyrighted abstracts/snippets)
  are written to a git-ignored ``raw/`` directory; only a manifest with
  identifiers, byte sizes and SHA256 hashes is committed.
* The judge is a single fixed Cursor SDK model. Cursor agent models do not
  expose temperature/top_p/seed, so those are recorded as ``not_supported``.
  Determinism across reruns is provided by caching raw judge outputs.
* Reruns reuse cached Paperclip responses and cached judge outputs, so the
  pipeline is rerunnable without changing results.

This module is import-safe: the pure helpers (build_query, truncate_top_k,
run_rubric_checks, validate_judge_json, extract_json) have no side effects and
are covered by tests/test_paperclip_tcr_inhibition_audit.py.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# --------------------------------------------------------------------------- #
# Constants (fixed by the task specification)
# --------------------------------------------------------------------------- #
SCRIPT_VERSION = "v1.0.0"
QUERY_SUFFIX = "TCR inhibition"
QUERY_TEMPLATE = "{TF} " + QUERY_SUFFIX
TOP_K = 10
PAPERCLIP_SOURCE = "pmc"  # source scope only; NOT a query term
PROMPT_REL = "prompts/paperclip_tcr_inhibition_judge_v1.txt"
SCHEMA_REL = "schemas/paperclip_tcr_inhibition_judge_v1.schema.json"
OUT_REL = "results/paperclip/v2_tcr_inhibition"
DRY_TFS = ["EGR2", "AHR", "KIF22"]
DEFAULT_JUDGE_MODEL = os.environ.get("PAPERCLIP_JUDGE_MODEL", "claude-4.5-sonnet")
JUDGE_MODEL_FALLBACK = "composer-2.5"


# --------------------------------------------------------------------------- #
# Small pure helpers (unit-tested)
# --------------------------------------------------------------------------- #
def build_query(tf: str) -> str:
    """Construct the single fixed Paperclip query for a candidate regulator.

    The symbol is used exactly as supplied; no aliases, no extra terms.
    """
    return QUERY_TEMPLATE.format(TF=tf)


def truncate_top_k(results: list, k: int = TOP_K) -> list:
    """Return at most ``k`` results, preserving Paperclip ranking order."""
    return list(results)[:k]


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(Path(path).read_bytes())


def utcnow() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def extract_json(text: str):
    """Extract the first well-formed top-level JSON object from LLM output.

    Handles ```json fenced blocks and leading/trailing prose. Returns a dict on
    success, or None if no parseable object is found.
    """
    if text is None:
        return None
    # Prefer a fenced ```json block.
    fence = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", text, re.DOTALL)
    candidates = []
    if fence:
        candidates.append(fence.group(1))
    # Fall back to the widest brace span.
    first, last = text.find("{"), text.rfind("}")
    if first != -1 and last != -1 and last > first:
        candidates.append(text[first:last + 1])
    for cand in candidates:
        try:
            obj = json.loads(cand)
            if isinstance(obj, dict):
                return obj
        except json.JSONDecodeError:
            continue
    return None


# --------------------------------------------------------------------------- #
# JSON-schema validation (draft-07 subset sufficient for our schema)
# --------------------------------------------------------------------------- #
def validate_judge_json(obj, schema: dict) -> list:
    """Validate ``obj`` against ``schema``. Returns a list of error strings.

    Uses jsonschema if installed; otherwise a compact built-in validator that
    covers the constructs used in paperclip_tcr_inhibition_judge_v1.schema.json
    (type, required, enum, additionalProperties, array items, minimum/minLength).
    """
    try:
        import jsonschema  # type: ignore

        v = jsonschema.Draft7Validator(schema)
        return [f"{'/'.join(map(str, e.path))}: {e.message}" for e in v.iter_errors(obj)]
    except Exception:
        return _mini_validate(obj, schema, "$")


_TYPE_MAP = {
    "object": dict, "array": list, "string": str,
    "integer": int, "number": (int, float), "boolean": bool,
}


def _mini_validate(obj, schema, path) -> list:
    errs = []
    t = schema.get("type")
    if t:
        py = _TYPE_MAP[t]
        # bool is a subclass of int; reject bool where integer/number expected.
        if t in ("integer", "number") and isinstance(obj, bool):
            return [f"{path}: expected {t}, got boolean"]
        if not isinstance(obj, py):
            return [f"{path}: expected {t}, got {type(obj).__name__}"]
    if "enum" in schema and obj not in schema["enum"]:
        errs.append(f"{path}: {obj!r} not in enum {schema['enum']}")
    if isinstance(obj, str) and "minLength" in schema and len(obj) < schema["minLength"]:
        errs.append(f"{path}: shorter than minLength {schema['minLength']}")
    if isinstance(obj, (int, float)) and not isinstance(obj, bool) and "minimum" in schema:
        if obj < schema["minimum"]:
            errs.append(f"{path}: below minimum {schema['minimum']}")
    if t == "object" and isinstance(obj, dict):
        for req in schema.get("required", []):
            if req not in obj:
                errs.append(f"{path}: missing required '{req}'")
        props = schema.get("properties", {})
        if schema.get("additionalProperties") is False:
            for k in obj:
                if k not in props:
                    errs.append(f"{path}: unexpected property '{k}'")
        for k, sub in props.items():
            if k in obj:
                errs += _mini_validate(obj[k], sub, f"{path}.{k}")
    if t == "array" and isinstance(obj, list) and "items" in schema:
        for i, item in enumerate(obj):
            errs += _mini_validate(item, schema["items"], f"{path}[{i}]")
    return errs


# --------------------------------------------------------------------------- #
# Deterministic rubric checks
# --------------------------------------------------------------------------- #
def run_rubric_checks(tf: str, judge: dict, retrieval_ids: list, exact_query: str,
                      papers_retrieved_log: int) -> tuple:
    """Return (rubric_valid: bool, reasons: list[str]).

    ``judge`` is the parsed judge object; ``retrieval_ids`` is the list of
    result IDs present in this TF's retrieval log; ``papers_retrieved_log`` is
    the row count for this TF in the retrieval log.
    """
    reasons = []
    if judge is None:
        return False, ["judge output missing or unparseable"]

    tier = str(judge.get("overall_evidence_tier", "")).lower()
    assessments = judge.get("paper_assessments", []) or []
    keys = judge.get("key_paperclip_result_ids", []) or []

    # exactly one tier
    if tier not in {"strong", "moderate", "weak", "none"}:
        reasons.append(f"invalid or missing overall tier: {tier!r}")

    # query must equal "<TF> TCR inhibition"
    if str(judge.get("query", "")) != exact_query:
        reasons.append(f"judge query {judge.get('query')!r} != exact query {exact_query!r}")

    # papers_retrieved must match the retrieval log
    if int(judge.get("papers_retrieved", -1)) != int(papers_retrieved_log):
        reasons.append(
            f"papers_retrieved={judge.get('papers_retrieved')} != retrieval log rows={papers_retrieved_log}")

    # supporting excerpts <= 20 words
    for a in assessments:
        exc = str(a.get("supporting_excerpt", "") or "")
        if len(exc.split()) > 20:
            reasons.append(
                f"supporting_excerpt >20 words for {a.get('paperclip_result_id')}")

    # only supplied result IDs may be cited (assessments + keys)
    known = set(retrieval_ids)
    for a in assessments:
        rid = a.get("paperclip_result_id")
        if rid and rid not in known:
            reasons.append(f"assessment cites unknown result id {rid!r}")
    for rid in keys:
        if rid not in known:
            reasons.append(f"key paper id {rid!r} not in retrieval log")

    def _relevant_primary():
        return [a for a in assessments
                if a.get("relevant") and a.get("study_type") == "primary experimental"]

    def _direct_causal_primary():
        return [a for a in _relevant_primary()
                if a.get("evidence_directness") == "direct causal or perturbational"]

    def _has_nonreview_primary():
        return len(_relevant_primary()) > 0

    # tier-specific structural requirements
    if tier in {"strong", "moderate"} and not keys:
        reasons.append(f"{tier} tier has no key_paperclip_result_ids")
    if tier == "strong" and not _direct_causal_primary():
        reasons.append("strong tier lacks a relevant direct-causal primary study")
    if tier == "moderate" and not _has_nonreview_primary():
        reasons.append("moderate tier lacks a relevant primary study")

    # review-only evidence cannot exceed weak
    relevant = [a for a in assessments if a.get("relevant")]
    if relevant and all(a.get("study_type") == "review" for a in relevant):
        if tier in {"strong", "moderate"}:
            reasons.append("review-only evidence cannot exceed weak")

    # zero retrieved papers must be none
    if papers_retrieved_log == 0 and tier != "none":
        reasons.append("zero retrieved papers but tier is not none")

    return (len(reasons) == 0), reasons


# --------------------------------------------------------------------------- #
# Repo / path helpers
# --------------------------------------------------------------------------- #
def find_repo_root(start: Path) -> Path:
    p = start.resolve()
    for cand in [p, *p.parents]:
        if (cand / "scripts").is_dir() and (cand / "results").is_dir() \
                and (cand / ".git").exists():
            return cand
    return p


def git_commit(repo: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


# --------------------------------------------------------------------------- #
# Paperclip CLI wrappers
# --------------------------------------------------------------------------- #
def paperclip_bin() -> str:
    b = shutil.which("paperclip") or os.path.expanduser("~/.local/bin/paperclip")
    if not Path(b).exists():
        raise RuntimeError("Paperclip CLI not found. This pipeline requires the "
                           "real Paperclip tool and will not simulate it.")
    return b


def paperclip_version(bin_: str) -> str:
    try:
        out = subprocess.run([bin_, "--version"], capture_output=True, text=True,
                             timeout=60).stdout.strip()
        return out or "unknown"
    except Exception:
        return "unknown"


def _run_paperclip(bin_: str, args: list, timeout: int = 180) -> subprocess.CompletedProcess:
    return subprocess.run([bin_, *args], capture_output=True, text=True, timeout=timeout)


def paperclip_search(bin_: str, query: str, n: int) -> tuple:
    """Run one search. Returns (result_id, raw_stdout, command_str)."""
    args = ["search", "-s", PAPERCLIP_SOURCE, query, "-n", str(n)]
    cmd = f"{Path(bin_).name} " + " ".join(
        (f'"{a}"' if " " in a else a) for a in args)
    cp = _run_paperclip(bin_, args)
    raw = (cp.stdout or "") + (("\n[stderr]\n" + cp.stderr) if cp.stderr else "")
    m = re.search(r"\[(s_[0-9a-fA-F]+)\]", cp.stdout or "")
    return (m.group(1) if m else None), raw, cmd


def paperclip_export(bin_: str, result_id: str, dest: Path) -> list:
    """Export a saved search result set to CSV and return ordered rows."""
    _run_paperclip(bin_, ["results", result_id, "--save", str(dest)])
    if not dest.exists():
        return []
    with dest.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def paperclip_meta(bin_: str, result_id: str) -> tuple:
    """Fetch /papers/<id>/meta.json. Returns (meta_dict_or_None, raw_text)."""
    cp = _run_paperclip(bin_, ["cat", f"/papers/{result_id}/meta.json"])
    raw = cp.stdout or ""
    # strip a trailing "[NNms]" timing line if present
    body = re.sub(r"\n\[\d+ms\]\s*$", "", raw).strip()
    try:
        return json.loads(body), raw
    except json.JSONDecodeError:
        return None, raw


# --------------------------------------------------------------------------- #
# Stage: freeze candidate union
# --------------------------------------------------------------------------- #
def stage_freeze(repo: Path, out: Path) -> dict:
    import pandas as pd

    tbl = repo / "results" / "tables"
    gr_path = tbl / "gremln_btla_vs_tcr_seed_excluded_tf_ranking.csv"
    g3_path = tbl / "genie3_btla_vs_tcr_seed_excluded_tf_ranking.csv"
    for p in (gr_path, g3_path):
        if not p.exists():
            raise FileNotFoundError(
                f"Canonical ranking file missing: {p}. Regenerate via nb03 first.")

    gr = pd.read_csv(gr_path)
    g3 = pd.read_csv(g3_path)
    gr_rank = dict(zip(gr["gene"], gr["gremln_csls_rank"]))
    g3_rank = dict(zip(g3["gene"], g3["genie3_rank"]))

    # GENIE3: unambiguous continuous score -> top 25 by rank.
    g3_derived = list(g3.sort_values("genie3_rank").reset_index(drop=True).head(25)["gene"])
    # GREmLN: integer seed-neighbour count with heavy ties at the top-25 boundary.
    # A naive head(25) is order-sensitive within a tied dense rank, so it is NOT the
    # authoritative membership.
    gr_derived = list(gr.head(25)["gene"])

    # Authoritative membership = the canonical materialised seed-excluded top-25
    # (top25_union_primary.csv), whose GREmLN boundary tie is resolved by the
    # summed-CSLS tie-break already applied in nb03. This guarantees the frozen
    # candidate set is identical to the published Table 3 / Table 4 candidates.
    canon = repo / "results" / "publication_data" / "top25_union_primary.csv"
    membership_source = "canonical top25_union_primary.csv"
    if canon.exists():
        cu = pd.read_csv(canon)
        gr_top = list(cu[cu["in_gremln_top25"] == True]["TF"])
        g3_top = list(cu[cu["in_genie3_top25"] == True]["TF"])
        # GENIE3 is unambiguous: derived head(25) must equal canonical.
        assert set(g3_top) == set(g3_derived), (
            "GENIE3 canonical top-25 differs from derived head(25): "
            f"canon-only={set(g3_top)-set(g3_derived)}, derived-only={set(g3_derived)-set(g3_top)}")
        gr_boundary = set(gr_derived) ^ set(gr_top)
        if gr_boundary:
            print(f"[freeze] GREmLN boundary tie resolved via canonical SumCSLS "
                  f"tie-break (differs from naive head(25) at: {sorted(gr_boundary)})")
    else:
        membership_source = "derived head(25) from ranking CSVs (canonical file absent)"
        gr_top, g3_top = gr_derived, g3_derived
        print("[freeze] WARNING: top25_union_primary.csv absent; using derived head(25).")

    union = sorted(set(gr_top) | set(g3_top))
    rows = []
    commit = git_commit(repo)
    gr_sha, g3_sha = sha256_file(gr_path), sha256_file(g3_path)
    combined_sha = sha256_text(gr_sha + g3_sha)
    for tf in union:
        in_gr = tf in gr_top
        in_g3 = tf in g3_top
        group = "shared" if (in_gr and in_g3) else ("gremln_only" if in_gr else "genie3_only")
        rows.append({
            "TF": tf,
            "gremln_rank": int(gr_rank[tf]) if tf in gr_rank else "",
            "genie3_rank": int(g3_rank[tf]) if tf in g3_rank else "",
            "in_gremln_top25": in_gr,
            "in_genie3_top25": in_g3,
            "candidate_group": group,
            "ranking_source_commit": commit,
            "ranking_file_sha256": combined_sha,
            "gremln_ranking_sha256": gr_sha,
            "genie3_ranking_sha256": g3_sha,
        })
    df = pd.DataFrame(rows)

    # ---- asserts ----
    assert int(df["in_gremln_top25"].sum()) == 25, \
        f"GREmLN must contribute 25, got {int(df['in_gremln_top25'].sum())}"
    assert int(df["in_genie3_top25"].sum()) == 25, \
        f"GENIE3 must contribute 25, got {int(df['in_genie3_top25'].sum())}"
    assert df["TF"].is_unique, "duplicate TF rows in candidate union"
    assert (df["in_gremln_top25"] | df["in_genie3_top25"]).all(), \
        "every union member must belong to at least one top-25 list"

    df.to_csv(out / "candidate_union.csv", index=False)
    print(f"[freeze] candidate_union.csv: {len(df)} TFs "
          f"(GREmLN 25, GENIE3 25, shared {int((df.candidate_group=='shared').sum())})")
    return {
        "candidate_count": len(df),
        "ranking_source_commit": commit,
        "gremln_ranking_sha256": gr_sha,
        "genie3_ranking_sha256": g3_sha,
        "ranking_file_sha256": combined_sha,
        "membership_source": membership_source,
        "tfs": list(df["TF"]),
    }


# --------------------------------------------------------------------------- #
# Stage: retrieval
# --------------------------------------------------------------------------- #
def stage_retrieve(repo: Path, out: Path, tfs: list, bin_: str, version: str,
                   use_cache: bool) -> dict:
    import pandas as pd

    raw_dir = out / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    per_tf_dir = raw_dir / "per_tf"
    per_tf_dir.mkdir(exist_ok=True)

    log_rows = []
    manifest_rows = []
    jsonl_path = raw_dir / "paperclip_raw_responses.jsonl"
    jsonl_fh = jsonl_path.open("w", encoding="utf-8")
    meta_cache: dict = {}
    queries_completed = 0
    tfs_with_zero = []
    tfs_under_10 = []

    for tf in tfs:
        query = build_query(tf)
        cache_file = per_tf_dir / f"{tf}.json"
        if use_cache and cache_file.exists():
            cached = json.loads(cache_file.read_text())
            result_id = cached["result_id"]
            raw_search = cached["raw_search"]
            command = cached["command"]
            rows = cached["results"]
            ts = cached["retrieval_timestamp_utc"]
            print(f"[retrieve] {tf}: using cached response ({len(rows)} results)")
        else:
            ts = utcnow()
            result_id, raw_search, command = paperclip_search(bin_, query, TOP_K)
            rows = []
            if result_id:
                tmp = raw_dir / f"_export_{tf}.csv"
                exported = truncate_top_k(paperclip_export(bin_, result_id, tmp), TOP_K)
                if tmp.exists():
                    tmp.unlink()
                for rank, r in enumerate(exported, start=1):
                    rid = (r.get("id") or "").strip()
                    if rid in meta_cache:
                        meta, meta_raw = meta_cache[rid]
                    else:
                        meta, meta_raw = paperclip_meta(bin_, rid)
                        meta_cache[rid] = (meta, meta_raw)
                    meta = meta or {}
                    rows.append({
                        "retrieval_rank": rank,
                        "result_id": rid,
                        "title": r.get("title", ""),
                        "authors": r.get("authors", ""),
                        "source": r.get("source", ""),
                        "date": r.get("date", ""),
                        "url": r.get("url", ""),
                        "snippet": r.get("abstract", ""),
                        "doi": meta.get("doi", ""),
                        "pmid": meta.get("pmid", ""),
                        "pmcid": meta.get("pmc_id", ""),
                        "journal": meta.get("journal", ""),
                        "publication_type": meta.get("article_type", ""),
                        "year": meta.get("pub_year", ""),
                        "abstract": meta.get("abstract", ""),
                        "meta_raw": meta_raw,
                    })
            cache_file.write_text(json.dumps({
                "tf": tf, "query": query, "result_id": result_id,
                "command": command, "raw_search": raw_search,
                "retrieval_timestamp_utc": ts, "results": rows,
            }, indent=2))
            print(f"[retrieve] {tf}: query={query!r} id={result_id} results={len(rows)}")
            command_ = command

        if len(rows) == 0:
            tfs_with_zero.append(tf)
        elif len(rows) < TOP_K:
            tfs_under_10.append(tf)
        queries_completed += 1

        # per-query raw search hash into manifest
        search_sha = sha256_text(raw_search)
        manifest_rows.append({
            "tf": tf, "query": query, "artifact": "search_stdout",
            "result_id": result_id or "", "bytes": len(raw_search.encode("utf-8")),
            "sha256": search_sha,
        })

        for r in rows:
            other_id = ""
            for cand in (r.get("pmid"), r.get("result_id")):
                if cand and cand != r.get("pmcid"):
                    other_id = cand
                    break
            meta_raw = r.get("meta_raw", "")
            raw_sha = sha256_text(meta_raw) if meta_raw else ""
            full_text = str(r.get("source", "")).lower() in {"pmc", "biorxiv", "medrxiv",
                                                             "arxiv", "biomedrxiv"}
            log_rows.append({
                "TF": tf,
                "exact_query": query,
                "retrieval_rank": r["retrieval_rank"],
                "result_id": r["result_id"],
                "relevance_score": "",  # not provided by the Paperclip CLI
                "title": r["title"],
                "authors": r["authors"],
                "year": r["year"],
                "journal": r["journal"],
                "doi": r["doi"],
                "pmcid": r["pmcid"],
                "other_source_id": other_id,
                "source_url": r["url"],
                "publication_type": r["publication_type"],
                "abstract_available": bool(str(r.get("abstract", "")).strip()),
                "full_text_available": full_text,
                "full_text_available_method": "source_heuristic(/papers=full-text)",
                "abstract": r["abstract"],
                "snippet": r["snippet"],
                "retrieval_timestamp_utc": ts,
                "raw_response_sha256": raw_sha,
                "paperclip_version": version,
                "paperclip_command": command if not (use_cache and cache_file.exists()) else command,
            })
            manifest_rows.append({
                "tf": tf, "query": query, "artifact": f"meta:{r['result_id']}",
                "result_id": r["result_id"],
                "bytes": len(str(r.get("meta_raw", "")).encode("utf-8")),
                "sha256": raw_sha,
            })
            jsonl_fh.write(json.dumps({
                "tf": tf, "query": query, "retrieval_rank": r["retrieval_rank"],
                "result_id": r["result_id"], "meta_raw": r.get("meta_raw", ""),
                "snippet": r.get("snippet", ""),
            }) + "\n")

    jsonl_fh.close()

    log_cols = ["TF", "exact_query", "retrieval_rank", "result_id", "relevance_score",
                "title", "authors", "year", "journal", "doi", "pmcid", "other_source_id",
                "source_url", "publication_type", "abstract_available", "full_text_available",
                "full_text_available_method", "abstract", "snippet",
                "retrieval_timestamp_utc", "raw_response_sha256", "paperclip_version",
                "paperclip_command"]
    log_df = pd.DataFrame(log_rows, columns=log_cols)
    log_df.to_csv(out / "paperclip_retrieval_log.csv", index=False)
    pd.DataFrame(manifest_rows).to_csv(out / "paperclip_raw_manifest.csv", index=False)

    unique_papers = log_df["result_id"].nunique() if len(log_df) else 0
    print(f"[retrieve] {len(tfs)} TFs, {len(log_df)} result rows, "
          f"{unique_papers} unique papers; zero-result TFs: {tfs_with_zero}")
    return {
        "queries_attempted": len(tfs),
        "queries_completed": queries_completed,
        "result_rows": len(log_df),
        "unique_papers": int(unique_papers),
        "tfs_zero_results": tfs_with_zero,
        "tfs_under_10": tfs_under_10,
    }


# --------------------------------------------------------------------------- #
# Stage: judge inputs
# --------------------------------------------------------------------------- #
def stage_inputs(out: Path, tfs: list) -> dict:
    import pandas as pd

    log = pd.read_csv(out / "paperclip_retrieval_log.csv", keep_default_na=False)
    jin_dir = out / "judge_inputs"
    jin_dir.mkdir(exist_ok=True)
    shas = {}
    for tf in tfs:
        sub = log[log["TF"] == tf].sort_values("retrieval_rank")
        papers = []
        for _, r in sub.iterrows():
            papers.append({
                "paperclip_result_id": r["result_id"],
                "rank": int(r["retrieval_rank"]),
                "title": r["title"],
                "year": r["year"],
                "publication_type": r["publication_type"],
                "abstract": r["abstract"],
                "relevant_passages": r["snippet"],
                "doi_or_pmcid": r["doi"] or r["pmcid"],
            })
        packet = {
            "tf": tf,
            "query": build_query(tf),
            "top_k": TOP_K,
            "papers_retrieved": len(papers),
            "papers": papers,
        }
        text = json.dumps(packet, indent=2, ensure_ascii=False)
        (jin_dir / f"{tf}.json").write_text(text, encoding="utf-8")
        shas[tf] = sha256_text(text)
    print(f"[inputs] built {len(tfs)} judge-input packets")
    return {"judge_input_sha256": shas}


# --------------------------------------------------------------------------- #
# Stage: judge (Cursor SDK)
# --------------------------------------------------------------------------- #
def _render_prompt(template: str, tf: str, n: int) -> str:
    # Only {TF} and {N} are placeholders; the JSON example braces stay literal.
    return template.replace("{TF}", tf).replace("{N}", str(n))


def _call_cursor_judge(prompt: str, model: str, api_key: str):
    """Return (text, run_id, agent_id, resolved_model, capture_method).

    The Cursor SDK judge is an *agentic* model with file-write tools. To keep
    the audit clean we run it in an isolated throwaway working directory (never
    the repo), so it cannot pollute the repository or read benchmark files. If
    the model routes its answer to a file tool instead of returning it as text,
    we recover that file's contents as the raw output (``capture_method`` records
    which path was used).

    Raises on startup failure or run status == error.
    """
    import tempfile
    from cursor_sdk import Agent, AgentOptions, LocalAgentOptions

    tmp = Path(tempfile.mkdtemp(prefix="paperclip_judge_"))
    capture = "returned_text"
    try:
        result = Agent.prompt(
            prompt,
            AgentOptions(api_key=api_key, model=model,
                         local=LocalAgentOptions(cwd=str(tmp), setting_sources=[])),
        )
        if getattr(result, "status", "finished") == "error":
            raise RuntimeError(f"judge run status=error id={getattr(result, 'id', '')}")
        text = getattr(result, "result", None) or getattr(result, "text", None) or ""
        # Fallback: if the agent wrote its JSON answer to a file instead of
        # returning it, recover the judge object from the isolated cwd.
        if extract_json(text) is None:
            best = None
            for jf in sorted(tmp.rglob("*.json")):
                try:
                    body = jf.read_text(encoding="utf-8")
                except Exception:
                    continue
                obj = extract_json(body)
                if isinstance(obj, dict) and "overall_evidence_tier" in obj:
                    best = body
            if best is not None:
                text = best
                capture = "agent_file_tool"
        return (text, getattr(result, "id", ""), getattr(result, "agent_id", ""),
                model, capture)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def stage_judge(repo: Path, out: Path, tfs: list, model: str, use_cache: bool) -> dict:
    import pandas as pd

    api_key = os.environ.get("CURSOR_API_KEY")
    prompt_template = (repo / PROMPT_REL).read_text(encoding="utf-8")
    schema = json.loads((repo / SCHEMA_REL).read_text(encoding="utf-8"))
    jin_dir = out / "judge_inputs"
    jraw_dir = out / "judge_raw_outputs"
    jraw_dir.mkdir(exist_ok=True)

    stats = {"attempted": 0, "completed": 0, "retries": 0, "parse_failures": [],
             "skipped_no_key": False, "model": model, "results": {}}

    if not api_key:
        stats["skipped_no_key"] = True
        print("[judge] CURSOR_API_KEY not set -> judge stage SKIPPED. "
              "Set CURSOR_API_KEY and rerun `--stage judge`.")
        return stats

    try:
        import cursor_sdk  # noqa: F401
    except Exception as e:  # pragma: no cover
        stats["skipped_no_key"] = True
        stats["import_error"] = str(e)
        print(f"[judge] cursor_sdk import failed: {e} -> judge SKIPPED.")
        return stats

    for tf in tfs:
        cache_file = jraw_dir / f"{tf}.json"
        packet = json.loads((jin_dir / f"{tf}.json").read_text())
        n = packet["papers_retrieved"]
        prompt = (_render_prompt(prompt_template, tf, n)
                  + "\n\nRETRIEVED PAPERS (JSON):\n"
                  + json.dumps(packet, ensure_ascii=False, indent=2))

        if use_cache and cache_file.exists():
            record = json.loads(cache_file.read_text())
            print(f"[judge] {tf}: using cached judge output")
        else:
            stats["attempted"] += 1
            attempts = []
            parsed = None
            for attempt in (1, 2):
                extra = "" if attempt == 1 else (
                    "\n\nYour previous response was not valid JSON. Return ONLY the "
                    "JSON object required by the schema, with no prose or code fences.")
                try:
                    text, run_id, agent_id, resolved, capture = _call_cursor_judge(
                        prompt + extra, model, api_key)
                except Exception as e:
                    attempts.append({"attempt": attempt, "startup_error": str(e)})
                    break
                obj = extract_json(text)
                attempts.append({"attempt": attempt, "run_id": run_id,
                                 "agent_id": agent_id, "raw": text,
                                 "capture_method": capture,
                                 "parsed_ok": obj is not None})
                if obj is not None:
                    parsed = obj
                    break
                stats["retries"] += 1
            record = {"tf": tf, "model": model, "timestamp_utc": utcnow(),
                      "attempts": attempts,
                      "parsed": parsed}
            cache_file.write_text(json.dumps(record, indent=2, ensure_ascii=False))
            if parsed is None:
                stats["parse_failures"].append(tf)

        parsed = record.get("parsed")
        if parsed is not None:
            stats["completed"] += 1
        stats["results"][tf] = {
            "parsed": parsed,
            "output_sha256": sha256_text(json.dumps(record, sort_keys=True)),
        }

    print(f"[judge] attempted={stats['attempted']} completed={stats['completed']} "
          f"retries={stats['retries']} parse_failures={stats['parse_failures']}")
    return stats


# --------------------------------------------------------------------------- #
# Stage: rubric + summary
# --------------------------------------------------------------------------- #
def stage_rubric_and_summary(repo: Path, out: Path, tfs: list, judge_stats: dict,
                             retr_stats: dict, freeze_stats: dict, model: str,
                             version: str, run_meta: dict) -> dict:
    import pandas as pd

    union = pd.read_csv(out / "candidate_union.csv")
    log = pd.read_csv(out / "paperclip_retrieval_log.csv", keep_default_na=False)
    schema = json.loads((repo / SCHEMA_REL).read_text(encoding="utf-8"))
    prompt_sha = sha256_file(repo / PROMPT_REL)
    schema_sha = sha256_file(repo / SCHEMA_REL)
    jin_dir = out / "judge_inputs"

    group = dict(zip(union["TF"], union["candidate_group"]))
    in_gr = dict(zip(union["TF"], union["in_gremln_top25"]))
    in_g3 = dict(zip(union["TF"], union["in_genie3_top25"]))

    tier_rows = []
    paper_rows = []
    tier_counts = {"strong": 0, "moderate": 0, "weak": 0, "none": 0, "missing": 0}
    rubric_fail = []

    for tf in tfs:
        exact_query = build_query(tf)
        sub = log[log["TF"] == tf]
        retrieval_ids = list(sub["result_id"])
        papers_retrieved_log = len(sub)
        jres = judge_stats.get("results", {}).get(tf, {})
        parsed = jres.get("parsed")

        judge_tier = str(parsed.get("overall_evidence_tier")).lower() if parsed else ""
        schema_errs = validate_judge_json(parsed, schema) if parsed else ["no judge output"]
        rubric_valid, reasons = run_rubric_checks(
            tf, parsed, retrieval_ids, exact_query, papers_retrieved_log)
        rubric_valid = rubric_valid and not schema_errs
        if schema_errs and parsed is not None:
            reasons = reasons + [f"schema: {e}" for e in schema_errs]

        final_tier = judge_tier if (parsed is not None and rubric_valid) else "missing"
        tier_counts[final_tier if final_tier in tier_counts else "missing"] += 1
        if not rubric_valid:
            rubric_fail.append(tf)

        assessments = (parsed or {}).get("paper_assessments", []) or []
        n_primary = sum(1 for a in assessments
                        if a.get("relevant") and a.get("study_type") == "primary experimental")
        n_direct = sum(1 for a in assessments
                       if a.get("relevant")
                       and a.get("study_type") == "primary experimental"
                       and a.get("evidence_directness") == "direct causal or perturbational")

        jin_sha = sha256_file(jin_dir / f"{tf}.json") if (jin_dir / f"{tf}.json").exists() else ""
        tier_rows.append({
            "TF": tf,
            "candidate_group": group.get(tf, ""),
            "in_gremln_top25": in_gr.get(tf, ""),
            "in_genie3_top25": in_g3.get(tf, ""),
            "exact_query": exact_query,
            "papers_retrieved": papers_retrieved_log,
            "papers_relevant": (parsed or {}).get("papers_relevant", "") if parsed else "",
            "primary_relevant_papers": n_primary if parsed else "",
            "direct_causal_papers": n_direct if parsed else "",
            "judge_tier": judge_tier,
            "rubric_valid": rubric_valid,
            "rubric_reasons": "; ".join(reasons),
            "final_usable_tier": final_tier,
            "key_paper_ids": ";".join((parsed or {}).get("key_paperclip_result_ids", []) or []),
            "tier_rationale": (parsed or {}).get("tier_rationale", "") if parsed else "",
            "confidence": (parsed or {}).get("confidence", "") if parsed else "",
            "judge_input_sha256": jin_sha,
            "judge_output_sha256": jres.get("output_sha256", ""),
            "model_id": model,
            "prompt_sha256": prompt_sha,
        })

        for a in assessments:
            paper_rows.append({
                "TF": tf,
                "paperclip_result_id": a.get("paperclip_result_id", ""),
                "retrieval_rank": a.get("retrieval_rank", ""),
                "relevant": a.get("relevant", ""),
                "study_type": a.get("study_type", ""),
                "experimental_system": a.get("experimental_system", ""),
                "evidence_directness": a.get("evidence_directness", ""),
                "phenotype": a.get("phenotype", ""),
                "direction": a.get("direction", ""),
                "supporting_excerpt": a.get("supporting_excerpt", ""),
                "rationale": a.get("rationale", ""),
            })

    tier_df = pd.DataFrame(tier_rows)
    tier_df.to_csv(out / "tf_evidence_tiers.csv", index=False)
    pd.DataFrame(paper_rows).to_csv(out / "paper_judgements.csv", index=False)

    # ---- legacy vs audited diagnostic (does NOT feed the manuscript) ----
    _write_legacy_diagnostic(repo, out, tier_df)

    # ---- run manifest ----
    manifest = {
        "repository_commit": freeze_stats.get("ranking_source_commit", git_commit(repo)),
        "ranking_input_hashes": {
            "gremln_seed_excluded": freeze_stats.get("gremln_ranking_sha256"),
            "genie3_seed_excluded": freeze_stats.get("genie3_ranking_sha256"),
            "combined": freeze_stats.get("ranking_file_sha256"),
        },
        "candidate_count": freeze_stats.get("candidate_count"),
        "exact_query_template": QUERY_TEMPLATE,
        "top_k": TOP_K,
        "paperclip_source_scope": PAPERCLIP_SOURCE,
        "paperclip_invocation": run_meta.get("paperclip_command_example", ""),
        "paperclip_version": version,
        "judge_model_identifier": model,
        "judge_parameters": {
            "temperature": "not_supported(cursor_sdk_agent)",
            "top_p": "not_supported(cursor_sdk_agent)",
            "seed": "not_supported(cursor_sdk_agent)",
            "max_output_tokens": "model_default",
            "runtime": "cursor_sdk local Agent.prompt",
            "determinism_note": "reruns reuse cached judge_raw_outputs",
        },
        "prompt_sha256": prompt_sha,
        "schema_sha256": schema_sha,
        "start_timestamp_utc": run_meta.get("start_ts"),
        "completion_timestamp_utc": utcnow(),
        "queries_attempted": retr_stats.get("queries_attempted"),
        "queries_completed": retr_stats.get("queries_completed"),
        "judge_calls_attempted": judge_stats.get("attempted"),
        "judge_calls_completed": judge_stats.get("completed"),
        "retries": judge_stats.get("retries"),
        "judge_skipped_no_key": judge_stats.get("skipped_no_key", False),
        "failures": {
            "judge_parse_failures": judge_stats.get("parse_failures", []),
            "rubric_validation_failures": rubric_fail,
            "tfs_zero_results": retr_stats.get("tfs_zero_results", []),
        },
        "script_version": SCRIPT_VERSION,
        "script_sha256": sha256_file(Path(__file__)) if Path(__file__).exists() else "",
    }
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2))

    # ---- audit summary ----
    _write_audit_summary(out, tier_df, log, retr_stats, judge_stats, tier_counts,
                         rubric_fail, model, version)
    print(f"[summary] tiers: {tier_counts}; rubric failures: {len(rubric_fail)}")
    return {"tier_counts": tier_counts, "rubric_failures": rubric_fail}


def _write_legacy_diagnostic(repo: Path, out: Path, tier_df) -> None:
    import pandas as pd

    legacy_path = repo / "results" / "paperclip" / "paperclip_union_top25_review.csv"
    if not legacy_path.exists():
        return
    legacy = pd.read_csv(legacy_path, keep_default_na=False)
    lmap = dict(zip(legacy["TF"], legacy["paperclip_evidence_tier"]))
    rows = []
    for _, r in tier_df.iterrows():
        tf = r["TF"]
        rows.append({
            "TF": tf,
            "candidate_group": r["candidate_group"],
            "legacy_tier": lmap.get(tf, "not_in_legacy"),
            "audited_judge_tier": r["judge_tier"],
            "audited_final_usable_tier": r["final_usable_tier"],
            "rubric_valid": r["rubric_valid"],
            "changed": str(lmap.get(tf, "")).lower() != str(r["final_usable_tier"]).lower(),
        })
    pd.DataFrame(rows).to_csv(out / "legacy_vs_audited_tiers.csv", index=False)


def _write_audit_summary(out, tier_df, log, retr_stats, judge_stats, tier_counts,
                         rubric_fail, model, version) -> None:
    import pandas as pd

    def _sm(mask):
        sub = tier_df[mask]
        return int((sub["final_usable_tier"].isin(["strong", "moderate"])).sum())

    n_tf = len(tier_df)
    unique_papers = log["result_id"].nunique() if len(log) else 0
    ft = int(log["full_text_available"].astype(str).str.lower().eq("true").sum()) if len(log) else 0
    ab_only = int(((log["abstract_available"].astype(str).str.lower() == "true") &
                   (log["full_text_available"].astype(str).str.lower() != "true")).sum()) if len(log) else 0
    neither = int(((log["abstract_available"].astype(str).str.lower() != "true") &
                   (log["full_text_available"].astype(str).str.lower() != "true")).sum()) if len(log) else 0

    gr = tier_df["in_gremln_top25"].astype(str).str.lower() == "true"
    g3 = tier_df["in_genie3_top25"].astype(str).str.lower() == "true"
    shared = tier_df["candidate_group"] == "shared"

    lines = [
        "# Paperclip TCR-inhibition evidence audit — summary",
        "",
        "> **Audit output only — manuscript and benchmark verdict not yet updated.**",
        "",
        f"- Judge model: `{model}`  |  Paperclip: `{version}`  |  query template: `{QUERY_TEMPLATE}`  |  top_k = {TOP_K}",
        "",
        "## Retrieval",
        f"- Candidate regulators (queries): **{n_tf}**",
        f"- Paperclip queries run: **{retr_stats.get('queries_completed')}**",
        f"- Total TF–paper result rows: **{len(log)}**",
        f"- Total unique papers (global): **{unique_papers}**",
        f"- Results with full text (source heuristic): **{ft}**",
        f"- Abstract only: **{ab_only}**",
        f"- Neither abstract nor full text: **{neither}**",
        f"- TFs with fewer than {TOP_K} results: **{len(retr_stats.get('tfs_under_10', []))}** "
        f"{retr_stats.get('tfs_under_10', [])}",
        f"- TFs with no results: **{len(retr_stats.get('tfs_zero_results', []))}** "
        f"{retr_stats.get('tfs_zero_results', [])}",
        "",
        "## Judge",
        f"- Judge parse failures: **{len(judge_stats.get('parse_failures', []))}** "
        f"{judge_stats.get('parse_failures', [])}",
        f"- Judge retries: **{judge_stats.get('retries', 0)}**",
        f"- Judge skipped (no CURSOR_API_KEY): **{judge_stats.get('skipped_no_key', False)}**",
        f"- Rubric validation failures: **{len(rubric_fail)}** {rubric_fail}",
        "",
        "## Tier distribution (final usable tier)",
        f"- strong: **{tier_counts['strong']}**",
        f"- moderate: **{tier_counts['moderate']}**",
        f"- weak: **{tier_counts['weak']}**",
        f"- none: **{tier_counts['none']}**",
        f"- missing (rubric-failed / no judge output): **{tier_counts['missing']}**",
        "",
        "## Provisional strong-or-moderate counts",
        f"- GREmLN top 25: **{_sm(gr)}** / {int(gr.sum())}",
        f"- GENIE3 top 25: **{_sm(g3)}** / {int(g3.sum())}",
        f"- Shared candidates: **{_sm(shared)}** / {int(shared.sum())}",
        f"- GREmLN-only candidates: **{_sm(gr & ~g3)}** / {int((gr & ~g3).sum())}",
        f"- GENIE3-only candidates: **{_sm(g3 & ~gr)}** / {int((g3 & ~gr).sum())}",
        "",
        "_All counts labelled: Audit output only — manuscript and benchmark verdict not yet updated._",
    ]
    (out / "audit_summary.md").write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo-root", default=None)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--all", action="store_true", help="Run the full candidate union.")
    g.add_argument("--dry", action="store_true",
                   help=f"Dry run on {DRY_TFS}.")
    g.add_argument("--tfs", default=None, help="Comma-separated TF subset.")
    ap.add_argument("--stage", default="all",
                    choices=["freeze", "retrieve", "inputs", "judge", "rubric",
                             "summary", "all"])
    ap.add_argument("--judge-model", default=DEFAULT_JUDGE_MODEL)
    ap.add_argument("--no-cache", action="store_true",
                    help="Ignore cached Paperclip/judge responses.")
    args = ap.parse_args(argv)

    here = Path(__file__).resolve().parent
    repo = Path(args.repo_root).resolve() if args.repo_root else find_repo_root(here)
    out = repo / OUT_REL
    out.mkdir(parents=True, exist_ok=True)
    use_cache = not args.no_cache
    start_ts = utcnow()

    bin_ = paperclip_bin()
    version = paperclip_version(bin_)

    # freeze always runs to define the candidate set (fast, deterministic).
    freeze_stats = stage_freeze(repo, out)
    all_tfs = freeze_stats["tfs"]
    if args.tfs:
        tfs = [t.strip() for t in args.tfs.split(",") if t.strip()]
    elif args.dry:
        tfs = [t for t in DRY_TFS if t in all_tfs] or DRY_TFS
    else:
        tfs = all_tfs

    run_meta = {"start_ts": start_ts,
                "paperclip_command_example": f"{Path(bin_).name} search -s {PAPERCLIP_SOURCE} \"<TF> {QUERY_SUFFIX}\" -n {TOP_K}"}

    retr_stats = {"queries_attempted": 0, "queries_completed": 0, "result_rows": 0,
                  "unique_papers": 0, "tfs_zero_results": [], "tfs_under_10": []}
    judge_stats = {"attempted": 0, "completed": 0, "retries": 0, "parse_failures": [],
                   "skipped_no_key": True, "results": {}}

    if args.stage in ("retrieve", "all"):
        retr_stats = stage_retrieve(repo, out, tfs, bin_, version, use_cache)
    if args.stage in ("inputs", "all"):
        stage_inputs(out, tfs)
    if args.stage in ("judge", "all"):
        judge_stats = stage_judge(repo, out, tfs, args.judge_model, use_cache)
    if args.stage in ("rubric", "summary", "all"):
        # rubric/summary need retrieval + judge context; reload if run standalone.
        if args.stage in ("rubric", "summary") and not retr_stats["result_rows"]:
            import pandas as pd
            log = pd.read_csv(out / "paperclip_retrieval_log.csv", keep_default_na=False)
            retr_stats["result_rows"] = len(log)
            retr_stats["queries_completed"] = log["TF"].nunique()
            retr_stats["queries_attempted"] = len(tfs)
        if args.stage in ("rubric", "summary") and not judge_stats["results"]:
            judge_stats = _reload_judge(out, tfs)
        stage_rubric_and_summary(repo, out, tfs, judge_stats, retr_stats,
                                 freeze_stats, args.judge_model, version, run_meta)

    print(f"\nDone. Outputs in {out}")
    return 0


def _reload_judge(out: Path, tfs: list) -> dict:
    """Reconstruct judge_stats from cached judge_raw_outputs for rubric/summary."""
    jraw_dir = out / "judge_raw_outputs"
    stats = {"attempted": 0, "completed": 0, "retries": 0, "parse_failures": [],
             "skipped_no_key": not jraw_dir.exists(), "results": {}}
    for tf in tfs:
        cf = jraw_dir / f"{tf}.json"
        if not cf.exists():
            continue
        record = json.loads(cf.read_text())
        parsed = record.get("parsed")
        if parsed is not None:
            stats["completed"] += 1
        else:
            stats["parse_failures"].append(tf)
        stats["results"][tf] = {
            "parsed": parsed,
            "output_sha256": sha256_text(json.dumps(record, sort_keys=True)),
        }
    return stats


if __name__ == "__main__":
    sys.exit(main())
