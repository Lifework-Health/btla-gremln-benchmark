"""Unit tests for the auditable Paperclip TCR-inhibition pipeline.

These cover the deterministic, side-effect-free parts of the pipeline:
query construction, top-10 truncation, JSON extraction, schema validation and
the post-judge rubric checks (including every rubric edge case required by the
task). They do not call Paperclip or the LLM judge.

Run:  python -m pytest tests/test_paperclip_tcr_inhibition_audit.py -q
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import run_paperclip_tcr_inhibition_audit as A  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
SCHEMA = json.loads((REPO / A.SCHEMA_REL).read_text())


# --------------------------------------------------------------------------- #
# query construction + truncation
# --------------------------------------------------------------------------- #
def test_build_query_exact():
    assert A.build_query("EGR2") == "EGR2 TCR inhibition"
    assert A.build_query("AHR") == "AHR TCR inhibition"
    assert A.build_query("KIF22") == "KIF22 TCR inhibition"


def test_build_query_no_alias_expansion():
    # symbol used exactly as supplied; no extra terms appended
    q = A.build_query("BHLHE40")
    assert q == "BHLHE40 TCR inhibition"
    assert "CD4" not in q and "checkpoint" not in q and "BTLA" not in q


def test_truncate_top_10():
    assert A.truncate_top_k(list(range(25)), 10) == list(range(10))


def test_truncate_fewer_than_10():
    assert A.truncate_top_k([1, 2, 3], 10) == [1, 2, 3]


def test_truncate_zero_results():
    assert A.truncate_top_k([], 10) == []


def test_duplicate_papers_across_tfs_not_deduped():
    # the pipeline keeps per-TF rows; simulate two TFs sharing a paper
    tf_a = A.truncate_top_k(["PMC1", "PMC2"], 10)
    tf_b = A.truncate_top_k(["PMC1", "PMC9"], 10)
    combined = [("A", p) for p in tf_a] + [("B", p) for p in tf_b]
    assert combined.count(("A", "PMC1")) == 1
    assert combined.count(("B", "PMC1")) == 1  # visible under each query
    assert len({p for _, p in combined}) == 3  # unique-global count


# --------------------------------------------------------------------------- #
# JSON extraction
# --------------------------------------------------------------------------- #
def test_extract_json_plain():
    assert A.extract_json('{"a": 1}') == {"a": 1}


def test_extract_json_fenced_with_prose():
    txt = "Here is the result:\n```json\n{\"tier\": \"weak\"}\n```\nThanks."
    assert A.extract_json(txt) == {"tier": "weak"}


def test_extract_invalid_json_returns_none():
    assert A.extract_json("not json at all") is None
    assert A.extract_json("{oops not valid}") is None


# --------------------------------------------------------------------------- #
# schema validation
# --------------------------------------------------------------------------- #
def _valid_judge(tf="EGR2", tier="none", assessments=None, keys=None,
                 retrieved=1, relevant=0, query=None):
    return {
        "tf": tf,
        "query": query if query is not None else A.build_query(tf),
        "papers_retrieved": retrieved,
        "papers_relevant": relevant,
        "paper_assessments": assessments if assessments is not None else [],
        "overall_evidence_tier": tier,
        "key_paperclip_result_ids": keys or [],
        "tier_rationale": "ok",
        "confidence": "low",
    }


def _assessment(rid="PMC1", rank=1, relevant=True,
                study_type="primary experimental",
                directness="direct causal or perturbational",
                direction="promotes inhibition", excerpt="short excerpt", phenotype="x"):
    return {
        "paperclip_result_id": rid, "retrieval_rank": rank, "relevant": relevant,
        "study_type": study_type, "experimental_system": "human primary T cells",
        "evidence_directness": directness, "phenotype": phenotype,
        "direction": direction, "supporting_excerpt": excerpt, "rationale": "because",
    }


def test_schema_accepts_valid():
    assert A.validate_judge_json(_valid_judge(), SCHEMA) == []


def test_schema_rejects_bad_enum():
    obj = _valid_judge(tier="excellent")
    assert A.validate_judge_json(obj, SCHEMA)


def test_schema_rejects_missing_field():
    obj = _valid_judge()
    del obj["confidence"]
    assert A.validate_judge_json(obj, SCHEMA)


# --------------------------------------------------------------------------- #
# rubric checks
# --------------------------------------------------------------------------- #
def test_rubric_strong_requires_direct_causal_primary():
    # strong but only indirect mechanistic -> fail
    a = _assessment(directness="indirect mechanistic")
    j = _valid_judge(tier="strong", assessments=[a], keys=["PMC1"], retrieved=1)
    ok, reasons = A.run_rubric_checks("EGR2", j, ["PMC1"], "EGR2 TCR inhibition", 1)
    assert not ok
    assert any("direct-causal primary" in r for r in reasons)


def test_rubric_strong_valid_with_direct_causal_primary():
    a = _assessment()
    j = _valid_judge(tier="strong", assessments=[a], keys=["PMC1"], retrieved=1)
    ok, reasons = A.run_rubric_checks("EGR2", j, ["PMC1"], "EGR2 TCR inhibition", 1)
    assert ok, reasons


def test_rubric_moderate_requires_primary():
    a = _assessment(study_type="review", directness="indirect mechanistic")
    j = _valid_judge(tier="moderate", assessments=[a], keys=["PMC1"], retrieved=1)
    ok, reasons = A.run_rubric_checks("EGR2", j, ["PMC1"], "EGR2 TCR inhibition", 1)
    assert not ok
    assert any("primary study" in r for r in reasons)


def test_rubric_review_only_cannot_exceed_weak():
    a = _assessment(study_type="review", directness="indirect mechanistic")
    j = _valid_judge(tier="strong", assessments=[a], keys=["PMC1"], retrieved=1)
    ok, reasons = A.run_rubric_checks("EGR2", j, ["PMC1"], "EGR2 TCR inhibition", 1)
    assert not ok
    assert any("review-only" in r for r in reasons)


def test_rubric_missing_cited_paper_id():
    a = _assessment(rid="PMC_UNKNOWN")
    j = _valid_judge(tier="strong", assessments=[a], keys=["PMC_UNKNOWN"], retrieved=1)
    ok, reasons = A.run_rubric_checks("EGR2", j, ["PMC1"], "EGR2 TCR inhibition", 1)
    assert not ok
    assert any("unknown result id" in r or "not in retrieval log" in r for r in reasons)


def test_rubric_zero_papers_must_be_none():
    j = _valid_judge(tier="weak", assessments=[], keys=[], retrieved=0)
    ok, reasons = A.run_rubric_checks("EGR2", j, [], "EGR2 TCR inhibition", 0)
    assert not ok
    assert any("zero retrieved papers" in r for r in reasons)


def test_rubric_zero_papers_none_is_valid():
    j = _valid_judge(tier="none", assessments=[], keys=[], retrieved=0)
    ok, reasons = A.run_rubric_checks("EGR2", j, [], "EGR2 TCR inhibition", 0)
    assert ok, reasons


def test_rubric_excerpt_over_20_words_fails():
    long_excerpt = " ".join(["w"] * 21)
    a = _assessment(excerpt=long_excerpt)
    j = _valid_judge(tier="strong", assessments=[a], keys=["PMC1"], retrieved=1)
    ok, reasons = A.run_rubric_checks("EGR2", j, ["PMC1"], "EGR2 TCR inhibition", 1)
    assert not ok
    assert any(">20 words" in r for r in reasons)


def test_rubric_wrong_query_fails():
    j = _valid_judge(tier="none", query="EGR2 exhaustion")
    ok, reasons = A.run_rubric_checks("EGR2", j, [], "EGR2 TCR inhibition", 0)
    assert not ok
    assert any("exact query" in r for r in reasons)


def test_rubric_papers_retrieved_mismatch_fails():
    a = _assessment()
    j = _valid_judge(tier="weak", assessments=[a], keys=[], retrieved=5)
    ok, reasons = A.run_rubric_checks("EGR2", j, ["PMC1"], "EGR2 TCR inhibition", 1)
    assert not ok
    assert any("retrieval log rows" in r for r in reasons)


def test_rubric_strong_without_key_ids_fails():
    a = _assessment()
    j = _valid_judge(tier="strong", assessments=[a], keys=[], retrieved=1)
    ok, reasons = A.run_rubric_checks("EGR2", j, ["PMC1"], "EGR2 TCR inhibition", 1)
    assert not ok
    assert any("no key_paperclip_result_ids" in r for r in reasons)


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
