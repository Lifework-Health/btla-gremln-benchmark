# Auditable Paperclip TCR-inhibition evidence (v2)

> **Audit output only — the manuscript, Results/Conclusions, Table 4, figures,
> literature counts and the superiority verdict are NOT updated by this build.**
> The legacy manual review in `../paperclip_union_top25_review.csv` is retained
> unchanged as a legacy analysis.

This directory holds a fully auditable, automated re-run of the literature
evidence layer for the BTLA GREmLN-vs-GENIE3 benchmark. Every candidate
regulator is judged from a **single fixed Paperclip query**:

```
<TF> TCR inhibition
```

with top-10 retrieval, per-paper provenance, and a fixed LLM judge whose raw
outputs are cached and validated against a JSON schema and a deterministic
rubric.

## How to reproduce

```bash
# from the repository root
export CURSOR_API_KEY="cursor_..."          # required for the judge stage only

# 1. dry run (EGR2, AHR, KIF22)
python scripts/run_paperclip_tcr_inhibition_audit.py --dry --stage all

# 2. full candidate union (43 TFs)
python scripts/run_paperclip_tcr_inhibition_audit.py --all --stage all

# stages can also be run individually: freeze | retrieve | inputs | judge | rubric | summary
# reruns reuse cached Paperclip responses and cached judge outputs (use --no-cache to refetch)
```

Unit tests (no network / no judge):

```bash
python -m pytest tests/test_paperclip_tcr_inhibition_audit.py -q
```

## Pipeline

| Stage | What it does |
|-------|--------------|
| `freeze`   | Candidate union = top-25 of each model from the canonical **seed-excluded** rankings (`results/tables/{gremln,genie3}_btla_vs_tcr_seed_excluded_tf_ranking.csv`). Membership is taken from the canonical `top25_union_primary.csv` so the frozen set is byte-identical to the published Table 3/4 candidates; the GREmLN top-25 boundary tie is resolved by the summed-CSLS tie-break already applied in nb03. |
| `retrieve` | One `paperclip search -s pmc "<TF> TCR inhibition" -n 10` per TF. No aliases, no extra terms, no second pass, no manual supplementation. Per-paper metadata is enriched from `/papers/<id>/meta.json`. Irrelevant results are retained. Papers are **not** de-duplicated across TFs. |
| `inputs`   | One evidence packet per TF (`judge_inputs/<TF>.json`) containing only Paperclip-supplied fields. |
| `judge`    | One fixed LLM judge call per TF via the **Cursor SDK** (`cursor_sdk`, local `Agent.prompt`). Prompt: `prompts/paperclip_tcr_inhibition_judge_v1.txt`. Raw outputs saved to `judge_raw_outputs/<TF>.json`; invalid JSON is retried once. |
| `rubric`   | Deterministic post-judge validation. Distinguishes `judge_tier`, `rubric_valid` and `final_usable_tier` (set to `missing` when the rubric fails — no manual correction). |
| `summary`  | Emits `tf_evidence_tiers.csv`, `paper_judgements.csv`, `run_manifest.json`, `audit_summary.md`, and the `legacy_vs_audited_tiers.csv` diagnostic. |

## Files

| File | Committed? | Contents |
|------|-----------|----------|
| `candidate_union.csv`         | yes | frozen 43-TF union + ranks, group, ranking provenance (commit + SHA256). |
| `paperclip_retrieval_log.csv` | yes | one row per TF–paper, full retrieval provenance. |
| `paperclip_raw_manifest.csv`  | yes | identifiers, byte sizes and SHA256 of every raw response. |
| `judge_inputs/<TF>.json`      | yes | exact judge input packets. |
| `judge_raw_outputs/<TF>.json` | yes | unmodified judge outputs (all attempts) + parsed object. |
| `paper_judgements.csv`        | yes | one row per TF–paper judge assessment. |
| `tf_evidence_tiers.csv`       | yes | one row per TF: judge_tier / rubric_valid / final_usable_tier + hashes. |
| `run_manifest.json`           | yes | full run provenance (models, hashes, counts, timestamps). |
| `audit_summary.md`            | yes | human-readable audit report. |
| `legacy_vs_audited_tiers.csv` | yes | diagnostic only; **not** used to update the manuscript. |
| `raw/`                        | **no (git-ignored)** | full raw Paperclip responses / abstracts / snippets and the per-TF response cache. Only the SHA256 manifest is committed. |

## Judge determinism note

Cursor SDK agent models do not expose `temperature`, `top_p` or `seed`, so these
are recorded as `not_supported` in `run_manifest.json`. Reproducibility is
therefore provided by **caching**: reruns reuse `judge_raw_outputs/<TF>.json`
(and the git-ignored Paperclip response cache) and do not change results unless
`--no-cache` is passed. The judge model identifier, prompt SHA256 and schema
SHA256 are all recorded.

## Source scope

The Paperclip query string is exactly `<TF> TCR inhibition`. `-s pmc` is a
required **source scope** flag (PubMed Central full-text papers), not a query
term; it is recorded in `run_manifest.json` as `paperclip_source_scope`.
