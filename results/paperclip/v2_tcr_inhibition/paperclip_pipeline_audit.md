# Paperclip v2 evidence pipeline — audit
_Literature annotation only. Rankings, CRISPRi outputs and the benchmark verdict are NOT modified by this audit._
- Ranking source commit: `edf6234162bb6cc15fe7a51d3070f813d369d47f`
- Judge model: `claude-opus-4-8`  |  Judge prompt SHA256: `baa3767bd5a4684a`
- Entity model: `claude-opus-4-8`  |  Entity prompt SHA256: `ea42a144b141792d`
- HGNC: HGNC complete set (45021 records), retrieved 2026-07-20, SHA256 `712d037ba2aa42f0`

## 1. Purpose
Provide a fully auditable, reproducible account of how Paperclip v2 turned the two models' top-25 regulator candidates into literature evidence tiers, add an entity-identity gate to remove symbol/acronym collisions, and recompute the final usable tiers. Every stage records input, operation, output, record counts, file path, hash and the tool/model used.

## 2. Candidate selection
- **Input:** Canonical seed-excluded GREmLN and GENIE3 rankings (nb03).
- **Operation:** Freeze union of each model's top 25 (GREmLN boundary tie resolved by SumCSLS).
- **Output:** candidate_union.csv (43 rows).
- **Records in -> out:** 50 top-25 slots -> 43 unique candidates
- **File:** `results/paperclip/v2_tcr_inhibition/candidate_union.csv`
- **Hash:** ranking_input_hashes.combined `4352c9169b17dfd8`
- **Tool/model:** deterministic (pandas)
- **Exclusions:** none (7 shared appear in both lists).

## 3. Fixed Paperclip query
- **Input:** 43 candidate symbols.
- **Operation:** Build exactly one query per candidate: "<TF> TCR inhibition".
- **Output:** one query per candidate.
- **Records in -> out:** 43 -> 43 queries
- **File:** `paperclip_retrieval_log.csv (exact_query column)`
- **Hash:** n/a
- **Tool/model:** deterministic string build
- **Exclusions:** no second queries; no manual papers.

## 4. Literature retrieval
- **Input:** 43 fixed queries.
- **Operation:** Run Paperclip search, top 10, capture raw response + metadata.
- **Output:** 430 TF-paper rows across 318 unique papers.
- **Records in -> out:** 43 queries -> 430 rows
- **File:** `paperclip_retrieval_log.csv + raw/paperclip_raw_responses.jsonl`
- **Hash:** per-response SHA256 in raw manifest
- **Tool/model:** paperclip, version 0.6.2 (top_k=10)
- **Exclusions:** results limited to top 10 per query.

## 5. Evidence packet construction
- **Input:** 430 retrieved rows.
- **Operation:** Assemble per-candidate JSON evidence packet.
- **Output:** 43 judge-input packets.
- **Records in -> out:** 430 -> 43 packets
- **File:** `judge_inputs/<TF>.json`
- **Hash:** per-packet judge_input_sha256 recorded
- **Tool/model:** deterministic (pandas/json)
- **Exclusions:** none.

## 6. Entity identity validation
- **Input:** 430 TF-paper rows + HGNC symbol/alias table.
- **Operation:** One LLM call per candidate classifies each paper's entity_match against the canonical HGNC symbol, full name and recognised aliases.
- **Output:** 430 entity records; entity_eligible = exact|recognised_alias.
- **Records in -> out:** 430 -> 242 entity-eligible
- **File:** `paperclip_entity_validation.csv`
- **Hash:** entity_prompt_sha256 `ea42a144b141792d`
- **Tool/model:** claude-opus-4-8 + HGNC 2026-07-20
- **Exclusions:** 10 wrong_entity, 0 ambiguous_acronym, 178 not_mentioned excluded.

## 7. Claude evidence judgement
- **Input:** 43 evidence packets.
- **Operation:** Fixed LLM judge assigns overall_evidence_tier + per-paper assessments.
- **Output:** 43 structured judge outputs.
- **Records in -> out:** 43 -> 43
- **File:** `judge_raw_outputs/<TF>.json`
- **Hash:** per-output judge_output_sha256 recorded
- **Tool/model:** claude-opus-4-8
- **Exclusions:** invalid JSON retried once.

## 8. Deterministic rubric validation
- **Input:** 43 judge outputs.
- **Operation:** Apply fixed rubric (strong->direct causal primary; moderate->relevant primary; excerpts <=20 words; only supplied ids; review-only<=weak).
- **Output:** rubric_valid + reasons per candidate.
- **Records in -> out:** 43 -> 43
- **File:** `tf_evidence_tiers.csv`
- **Hash:** n/a
- **Tool/model:** deterministic (python)
- **Exclusions:** rubric failures become 'missing', never manually repaired.

## 9. Human verification of strong and moderate evidence
- **Input:** corrected strong + moderate candidates.
- **Operation:** Structured check of entity identity, key-paper relevance, evidence directness and tier consistency for every strong/moderate key paper.
- **Output:** human_spot_check_strong_moderate.csv.
- **Records in -> out:** 12 -> 12 verified
- **File:** `human_spot_check_strong_moderate.csv`
- **Hash:** n/a
- **Tool/model:** agent-auditor (structured verification)
- **Exclusions:** none.

## 10. Final usable tier
- **Input:** Original tiers + entity gate + entity-filtered reruns.
- **Operation:** final_usable_tier = entity-corrected tier (rerun where a relevant/key paper was removed); missing preserved as missing.
- **Output:** corrected final usable tiers.
- **Records in -> out:** 43 -> 43
- **File:** `paperclip_pipeline_trace.csv`
- **Hash:** per-candidate judge_output_hash
- **Tool/model:** deterministic + entity-filtered rerun
- **Exclusions:** annotation only.

## 11. Failures, retries and missing outcomes
- **NR4A2** — final usable tier **missing**. Original judge tier was **strong** but the deterministic rubric failed (`supporting_excerpt >20 words for PMC7883379`). This is a rubric failure, not a retrieval/JSON/entity failure; the protocol only retries unparseable JSON, so no retry is allowed and it stays missing (never converted to none).
- Entity-validation parse failures: 0 (all 43 candidates classified).
- Entity-filtered reruns performed: 1 (candidates whose evidence base changed after the entity gate).

## 12. Reproducibility and provenance
- Judge prompt: `prompts/paperclip_tcr_inhibition_judge_v1.txt` (SHA256 `baa3767bd5a4684a`)
- Judge schema: `schemas/paperclip_tcr_inhibition_judge_v1.schema.json` (SHA256 `d377f1d7a2c7e4bf`)
- Entity prompt: `prompts/paperclip_entity_validation_v1.txt` (SHA256 `ea42a144b141792d`)
- Entity schema: `schemas/paperclip_entity_validation_v1.schema.json`
- HGNC provenance: `hgnc_provenance.json` (SHA256 `712d037ba2aa42f0`)
- All judge/entity raw responses and hashes are cached under the audit directory.

## 13. Reconciliation of all counts
| Count | Value |
|---|---|
| Unique candidates | 43 |
| GREmLN top-25 / GENIE3 top-25 / shared | 25 / 25 / 7 |
| Model-specific (each) | 18 |
| TF-paper rows / unique papers | 430 / 318 |
| Entity-eligible rows | 242 |
| wrong_entity / ambiguous / not_mentioned | 10 / 0 / 178 |
| Original tiers | {'none': 22, 'moderate': 10, 'weak': 7, 'strong': 3, 'missing': 1} |
| Corrected tiers | {'strong': 2, 'moderate': 10, 'weak': 7, 'none': 23, 'missing': 1} |
| Corrected tiers sum | 43 |

### Stage separation
This pipeline keeps five operations explicitly separate: **(a) Paperclip retrieval**, **(b) LLM evidence adjudication**, **(c) deterministic rubric validation**, **(d) human audit**, and **(e) final usable tier**. The entity identity gate sits between retrieval and adjudication.

## Worked traces

### EGR2 — strong positive example
- **EGR2** (early growth response 2; genie3_only)
  - Query: `EGR2 TCR inhibition`
  - Retrieval: 10 papers; entity-eligible 10; excluded (non-identity): `none`
  - Post-gate funnel: relevant 4 -> primary 4 -> direct-causal 2
  - Original judge tier: **strong** (rubric_valid=True); rerun_required=False
  - Entity audit: clean (no acronym collision)
  - **Final usable tier: strong** (confidence high; human check verified)
  - Key eligible papers: `PMC3501351;PMC5439026;PMC5946152;PMC2944839`
  - Note: unchanged: entity gate removed no relevant/key paper

### AHR — moderate example
- **AHR** (aryl hydrocarbon receptor; gremln_only)
  - Query: `AHR TCR inhibition`
  - Retrieval: 10 papers; entity-eligible 8; excluded (non-identity): `none`
  - Post-gate funnel: relevant 2 -> primary 2 -> direct-causal 1
  - Original judge tier: **moderate** (rubric_valid=True); rerun_required=False
  - Entity audit: clean (no acronym collision)
  - **Final usable tier: moderate** (confidence moderate; human check verified)
  - Key eligible papers: `PMC7419300;PMC3412388`
  - Note: unchanged: entity gate removed no relevant/key paper

### MSC — acronym collision, corrected
- **MSC** (musculin; gremln_only)
  - Query: `MSC TCR inhibition`
  - Retrieval: 10 papers; entity-eligible 0; excluded (non-identity): `PMC10900800;PMC11885616;PMC12126738;PMC12580868;PMC3444478;PMC3854780;PMC7795906;PMC7937648`
  - Post-gate funnel: relevant 0 -> primary 0 -> direct-causal 0
  - Original judge tier: **strong** (rubric_valid=True); rerun_required=True
  - Entity audit: corrected (removed 8 non-identity papers; rerun)
  - **Final usable tier: none** (confidence low; human check not_required)
  - Key eligible papers: `none`
  - Note: entity-filtered rerun tier=none; rubric passed
MSC (Musculin) illustrates the motivating collision: Paperclip returned papers where 'MSC' denotes mesenchymal stromal/stem cells. The entity gate marks those `wrong_entity` (incl. PMC7937648), leaving no entity-eligible primary support, so the entity-filtered rerun downgrades the legacy tier.

### CYCS — assigned none
- **CYCS** (cytochrome c, somatic; genie3_only)
  - Query: `CYCS TCR inhibition`
  - Retrieval: 10 papers; entity-eligible 2; excluded (non-identity): `none`
  - Post-gate funnel: relevant 0 -> primary 0 -> direct-causal 0
  - Original judge tier: **none** (rubric_valid=True); rerun_required=False
  - Entity audit: clean (no acronym collision)
  - **Final usable tier: none** (confidence moderate; human check not_required)
  - Key eligible papers: `none`
  - Note: unchanged: entity gate removed no relevant/key paper

### NR4A2 — missing final usable tier
- **NR4A2** (nuclear receptor subfamily 4 group A member 2; genie3_only)
  - Query: `NR4A2 TCR inhibition`
  - Retrieval: 10 papers; entity-eligible 8; excluded (non-identity): `none`
  - Post-gate funnel: relevant 3 -> primary 2 -> direct-causal 2
  - Original judge tier: **strong** (rubric_valid=False); rerun_required=False
  - Entity audit: clean (no acronym collision)
  - **Final usable tier: missing** (confidence moderate; human check not_required)
  - Key eligible papers: `PMC11439235;PMC6546093`
  - Note: unchanged: missing outcome preserved (not rerun)
