# Paperclip v2 audit summary (corrected)

Literature annotation only. Manuscript, Table 4 and benchmark verdict unchanged.

## Entity identity gate (mutually exclusive)

| entity_match | n |
|---|---|
| exact | 232 |
| recognised_alias | 10 |
| wrong_entity | 10 |
| ambiguous_acronym | 0 |
| not_mentioned | 178 |
| **total** | **430** |
| entity eligible (exact + recognised_alias) | 242 |
| excluded (wrong_entity + ambiguous_acronym + not_mentioned) | 188 |

Example collision: MSC = Musculin, not mesenchymal stromal cells.

## Corrected final usable tiers

| Tier | n |
|---|---|
| strong | 2 |
| moderate | 10 |
| weak | 7 |
| none | 23 |
| missing | 1 |
| **sum** | **43** |

## Missing GENIE3 specific candidate

- **NR4A2** (genie3_only): rubric failure (`supporting_excerpt >20 words for PMC7883379`). Rerun allowed: **no**. Remains **missing**.

## Human verification

- Status: **pending**. Named human reviewer has not examined the strong or moderate candidates. See `human_verification.csv`.
- LLM / automated structural checks are **not** human verification.

## Entity-filtered LLM judgements

- Claim that all 43 final judgements received only entity-eligible papers: **not confirmed** (1/43 used an eligible-only packet under the selective-rerun protocol).
- Zero entity-eligible TFs (7): GTF2A2, HIVEP3, HMGN3, MSC, SMAP2, SNAPC4, ZNF121.
- Details: `entity_filtered_judgement_status.csv`.
