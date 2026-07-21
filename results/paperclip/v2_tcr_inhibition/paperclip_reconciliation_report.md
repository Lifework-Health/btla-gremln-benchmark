# Paperclip v2 — reconciliation report
_Literature annotation only; the benchmark verdict is NOT changed in this task._

## Original vs corrected tier distribution
| Tier | Original | Corrected (after entity validation) |
|---|---|---|
| strong | 3 | 2 |
| moderate | 10 | 10 |
| weak | 7 | 7 |
| none | 22 | 23 |
| missing | 1 | 1 |
| **sum** | 43 | 43 |

## Candidates whose tier changed
| TF | group | legacy | corrected | reason |
|---|---|---|---|---|
| MSC | gremln_only | strong | none | entity-filtered rerun tier=none; rubric passed |

## Missing candidate
- **NR4A2** (genie3_only).
  - Original judge tier: **strong**; deterministic rubric FAILED: `supporting_excerpt >20 words for PMC7883379`.
  - Cause class: **rubric failure** (a supporting excerpt exceeded the 20-word limit) — not retrieval failure, not invalid JSON, not a missing entity match. Entity audit is clean (entity-eligible papers 8, direct-causal 2).
  - Retry under existing protocol: **not allowed**. The protocol retries once only on unparseable JSON; a content rubric violation is not a retry trigger, so re-running would change the frozen protocol.
  - Outcome: **remains missing** after the audit; not converted to none.

## Corrected strong+moderate counts
| Group | strong+moderate |
|---|---|
| GREmLN top 25 | 7 |
| GENIE3 top 25 | 8 |
| GREmLN specific (n=18) | 4 |
| GENIE3 specific (n=18) | 5 |
| Shared (n=7) | 3 |

## Literature comparison outcome
- Model-specific strong+moderate: GREmLN 4 vs GENIE3 5 -> literature support **leans GENIE3**.
- The benchmark corroboration rule compares model-specific candidates (shared candidates count for both). This literature annotation is not a predictive input and does not, by itself, change the CRISPRi-anchored verdict.

## Would the benchmark verdict change?
- Not in this task. The verdict is anchored on CRISPRi causal validation; Paperclip tiers are annotation only. The entity-corrected literature picture leans GENIE3 on model-specific candidates, which is recorded here for transparency but is explicitly not applied to the verdict, Table 4, Results or Conclusions.
