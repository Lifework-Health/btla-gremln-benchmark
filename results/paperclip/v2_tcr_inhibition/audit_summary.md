# Paperclip TCR-inhibition evidence audit — summary

> **Audit output only — manuscript and benchmark verdict not yet updated.**

- Judge model: `claude-opus-4-8`  |  Paperclip: `paperclip, version 0.6.2`  |  query template: `{TF} TCR inhibition`  |  top_k = 10

## Retrieval
- Candidate regulators (queries): **43**
- Paperclip queries run: **43**
- Total TF–paper result rows: **430**
- Total unique papers (global): **318**
- Results with full text (source heuristic): **400**
- Abstract only: **30**
- Neither abstract nor full text: **0**
- TFs with fewer than 10 results: **0** []
- TFs with no results: **0** []

## Judge
- Judge parse failures: **0** []
- Judge retries: **0**
- Judge skipped (no CURSOR_API_KEY): **False**
- Rubric validation failures: **1** ['NR4A2']

## Tier distribution (final usable tier)
- strong: **3**
- moderate: **10**
- weak: **7**
- none: **22**
- missing (rubric-failed / no judge output): **1**

## Provisional strong-or-moderate counts
- GREmLN top 25: **8** / 25
- GENIE3 top 25: **8** / 25
- Shared candidates: **3** / 7
- GREmLN-only candidates: **5** / 18
- GENIE3-only candidates: **5** / 18

_All counts labelled: Audit output only — manuscript and benchmark verdict not yet updated._
