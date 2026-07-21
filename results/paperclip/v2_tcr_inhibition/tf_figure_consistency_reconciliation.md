# TF-level Paperclip figure consistency reconciliation

## Why the previous GREmLN ordering differed

The previous TF-level figure ordered GREmLN candidates by row order in `gremln_btla_vs_tcr_seed_excluded_tf_ranking.csv` (neighbour-count / CSLS score order) after filtering to the top 25 membership set. That ignored the canonical within-tie SumCSLS procedure used for Table 3: sort by `gremln_dense_rank` ascending, then `gremln_csls_seed_sum` descending. Consequently ties at the same seed-neighbour count were shown in an order that did not match the publication table (for example KIF22 ahead of EZH2 by SumCSLS in Table 3, but EZH2 first in the previous figure).

GENIE3 previously used ranking-file order; the canonical residual-tie rule is alphabetical within equal `genie3_dense_rank`. In this dataset GENIE3 top 25 order was already nearly identical.

## Corrected canonical display order (top 5)

- GREmLN: KIF22, EZH2, HMGN3, HIRIP3, STAT4
- GENIE3: EGR2, JUNB, FOXP1, HIVEP3, IRF8

## Candidates whose displayed position changed

### GREmLN
- KIF22: 2 → 1
- EZH2: 1 → 2
- HMGN3: 5 → 3
- HIRIP3: 8 → 4
- STAT4: 10 → 5
- SNAPC4: 7 → 6
- RBPJ: 12 → 7
- CREM: 6 → 8
- STAU2: 3 → 9
- MLX: 13 → 10
- AHR: 9 → 12
- NFIL3: 4 → 13
- MSC: 19 → 15
- GTF3C2: 25 → 16
- ZBTB32: 18 → 17
- TFDP1: 22 → 18
- BCL3: 20 → 19
- ZNF706: 16 → 20
- ZNF121: 24 → 21
- JUNB: 23 → 22
- PRDM1: 17 → 23
- IRF8: 15 → 24
- PLAGL2: 21 → 25

### GENIE3
- none

## Candidates whose supporting key paper count changed

### GREmLN
- CREM: 2 → 1 (now supporting = eligible ∩ relevant)
- IRF8: 3 → 0 (now supporting = eligible ∩ relevant)

### GENIE3
- IRF8: 3 → 0 (now supporting = eligible ∩ relevant)
- CREM: 2 → 1 (now supporting = eligible ∩ relevant)

## Candidates whose displayed direct (5th) count changed

Previous 5th number was broad direct causal; now it is qualifying direct (strong-tier phenotype).

### GREmLN
- EZH2: displayed 1 → 0 (broad direct retained = 1)
- RBPJ: displayed 2 → 0 (broad direct retained = 2)
- NFIL3: displayed 1 → 0 (broad direct retained = 1)
- BCL3: displayed 3 → 0 (broad direct retained = 3)

### GENIE3
- EGR2: displayed 2 → 1 (broad direct retained = 2)
- FOXP1: displayed 3 → 0 (broad direct retained = 3)
- BHLHE40: displayed 5 → 1 (broad direct retained = 5)
- ID2: displayed 1 → 0 (broad direct retained = 1)
- MAF: displayed 1 → 0 (broad direct retained = 1)

## Evidence tiers

- No audited provisional tier changed in this correction pass.

## NR4A2 exact failure

- TF: **NR4A2** (GENIE3 specific)
- Original judge tier: **strong**
- Entity eligible / relevant / primary / broad direct / qualifying direct: 8 / 3 / 2 / 2 / 2
- Supporting key papers: `PMC11439235;PMC6546093`
- Assessment-cited papers: `PMC12889511;PMC11439235;PMC6546093`
- Exact failed assertion: `the supporting excerpt for PMC7883379 exceeded the deterministic 20-word limit`
- Protocol-compliant retry: **not allowed** (retry is only for unparseable JSON; content rubric violations are frozen as missing).
- Figure audit label: **Supporting excerpt exceeded word limit**
- Final tier: remains **missing**.

## Focus TF direct-evidence audit

Whether previously displayed direct counts were broad causal rather than strong-tier qualifying:

- **BHLHE40** (weak): broad direct = 5, qualifying direct = 1 — previously displayed broad causal as if strong-qualifying.
- **RBPJ** (weak): broad direct = 2, qualifying direct = 0 — previously displayed broad causal as if strong-qualifying.
- **ID2** (weak): broad direct = 1, qualifying direct = 0 — previously displayed broad causal as if strong-qualifying.
- **MAF** (weak): broad direct = 1, qualifying direct = 0 — previously displayed broad causal as if strong-qualifying.
- **JUNB** (moderate): broad direct = 2, qualifying direct = 2 — counts agree.
- **EZH2** (moderate): broad direct = 1, qualifying direct = 0 — previously displayed broad causal as if strong-qualifying.
- **NFIL3** (moderate): broad direct = 1, qualifying direct = 0 — previously displayed broad causal as if strong-qualifying.
- **AHR** (moderate): broad direct = 1, qualifying direct = 1 — counts agree.
- **BCL3** (moderate): broad direct = 3, qualifying direct = 0 — previously displayed broad causal as if strong-qualifying.
- **FOXP1** (moderate): broad direct = 3, qualifying direct = 0 — previously displayed broad causal as if strong-qualifying.

## Aggregate chart reconciliation

- GREmLN: {'strong': 1, 'moderate': 6, 'weak': 3, 'none': 15, 'missing': 0} (sum=25)
- GENIE3: {'strong': 2, 'moderate': 6, 'weak': 6, 'none': 10, 'missing': 1} (sum=25)

Matches `paperclip_evidence_summary_counts.csv`.
