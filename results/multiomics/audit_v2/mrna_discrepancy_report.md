# mRNA discrepancy report — old 08e GENIE3 chart vs audit_v2

_Evidence reconciliation only. Manuscript, Table 4 and benchmark verdict
unchanged. Restricted log2FC / adj-p values are in
`local_restricted/mrna_reconciliation_values.RESTRICTED.csv`; this report and
`mrna_old_vs_current_reconciliation.csv` carry statuses only._

## 1. Exact reconstruction of the old chart

| item | finding |
|---|---|
| Generating code | `notebooks/BTLA_CSLS_MultiOmics/08e_BTLA_candidate_multiomics_evidence.ipynb`, functions `build_multiomics_matrix` + `plot_multiomics_matrix`, via `_compact_layer_evidence(gene, "transcriptomics", contrast=None)` |
| Candidate list | `results/btla_csls_multiomics/btla_tf_genie3_top25.csv` (25 rows) |
| Seed inclusion | **Seed-inclusive**. Contains 5 BTLA_vs_TCR seeds (BACH2, NR4A3, POU2AF1, REL, TBX21). The chart marks seeds with `*`. |
| Ranking version / date | 08-series GENIE3 top-25 (columns include `genie3_score_hybrid`, `in_GREmLN_top25`); file mtime 2026-06-24, committed 2026-06-28 (`c585c6f`, "Add BTLA CSLS Multi-Omics notebook suite 08a–08g"). This predates the canonical seed-excluded ranking used by the benchmark. |
| Source workbook / sheet | `data/multiomics/Transcriptomics/@Transcriptomics/Tables/DEG_wide_statistics.xlsx`, `Sheet 1` |
| Contrast(s) used | **Both** `BTLAvsTCR` and `TCRvsUC` — the mRNA cell was built with `contrast=None`, so all transcriptomics rows for a gene (both contrasts) were pooled. |
| Timepoints | 1h, 4h, 24h (all pooled) |
| Columns selected | `BTLAvsTCR_{tp}_log2FoldChange`/`_padj` and `TCRvsUC_{tp}_log2FoldChange`/`_padj` |
| Thresholds | `support_from_p`: qualifying = padj<0.05 AND |log2FC|≥log2(1.5); inner "strong" = padj<0.01 AND |log2FC|≥1.5·log2(1.5). **Same source-native cut as audit_v2** (not a threshold difference). |
| Glyph meaning | `++` = strong tier; `·` = moderate tier (both qualifying); `ns` = measured but not significant; `—` = no row / not measured. Label also appended timepoint + ↑/↓ direction. |
| Selection across contrasts/timepoints | **Yes** — `hits` = all qualifying rows (both contrasts, all timepoints), sorted by tier then adjusted p, and the single strongest was displayed. |

## 2. Root cause

**The apparent "10–15 mRNA signals" are 12 positive glyphs, and every one of
them came from `TCRvsUC` (generic TCR-vs-unstimulated activation), not the
BTLA-specific `BTLAvsTCR` contrast.** In the frozen long-evidence table the
transcriptomics qualifying rows are: `BTLAvsTCR` = 8 (3 moderate + 5 strong,
spread over genes/timepoints), `TCRvsUC` = 40 (9 moderate + 31 strong). Because
the chart pooled both contrasts and displayed the strongest, the BTLA column was
populated by activation-driven `TCRvsUC` changes.

Old positive calls (all 12, all `TCR_vs_UC`):
BACH2, EGR2, FOSL2, HIVEP3, ID2, JUN, JUNB, NR4A3, POU2AF1, PRDM1, REL, TBX21.

## 3. Decomposition of the difference

- **Contrast pooling (explanations A + B): dominant cause.** 12/12 old positives
  are `TCRvsUC`; 0/12 are `BTLAvsTCR`. Restricting to the BTLA-specific contrast
  alone removes all 12.
- **Seed-inclusive candidate universe (explanation C): secondary.** The old list
  had 25 seed-inclusive candidates; 9 are absent from the current 43 seed-excluded
  union — BACH2, NR4A3, POU2AF1, REL, TBX21 (BTLA_vs_TCR seeds) plus JUN, FOSL2,
  SATB1, ARID5B. **7 of the 12 old positives** sit on candidates not in the current
  union (BACH2, FOSL2, JUN, NR4A3, POU2AF1, REL, TBX21). The other 5 (EGR2, HIVEP3,
  ID2, JUNB, PRDM1) are in the union but their positives were `TCRvsUC`.
- **Threshold (D, E): not a factor.** Old chart and audit_v2 use the same
  source-native cut (padj<0.05, |log2FC|>log2(1.5)).
- **audit_v2 bug (F): none.** audit_v2's BTLA-specific mRNA result was
  independently recomputed from the raw workbook (section 3 below) and matches.
- **Old-chart bug (G): not a parsing bug, a scoping flaw.** Gene mapping and
  column parsing were correct; the defect is conceptual — pooling `TCRvsUC` into a
  panel presented as BTLA multiomics evidence, so generic activation was mislabelled
  as candidate mRNA support for the BTLA contrast.

## 4. Focus candidates (shared between old chart and current union)

All 15 are `not_seed` and in the current union. Old mRNA glyph → current
BTLA-specific result:

| TF | old glyph | old contrast | current BTLA-specific | note |
|---|---|---|---|---|
| JUNB | `++` 1h | TCR_vs_UC | no qualifying | old positive was TCRvsUC |
| HIVEP3 | `·` 24h | TCR_vs_UC | no qualifying (nearest miss, gap 0.01) | old positive was TCRvsUC |
| EGR2 | `++` 1h | TCR_vs_UC | no qualifying | old positive was TCRvsUC |
| ID2 | `++` 4h | TCR_vs_UC | no qualifying | old positive was TCRvsUC |
| PRDM1 | `++` 4h | TCR_vs_UC | no qualifying | old positive was TCRvsUC |
| IRF8, STAT1, CREM, MYC, ZBTB32, MDM2, TRAF4, ZEB2, BCL3, MAF | `ns` / `—` | — | no qualifying | already non-significant in old chart too |

(Exact log2FC/padj per contrast/timepoint for these are in the restricted values file.)

## 5. Verification of the current 43-candidate BTLA-specific mRNA (recomputed from raw)

Source rule: `padj < 0.05 AND |log2FC| > log2(1.5)`, `BTLAvsTCR` only, per timepoint 1h/4h/24h.

| metric | count (of 43) |
|---|---|
| measured (≥1 timepoint with effect+padj) | 41 |
| not measured | 0 |
| absent from assay universe | 2 (AHR, MSC) |
| adjusted-p significant (any timepoint) | 1 |
| effect-threshold passed (any timepoint) | 8 |
| **passing both (qualifying)** | **0** |

Nearest miss: **HIVEP3**, gap 0.01 in |log2FC| at its padj-significant timepoint
(the single padj-significant candidate; effect just below log2(1.5)). No other
candidate is both padj-significant and effect-near-threshold.

Assertions (all pass): all 43 accounted for; no `TCRvsUC` column enters the
BTLA-specific result; the 2 absent candidates are marked
`absent_from_assay_universe` (not "no signal"); every call traces to a source
row; all zero/no-qualifying calls were measured or explicitly marked
not-measured/absent. Two candidates additionally show `failed_qc` at some
timepoints (effect present but padj `NaN` from DESeq2 independent filtering) —
these are recorded, not treated as signal.

## 6. Conclusion

**Zero qualifying BTLA-specific candidate mRNA is correct.** The old chart's
mRNA counts were produced by (i) pooling the `TCRvsUC` activation contrast into a
BTLA panel and (ii) a seed-inclusive candidate set; neither reflects a
BTLA-specific candidate mRNA response. No bug exists in audit_v2. This is
biologically expected: the ranked regulators were nominated as upstream
explanations of the BTLA transcriptional response, not selected for changes in
their own transcript.
