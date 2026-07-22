# Candidate mRNA across the full experimental sequence (audit_v2 refinement)

Candidate mRNA is **contextual** evidence and does **not** enter the independent
molecular corroboration count. This refinement describes transcriptional
responsiveness of the 43 candidates across three prespecified contrasts, kept
separate and never collapsed into one strongest result.

## Contrasts

| # | Key | Biological contrast | Source | Status |
|---|-----|---------------------|--------|--------|
| C1 | `TCR_activation` | TCR vs unstimulated control | `TCRvsUC` (DEG workbook) | available |
| C2 | `Incremental_BTLA` | BTLA+TCR vs TCR | `BTLAvsTCR` (DEG workbook) | available |
| C3 | `Combined_BTLA_TCR_vs_UC` | BTLA+TCR vs unstimulated control | — | `source_unavailable_pending_model` |

Source-native transcriptomics threshold (Methods.docx):
**adjusted p < 0.05 AND |log2 fold change| > log2(1.5)**.

## Does contrast 3 exist? (confirmation step)

**No.** Both source workbooks —
`DEG_statistics.xlsx` (DESeq2 long-form, columns `baseMean, log2FoldChange,
lfcSE, stat, pvalue, padj`) and `DEG_wide_statistics.xlsx` — contain only
`BTLAvsTCR` and `TCRvsUC`, each at 1h/4h/24h. The Methods state only pairwise
per-timepoint designs were fitted ("Experimental designs were generated for each
stimulation condition combination for each time point, e.g. TCS stimulation vs
TCS+BTLA stimulation at 1 hour"). A combined BTLA+TCR-vs-UC Wald contrast was
never produced.

**Can it be regenerated?** Not currently. Regenerating a proper contrast requires
the fitted DESeq2 `dds` object (dispersion estimates) or the raw counts + coldata
/ design matrix. An exhaustive search of both data roots (gmorgan and Viet) found
no `.rds`/`.RData`, no count matrix, no coldata, and no DE R script; `pydeseq2`
is not installed. Contrast 3 is therefore represented by the explicit,
non-fabricated state `source_unavailable_pending_model` and is **never**
synthesised by adding fold changes or combining p-values from C1 and C2.

> **Data request:** to complete C3, provide the DESeq2 `dds` object (or the raw
> count matrix + sample/design table) for the RNA-seq timecourse so the
> BTLA+TCR-vs-unstimulated contrast can be fitted with the original design.

## Counts (of 43 candidates)

| Metric | n |
|--------|---|
| Qualify in C1 (TCR activation) | **5** |
| Qualify in C2 (incremental BTLA) | **0** |
| Qualify in C3 (combined) | pending (not measurable) |
| Qualify in more than one available contrast | **0** |
| Qualify in any available context | **5** |

C1 qualifiers (all up-regulated): **EGR2** (1h;4h;24h), **HIVEP3** (4h;24h),
**ID2** (4h), **JUNB** (1h;24h), **PRDM1** (4h).

By group (n; C1; C2; any): GREmLN-specific (18; 0; 0; 0); shared (7; 2; 0; 2);
GENIE3-specific (18; 3; 0; 3).

mRNA measured denominator (in DEG universe): GREmLN-specific 16/18,
shared 7/7, GENIE3-specific 18/18.

## Preserved per-candidate fields

`mrna_contrast_layers_summary_v2.csv` and `mrna_contrast_layers_long_v2.csv` keep,
for each candidate and contrast: qualifying contrast, qualifying timepoint(s),
direction (sign), and measured / no-signal / not-measured / absent / failed-QC /
pending status. The descriptive flag `transcriptionally_responsive_any_context`
is always accompanied by `responsive_contexts` (contrast@timepoint:direction).

## Redaction

Committed files carry statuses, qualifying contrast/timepoint and direction sign
only. Numeric log2FC / adjusted p are written only to
`local_restricted/mrna_contrast_layers_values.RESTRICTED.csv` (git-ignored).

_No manuscript or benchmark-verdict change was made._
