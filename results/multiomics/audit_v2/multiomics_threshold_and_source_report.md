# Threshold sensitivity and source-version reconciliation (sections 4-5)

_Evidence checks only. Manuscript, Table 4 and verdict unchanged._

## Section 4 — protein / phosphosite threshold sensitivity

**Primary (source-native) thresholds confirmed:**
- transcriptomics: `padj < 0.05` and `|log2FC| > log2(1.5) ≈ 0.585`
- protein and phosphoproteomics: `padj < 0.05` and `|log2FC| > log2(1.3) ≈ 0.379`

The seven reported independent molecular associations use the **protein/phospho
source-specific** threshold (0.379): MYC, SMAP2, STAT5B (protein: SMAP2, STAT5B;
phospho: MYC, GTF3C2, HIRIP3, STAU2, TFDP1).

**Harmonised sensitivity re-run** (stricter, `|log2FC| >= 0.585` for protein and
phospho):

| | primary (log2(1.3)) | harmonised (0.585) |
|---|---|---|
| candidates with independent association | 7 | 7 |
| retained | — | GTF3C2, HIRIP3, MYC, SMAP2, STAT5B, STAU2, TFDP1 |
| lost | — | **none** |
| gained | — | none |
| GREmLN-specific | 4 (GTF3C2, HIRIP3, STAU2, TFDP1) | 4 |
| shared | 0 | 0 |
| GENIE3-specific | 3 (MYC, SMAP2, STAT5B) | 3 |

**Interpretation:** all seven associations are robust — every qualifying
protein/phosphosite item already exceeds the stricter 0.585 cut, so no candidate
changes. Per-candidate layer/timepoint detail is in
`protein_phospho_threshold_sensitivity.csv` (site residues not shown). The
harmonised analysis is reported for sensitivity only and does **not** replace the
source-specific primary analysis.

## Section 5 — source-version reconciliation (gmorgan vs Viet)

Both data roots were compared for the three differential-expression workbooks:

| layer | gmorgan file | Viet file | SHA256 | schema | coverage (43) | qualifying diff |
|---|---|---|---|---|---|---|
| transcriptomics | `Transcriptomics/@Transcriptomics/Tables/DEG_wide_statistics.xlsx` | `DEG_wide_statistics.xlsx` | **identical** (`adab2a38…`) | identical | 41/41 both (universe 16725) | n/a |
| proteomics | `Proteomics/@Proteomics/Tables/WP_DEP_wide_statistics.xlsx` | `WP_DEP_wide_statistics.xlsx` | **identical** (`47f6e633…`) | identical | 16/16 both | none (both: SMAP2, STAT5B) |
| phosphoproteomics | `Phosphoproteomics/@Phosphoproteomics/@Tables/DEPPs_All_statistics.xlsx` | `DEPPs_All_statistics.xlsx` | **identical** (`47d36cd9…`) | identical | 5/5 both | none (both: GTF3C2, HIRIP3, MYC, STAU2, TFDP1) |

All three raw workbooks are **byte-identical** across the two roots (same SHA256,
sheets, columns, contrasts, timepoints, coverage, qualifying calls). The Viet
root stores them in a flat directory; the gmorgan root uses the nested
`@`-prefixed layout, but the file contents are the same release.

Methods and thresholds are shared (single `Documents/Documents/Methods.docx`:
limma for mass-spec, DESeq2 for transcriptomics; protein/phospho FC>1.3, RNA
FC>1.5, padj<0.05). The earlier "70 vs 88 gene" difference concerned downstream
summary CSVs / candidate lists, **not** these raw DE sources.

**Decision:** the audit uses the **gmorgan** copies already frozen in
`run_manifest_v2.json` (`raw_source_sha256`). Because the Viet copies are
byte-identical, no source switch is required and none was made. The selected
release and hashes remain frozen in the manifest; `source_version_freeze.json`
records the side-by-side comparison and the selection rationale.
