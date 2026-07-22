# BTLA multiomics — source specification (approved before implementation)

_Evidence construction only. Manuscript, Table 4, benchmark verdict, model
rankings, CRISPRi and Paperclip tiers are untouched. Rankings are read-only
inputs; multiomics evidence cannot alter any rank._

This document records the ground truth recovered from the raw source tables and
the study Methods, and the rules approved for `audit_v2`. It supersedes
`audit_v1`, which was built from the 08e long-form derivative and inherited that
notebook's gene gate and coIP handling.

## 1. Source files and exact columns (gmorgan vintage, consumed by the benchmark)

All under `DATA_ROOT/data/multiomics/`. SHA256 in `run_manifest_v2.json`.

| Layer | File | Key columns |
|---|---|---|
| Transcriptomics | `Transcriptomics/@Transcriptomics/Tables/DEG_wide_statistics.xlsx` | `gene`; `BTLAvsTCR_{1h,4h,24h}_log2FoldChange`; `BTLAvsTCR_{...}_padj`; `TCRvsUC_{...}` |
| Proteomics | `Proteomics/@Proteomics/Tables/WP_DEP_wide_statistics.xlsx` | `protein` (UniProt_Gene), `Gene`; `Btla_{2mn,20mn,4h,24h}_logFC`; `Btla_{...}_adj.P.Val`; `TCRvsUC_{...}` |
| Phosphoproteomics | `Phosphoproteomics/@Phosphoproteomics/@Tables/DEPPs_All_statistics.xlsx` | sheets `Btla_{2mn,20mn,4h,24h}`, `TCR_{...}`; per-row `rn`, `logFC`, `t`, `P.Value`, `adj.P.Val`, `Gene_STY_canonical` (gene+residue), `Gene` |
| coIP | `coIP/@coIP/@Tables/@BTLA_coIP_signif_wide.xlsx` | `protein`, `Gene`; intensity columns `TcrBtla_*`, `Tcr_*`, `uc_*` (no stats) |
| TF activity (contextual) | `Transcriptomics/@Transcriptomics/Tables/TFactivitiesALL.xlsx` | `source`; `{BTLAvsTCR,TCRvsUC,BTLA_vsUC}_{1h,4h,24h}_{score,p_value}` |

## 2. Contrast meanings
- **`BTLAvsTCR` / `Btla_*`** — anti-BTLA crosslink + TCR **vs** TCR alone: the BTLA-specific effect. **Primary contrast used for all evidence.**
- **`TCRvsUC` / `TCR_*`** — TCR **vs** unstimulated: generic activation; contextual only, not used for BTLA support.
- **`BTLA_vsUC`** (TF activity only) — BTLA+TCR vs unstimulated.

## 3. Fold-change scale
log2 throughout (MaxLFQ intensities log2-transformed; limma `logFC` and DESeq2
`log2FoldChange` are both log2). Timepoints: RNA = 1h/4h/24h; protein & phospho
= 2min/20min/4h/24h. **4h** is the time-matched primary view; acute (2/20min)
exists only for protein and phospho.

## 4. Thresholds — recovered from `Documents/Documents/Methods.docx`
> "The threshold to identify significantly differentially regulated proteins and
> phosphoproteins was adjusted to p-value < 0.05 and absolute-fold-change > 1.3.
> For transcriptomics, the threshold was set to adjusted p-value 0.05 and
> absolute-fold-change > 1.5." (Differential analysis: limma for mass-spec,
> DESeq2 for transcriptomics.)

Applied source-native rules (primary):
- **Transcriptomics:** `adj.p < 0.05` AND `|log2FC| > log2(1.5) = 0.585`
- **Proteomics & phosphoproteomics:** `adj.p < 0.05` AND `|log2FC| > log2(1.3) = 0.378`

Notes:
- These match the values 08e used, so the earlier provisional thresholds were in
  fact source-correct for the qualifying cut. **However** 08e's extra `strong`
  sub-tier (`adj.p < 0.01 & |FC| >= 1.5x`) is **not** in the Methods and has been
  **dropped**; each item is simply qualifying or not.
- "absolute fold-change > 1.3/1.5" is read as **linear** FC (log2 0.378 / 0.585),
  the standard reading. The stricter alternative (`|log2FC| > 1.3/1.5`) was not used.
- A harmonised sensitivity flag (`|log2FC| >= 0.585` for all layers) is also
  recorded per item but is **not** the primary rule (it is stricter than source
  for protein/phospho).

## 5. Evidence roles (approved)
- **mRNA — contextual only.** We report whether the regulator's own transcript
  changes at 1h/4h/24h in the BTLA contrast; it is never independent
  corroboration. (Per clarification: a candidate is an upstream regulator; lack
  of its own mRNA change, especially at 4h, is expected and is not evidence
  against it.)
- **Protein abundance — independent observed molecular evidence;** timepoints
  analysed separately; 4h primary; acute and late retained.
- **Phosphosites — independent observed molecular evidence;** every site and
  timepoint separate; acute (2/20min) vs sustained (4/24h) distinguished; sites
  never averaged; phosphorylation direction is **not** read as activation.
- **TF activity — contextual** (decoupleR ULM derived from the transcriptome).
- **coIP — EXCLUDED from corroboration.** The available file is an intensity
  matrix (`TcrBtla`/`Tcr`/`uc`, values ~3–16 log2) with no statistical BTLA
  contrast; the proper **BTLA-IP-vs-IgG** differential named in the Methods is
  **not present** in this data tree. Presence is recorded but never scored.

## 6. Three separated concepts
1. **measured molecular change** — a qualifying protein/phospho item exists.
2. **independent molecular association with the BTLA contrast** — ≥1 qualifying
   protein or phosphosite item in `BTLAvsTCR`; **direction-agnostic**.
3. **directionally interpretable mechanistic support** — left **`unresolved`**
   because model rankings are unsigned and no site-specific functional model was
   supplied. A qualifying item establishes association even when its direction
   cannot be called supportive or opposing.

## 7. Six-state status vocabulary (per candidate × assay)
`measured_qualifying` · `measured_no_qualifying_signal` · `not_measured` ·
`absent_from_assay_universe` · `ambiguous_mapping` · `failed_qc`
(plus `excluded_no_statistical_contrast` for coIP). **`not_measured` and
`absent_from_assay_universe` are never treated as "no evidence".**

## 8. Coverage (raw, all 43 — bypassing the old 08e gene gate)
| Layer | in universe | absent | with qualifying signal |
|---|---|---|---|
| Transcriptomics (context) | 41 | 2 (AHR, MSC) | 0 |
| Proteomics (independent) | 16 | 27 | 2 |
| Phosphoproteomics (independent) | 5 | 38 | 5 |
| TF activity (context) | 28 | 15 | 10 |
| coIP (excluded) | 43 detected-or-not | — | 0 (excluded) |

Phosphosites per present candidate: HIRIP3 56, GTF3C2 32, MYC 28, STAU2 20,
TFDP1 8.

## 9. Mapping / QC issues
- coIP is not a statistical contrast (see §5); only GAR1 (of 43) is even detected.
- Proteomics `Gene` maps cleanly for the 43 (no multi-protein-per-gene collisions);
  UniProt IDs carry isoform suffixes (e.g. `-2`). Any future collision →
  `ambiguous_mapping`.
- Phospho maps via `Gene` + residue (`Gene_STY_canonical`, format `GENE_<residue>`,
  residues restricted); multiple sites/gene kept separate.
- mRNA was labelled `independent` in 08e; corrected here to contextual.

## 10. Correction vs audit_v1
audit_v1 reported only MYC with independent signal and 19 candidates as
"source unavailable". That was because 08e iterated only over
`btla_tf_candidate_union.csv` (an older list). Reading raw sources directly,
**7** candidates have an independent molecular association
(protein: SMAP2, STAT5B; phospho: MYC, GTF3C2, HIRIP3, STAU2, TFDP1); the
"source unavailable" class is replaced by proper per-assay `absent_from_assay_universe`
vs `not_measured` statuses.
