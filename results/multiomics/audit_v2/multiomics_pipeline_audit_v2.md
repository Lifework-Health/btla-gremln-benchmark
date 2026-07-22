# BTLA multiomics evidence audit v2 — pipeline audit report
_Generated: 2026-07-22T12:48:25Z. Supersedes audit_v1._
**Scope.** Evidence construction only; manuscript, Table 4, verdict, rankings, CRISPRi and Paperclip tiers untouched. See `multiomics_source_specification.md` for the approved spec.
## 1. Headline
- Rebuilt directly from the raw assay tables (not the 08e derivative), so coverage is no longer limited by 08e's old gene gate.
- **7 of 43** candidates show an independent molecular association with the BTLA contrast (qualifying protein or phosphosite item): MYC, SMAP2, STAT5B, GTF3C2, HIRIP3, STAU2, TFDP1.
  - Protein-qualifying: SMAP2 (4h, time-matched primary), STAT5B (20min, acute).
  - Phospho-qualifying: MYC, GTF3C2, HIRIP3, STAU2, TFDP1 (site/timepoint detail restricted).
- **No shared candidate** has an independent association; all 7 are model-specific (GREmLN-specific: GTF3C2, HIRIP3, STAU2, TFDP1; GENIE3-specific: MYC, SMAP2, STAT5B).
- **mRNA context change qualifying = 0** across all 43 in the BTLA contrast. Consistent with the clarification that these are upstream regulators whose own transcript need not move (lack of mRNA change is NOT evidence against a candidate).
- Direction is left **unresolved** for every independent item (rankings unsigned; no site-specific functional model). Association is established; supportive/opposing is not claimed.
## 2. Thresholds (source-native, recovered from Methods.docx)
- Transcriptomics: padj < 0.05 AND |log2FC| > log2(1.5) = 0.585
- Proteomics & phosphoproteomics: padj < 0.05 AND |log2FC| > log2(1.3) = 0.379
- Harmonised sensitivity flag (not primary): |log2FC| >= 0.585.
- 08e's non-source `strong` sub-tier dropped: items are qualifying or not.
## 3. Coverage (raw, all 43)
| layer | role | in universe | absent | qualifying |
|---|---|---|---|---|
| transcriptomics | contextual | 41 | 2 | 0 |
| proteomics | independent | 16 | 27 | 2 |
| phosphoproteomics | independent | 5 | 38 | 5 |
| tf_activity | contextual | 28 | 15 | 10 |
| coIP | excluded | 43 | 0 | 0 |
## 4. Model summaries
| view | n | independent assoc. | protein qual | phospho qual | mRNA context | TF-act context |
|---|---|---|---|---|---|---|
| GREmLN_top25 | 25 | 4 | 0 | 4 | 0 | 5 |
| GENIE3_top25 | 25 | 3 | 2 | 1 | 0 | 9 |
| shared | 7 | 0 | 0 | 0 | 0 | 4 |
| GREmLN_specific | 18 | 4 | 0 | 4 | 0 | 1 |
| GENIE3_specific | 18 | 3 | 2 | 1 | 0 | 5 |
| union_43 | 43 | 7 | 2 | 5 | 0 | 10 |
## 5. Three separated concepts
1. measured molecular change · 2. independent molecular association (BTLA contrast, direction-agnostic) · 3. directionally interpretable support (**unresolved** for all). A qualifying protein/phosphosite change establishes (1) and (2) without implying (3).
## 6. coIP
Excluded from corroboration: the available file is an intensity matrix with no statistical BTLA contrast, and the BTLA-IP-vs-IgG differential named in the Methods is absent from this tree. Only GAR1 (of 43) is detected. Presence is recorded, never scored.
## 7. Assertions
All 16/16 checks pass:
- [x] union == 43 (43.0)
- [x] summary rows == 43 (43.0)
- [x] GREmLN view == 25 (25.0)
- [x] GENIE3 view == 25 (25.0)
- [x] shared == 7 (7.0)
- [x] GREmLN_specific + shared == 25 (nan)
- [x] GENIE3_specific + shared == 25 (nan)
- [x] mRNA rows are all contextual role (nan)
- [x] tf_activity rows are all contextual role (nan)
- [x] coIP excluded from corroboration (nan)
- [x] every independent association backed by a qualifying protein/phospho item (nan)
- [x] not_measured never qualifying (nan)
- [x] absent_from_assay_universe never qualifying (nan)
- [x] RNA threshold = log2(1.5) (0.5849625007211562)
- [x] protein/phospho threshold = log2(1.3) (0.3785116232537298)
- [x] independent count reconciles (union) (nan)
## 8. Unresolved questions for the data owner
1. Provide the BTLA-IP-vs-IgG coIP differential file to enable an interaction layer.
2. Provide a site-specific functional model so phosphosite direction can move from `unresolved` to supportive/opposing.
3. Confirm the linear-fold-change reading of the Methods thresholds.
4. Confirm the gmorgan multiomics vintage is the intended source (vs Viet's newer copy).
## 9. Files
- `multiomics_source_specification.md`
- `run_manifest_v2.json`
- `multiomics_evidence_long_v2.csv (redacted)`
- `multiomics_candidate_summary_v2.csv`
- `multiomics_model_summaries_v2.csv`
- `multiomics_coverage_v2.csv`
- `multiomics_audit_checks_v2.csv`
- `figures/fig_multiomics_v2_candidate_matrix.png`
- `figures/fig_multiomics_v2_model_summary.png`
- `figures/fig_multiomics_v2_coverage.png`
- `local_restricted/*.FULL.csv` — git-ignored (effect sizes, adjusted p, phosphosite residues, per-target direction).
