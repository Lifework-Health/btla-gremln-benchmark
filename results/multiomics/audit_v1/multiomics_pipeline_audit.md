# BTLA multiomics evidence audit (audit_v1)
_Generated: 2026-07-22T09:59:27Z_
**Scope.** Evidence construction and validation only. This audit does NOT modify the manuscript, Results/Conclusions, Table 4, the benchmark verdict, the GREmLN/GENIE3 rankings, the CRISPRi analysis, or the Paperclip tiers. Rankings are read-only inputs; evidence joins cannot alter any rank.
## 1. Headline audit findings
- **Coverage gap:** only **24/43** union candidates were ever scored by the source `08e` pipeline. The remaining **19** are recorded as `source_unavailable_not_evaluated_by_08e_pipeline` — NOT as "no evidence". They were outside the gene set the 08e notebook iterated over and must be re-run against the frozen 43-candidate union before any absence claim.
  - Not evaluated: CYCS, EZH2, FOXP1, GAR1, GTF2A2, GTF3C2, HIRIP3, HMGN3, MLX, NR4A2, PLAGL2, SMAP2, SNAPC4, SRP9, STAT5B, STAU2, TFDP1, ZNF121, ZNF706
- **Independent molecular evidence is rare:** exactly **one** candidate (**MYC**) has any class-A (independent observed molecular) positive item that passes the documented significance and effect thresholds — MYC, via 3 phosphosites.
- **Direction is not asserted:** for magnitude-only molecular layers (RNA, protein, phosphosite, TF/kinase activity) the audit records `positive_unresolved_direction`. We do not assume an increase or decrease is supportive of the BTLA inhibitory phenotype without a documented per-target functional model. Consequently the strict decision-rule field `independent_molecular_support_for_decision_rule` is **True for 0 candidates** (MYC included) pending a directional rule from the data owner.
- **Legacy over-crediting:** 17 candidates the legacy summary labelled `moderate`/`strong` multiomics support are reclassified, because the legacy field treated any non-empty support string (including class-C contextual TF-activity/BIONIC) as support: EGR2, FOXM1, HIVEP3, ID2, MYC, TGIF1, AHR, BCL3, NFIL3, RBPJ, STAT4, CREM, IRF8, JUNB, PRDM1, STAT1, ZBTB32.
## 2. Independence classification (why MYC ≠ the others)
| class | meaning | layers | counts toward decision rule |
|---|---|---|---|
| A | independent observed molecular assay | proteomics, phosphoproteomics, coIP | yes (if supportive) |
| B | orthogonal assay-derived computational | kinase_activity | no |
| C | transcript-/model-derived contextual | transcriptomics, tf_activity, bionic_gnn | no |
| D | curated / unclear lineage | early_synapse_trafficking | no (excluded) |
`tf_activity` and `transcriptomics` are **class C** because they are computed from (or are) the transcript signal that drives candidate ranking; crediting them as independent confirmation would double-count the ranking input.
## 3. Missingness vocabulary (never silently converted to 'no support')
- `not_measured` — entity not present in an assay that was run.
- `measured_no_qualifying_signal` — measured but did not pass thresholds.
- `source_unavailable_not_evaluated_by_08e_pipeline` — the assay/candidate was not processed by the source pipeline at all.
- `not_applicable` — curated hypothesis annotation, not a measurement.
## 4. Directional interpretation audit
Every positive molecular item carries `relation_to_BTLA_hypothesis`. For all magnitude-only layers this is `unresolved` by policy. coIP is `context only` (pulldown membership does not resolve a direct interaction; the CRAPome file is present but was NOT applied in 08e). No item is scored `supportive` without a documented directional model, so the decision rule cannot inflate.
## 5. Layer registry
See `multiomics_layer_registry.csv` (8 layers). Key provenance: RNA thresholds padj<0.05 & |log2FC|≥log2(1.5); protein/phospho padj<0.05 & |logFC|≥log2(1.3); coIP |t|>2 on a pre-filtered significant table; TF/kinase activity decoupleR ULM raw p<0.05 with no effect floor.
## 6. Reconciliation with legacy tables
`legacy_vs_audited_multiomics.csv`: 43 rows. 7 unchanged, 17 reclassified downward, 19 added to audit scope (absent from legacy summary).
## 7. Assertions
All **18/18** audit checks pass:
- [x] union has exactly 43 candidates (43)
- [x] GREmLN contributes 25 (25)
- [x] GENIE3 contributes 25 (25)
- [x] shared == 7 (7)
- [x] summary has each candidate exactly once (43)
- [x] every positive call has source_row_hash (82)
- [x] every positive statistical call passed documented thresholds (79)
- [x] every phosphosite row has a site id (28)
- [x] every coIP row has a bait (0)
- [x] class A only from observed molecular assays ({'phosphoproteomics', 'proteomics'})
- [x] decision field requires >=1 class A supportive item (0)
- [x] no class C/D layer feeds decision field (enforced by construction (class A only))
- [x] not_measured never positive (8)
- [x] opposing not counted as supportive (structural)
- [x] shared candidates single consistent summary row (single-view summary)
- [x] GREmLN view has 25 rows (25)
- [x] GENIE3 view has 25 rows (25)
- [x] ranks sourced only from frozen union (structural)
## 8. Unresolved questions for the data owner
1. Re-run the 08e pipeline over the frozen 43-candidate union so the 19 `source_unavailable` candidates get real statuses.
2. Provide a documented per-layer/per-target directional model (does an increase in a given phosphosite/protein support the BTLA inhibitory phenotype?) so `relation_to_BTLA_hypothesis` can move from `unresolved` to supportive/opposing.
3. Decide whether CRAPome filtering should be applied to coIP before any interaction claim.
4. Confirm MYC phosphosite residues/effect sizes (restricted) before external use.
## 9. Files
- `source_inventory.csv`
- `run_manifest.json`
- `candidate_union.csv`
- `multiomics_layer_registry.csv`
- `multiomics_evidence_long_audited.csv`
- `multiomics_candidate_summary_audited.csv`
- `legacy_vs_audited_multiomics.csv`
- `multiomics_audit_checks.csv`
- `figures/fig_multiomics_audit_pipeline.png`
- `figures/fig_multiomics_candidate_matrix.png`
- `figures/fig_multiomics_model_summary.png`
- `local_restricted/multiomics_evidence_long_audited.FULL.csv` — **git-ignored**; effect sizes + phosphosite residues.
