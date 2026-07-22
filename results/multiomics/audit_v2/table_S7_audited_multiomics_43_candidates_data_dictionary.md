# Table S7 data dictionary

Source: final audit_v2 outputs (multiomics_candidate_summary_v2.csv, mrna_contrast_layers_summary_v2.csv, multiomics_evidence_long_v2.csv) and the frozen union results/publication_data/top25_union_primary.csv. No values were rerun or reinterpreted.

| Column | Description |
| --- | --- |
| Candidate regulator | Canonical gene symbol. |
| Candidate group | GREmLN specific, Shared, or GENIE3 specific. |
| GREmLN rank | Canonical GREmLN ordinal (dense) rank if in the GREmLN top 25, else em dash. |
| GENIE3 rank | Canonical GENIE3 ordinal (dense) rank if in the GENIE3 top 25, else em dash. |
| TCR activation mRNA | Contextual. TCR vs unstimulated control. Increased/Decreased with qualifying timepoint(s), or a controlled status (Measured, no qualifying change / Absent from assay universe / Ambiguous mapping / Failed quality control). |
| Incremental BTLA mRNA | Contextual. BTLA+TCR vs TCR. Same vocabulary; never collapsed with the TCR activation contrast. |
| Protein abundance | Independent. BTLA+TCR vs TCR. Qualifying association with timepoint(s) (2 min, 20 min, 4 h, 24 h) or a controlled status. |
| Phosphosite abundance | Independent. BTLA+TCR vs TCR. Qualifying association with timepoint(s) and the number of distinct qualifying sites; residues and site identifiers are not exposed and sites are not averaged. |
| Independent molecular association | Yes — protein / phosphosite / protein and phosphosite; No qualifying association; Not assessable because not measured. True only for a qualifying directly measured protein or phosphosite change in BTLA+TCR vs TCR. |
| Interpretation note | Concise contextualisation (e.g. Activation-associated mRNA only; Independent protein/phosphosite association; No qualifying signal in measured layers; Sparse assay coverage). |

Thresholds: transcriptomics adjusted p<0.05 and |log2FC|>log2(1.5); protein and phosphosite adjusted p<0.05 and |log2FC|>log2(1.3). Harmonised sensitivity threshold |log2FC|>=0.585 retained the same seven independent associations. Restricted effect sizes, adjusted p values and phosphosite residues/site identifiers remain only in git-ignored local_restricted/ files.
