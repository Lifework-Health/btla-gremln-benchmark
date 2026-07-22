# Commit-safety report — multiomics audit_v2
_Generated: 2026-07-22T12:48:25Z_
## Safe to commit (statuses / counts / hashes)
- [SAFE] multiomics_source_specification.md
- [SAFE] run_manifest_v2.json (hashes + counts + threshold provenance)
- [SAFE] multiomics_evidence_long_v2.csv (statuses; effect_log2fc, adjusted_p, phosphosite residue and direction REDACTED; phospho site id -> REDACTED_site)
- [SAFE] multiomics_candidate_summary_v2.csv (counts + statuses; phospho site-timepoint list -> REDACTED_sites)
- [SAFE] multiomics_model_summaries_v2.csv
- [SAFE] multiomics_coverage_v2.csv
- [SAFE] multiomics_audit_checks_v2.csv
- [SAFE] multiomics_pipeline_audit_v2.md
- [SAFE] figures/*.png / *.svg

## Must NOT be committed
- [RESTRICTED] `local_restricted/multiomics_evidence_long_v2.FULL.csv` and `local_restricted/multiomics_candidate_summary_v2.FULL.csv` — effect sizes, adjusted p, phosphosite residues, per-target direction, site-timepoint lists.
- [RESTRICTED] raw `.xlsx` sources under DATA_ROOT/data/multiomics/ (only SHA256 recorded).

## Pre-commit checklist
- [ ] no file under local_restricted/ staged
- [ ] no .xlsx staged
- [ ] committed CSVs contain no >=3dp decimal effect values
