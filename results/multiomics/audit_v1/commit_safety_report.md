# Commit-safety report — multiomics audit_v1
_Generated: 2026-07-22T09:59:27Z_
## Safe to commit (statuses / hashes / aggregate counts only)
- [SAFE] source_inventory.csv (redacted relative paths + SHA256, no effect sizes)
- [SAFE] run_manifest.json (hashes + counts)
- [SAFE] candidate_union.csv (membership + ranks, already published)
- [SAFE] multiomics_layer_registry.csv (schema/provenance)
- [SAFE] multiomics_evidence_long_audited.csv (statuses only; effect_value, adjusted_p, phosphosite id/residue and per-target direction REDACTED)
- [SAFE] multiomics_candidate_summary_audited.csv (class counts + statuses)
- [SAFE] legacy_vs_audited_multiomics.csv (status-level reconciliation)
- [SAFE] multiomics_audit_checks.csv
- [SAFE] multiomics_pipeline_audit.md
- [SAFE] figures/*.png / *.svg (status colours only, no numeric values)

## Must NOT be committed (publication gate)
- [RESTRICTED] `local_restricted/multiomics_evidence_long_audited.FULL.csv` — log2FC/logFC/t-statistics, adjusted p-values, phosphosite residues, per-target direction. Kept under `local_restricted/` and added to `.gitignore`.
- [RESTRICTED] all raw Excel sources under `DATA_ROOT/data/multiomics/` (never staged; only SHA256 recorded).

## Redaction method
- Absolute paths replaced by paths relative to `BTLA_BENCH_DATA_ROOT`.
- Effect sizes and phosphosite residues dropped from committed tables and replaced with `REDACTED_publication_gate`.
- Per-target direction (sign of a significant change) removed from committed tables; only the direction-neutral `positive_unresolved_direction` status remains.

## Pre-commit checklist
- [ ] `git status` shows no file under `local_restricted/`.
- [ ] no `.xlsx` staged.
- [ ] grep committed CSVs for numeric effect columns returns nothing.
