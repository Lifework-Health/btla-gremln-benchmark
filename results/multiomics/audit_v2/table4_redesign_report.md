# Table 4A/4B publication redesign report

Presentation-only redesign. Evidence values unchanged from audited CSVs
(`table4a_gremln_evidence.csv`, `table4b_genie3_evidence.csv`).

## Layout

- Page size: **14 × 8.5 in landscape** (one table per page; no split required)
- PNG: 300 dpi (~4050 × 2568 px)
- SVG: primary vector output
- Title: **15.5 pt** bold
- Headers: **10.75 pt** bold (two-level; grouped “BTLA experimental evidence”)
- Body: **10.25 pt** (minimum enforced 9.5 pt)
- Notes: **8.5 pt**, two columns
- Max body lines in any cell: **1** (no cell wraps to 2 lines)
- Journal style: dark top/header/bottom rules; no vertical rules; light zebra rows; no pills

## Compact display mapping (presentation only)

| Source CSV language | Display label |
|---|---|
| Detected — BTLA-concordant (x/y) | Concordant (x/y) |
| Detected — BTLA-anti-concordant (x/y) | Anti-concordant (x/y) |
| No detected response (0/y) | No response (0/y) |
| Not represented in screen | Not screened |
| Failed on-target quality control | Failed QC |
| TCR activation ↑: … | TCR ↑ … |
| No qualifying transcriptional response | No qualifying change |
| Absent from transcriptomic assay universe | Not measured |
| Phosphosite association: … | Phospho: … |
| Protein association: … | Protein: … |
| No qualifying association in measured layers | No qualifying signal |
| Not measured in protein or phosphoproteomic assays | Not measured |

## Reconciliation (preserved)

- GREmLN: CRISPRi 3; lit strong/moderate 7; mRNA 2; independent 4
- GENIE3: CRISPRi 6; lit strong/moderate 9; mRNA 5; independent 3
- Shared independent: 0

## Assertions

All typography, wrap, order, and scientific assertions passed. Generation would
fail rather than reduce body text below 9.5 pt.

## Outputs

- `results/multiomics/audit_v2/table4a_gremln_evidence.svg/.png`
- `results/multiomics/audit_v2/table4b_genie3_evidence.svg/.png`
- Source CSVs retained unchanged beside these renders.
