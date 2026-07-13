<!--
COLLABORATIVE DRAFT — BTLA / GREmLN vs GENIE3 benchmark report
Target length: 2000–2500 words.

Angle: transparency-first benchmark that uses orthogonal CD4 CRISPRi Perturb-seq
perturbational validation to compare a graph-aware foundation model (GREmLN) against
classical GRN inference (GENIE3) for BTLA-response transcription-factor prioritisation.

Anchored on the definitive comparison in notebook 10
(GREmLN built on the raw masked CD4 GENIE3 50k-edge graph as its prior vs that same
masked-raw GENIE3 graph). Validators = (1) orthogonal CD4 CRISPRi perturbational
validation, (2) Paperclip literature, (3) BTLA multi-omics. Held-out recovery is
SUPPLEMENTARY only (pending a leakage/asymmetry audit). PRIMARY nomination comparison =
seed-EXCLUDED top-25; seed-INCLUSIVE rankings are a biological-context analysis only.

Data note: BTLA transcriptomic/multi-omic data-publication permission is unresolved; do
NOT state that code/results are open until the separate publication repository is live and
data clearance is confirmed.

How we edit together:
- Leave TODO/❓ markers inline where you want input or a number verified.
- Export later via: pandoc btla_gremln_benchmark.md -o paper.docx
-->

# Benchmarking GREmLN against GENIE3 for BTLA-response transcription-factor prioritisation using orthogonal CRISPRi validation

<!-- TODO: author list + affiliations -->

## Abstract

Identifying the transcription factors (TFs) that shape co-inhibitory checkpoint signaling is a prerequisite for engineering T-cell states, yet regulator-inference methods are rarely compared against independent perturbation data. We built one definitive graph-aware single-cell foundation model, GREmLN, using a masked CD4 GENIE3 co-expression graph as its regulatory prior, and benchmarked it against that same masked-raw GENIE3 graph for BTLA-response TF prioritisation in human CD4⁺ T cells. Candidate TFs were ranked by graph proximity to a BTLA-specific differential-expression (DEG) panel (BTLA+TCR vs TCR, 4 h, 249 genes) using CSLS embedding neighbourhoods (GREmLN) or summed outgoing TF→BTLA-DEG edge weight (GENIE3), over a common 334-TF universe scorable by both methods. The **primary regulator-nomination comparison excludes the BTLA panel genes themselves** (seed-excluded top-25); a seed-inclusive ranking is reported separately as biological context. Candidates were then examined with three orthogonal, non-predictive evidence layers: CD4 CRISPRi Perturb-seq perturbational validation, Paperclip literature review, and BTLA multi-omics. Over the common universe the two rankings were moderately correlated (Spearman ρ ≈ 0.46), sharing 7 of 25 primary TFs. CRISPRi functional target-response support was modest and overlapping (mean confirmed-target fraction 0.049 for GREmLN [95% CI 0.00–0.13] and 0.054 for GENIE3 [0.01–0.11]); **neither method met the superiority decision rule specified for the definitive comparison.** Functional target-response support concentrated in GENIE3-side activation regulators (EGR2, BHLHE40, JUNB) and the seed TFs REL and TBX21, whereas GREmLN-specific candidates (AHR, NFIL3, PRDM1, RBPJ, STAT4) were supported mainly by literature and derived multi-omics and largely **lacked orthogonal CRISPRi functional support** (AHR showing contradictory-direction target movement; NFIL3/PRDM1 absent from the screen). Because both scoring schemes are unsigned, target movement after knockdown does not by itself establish regulatory direction. We conclude that a graph foundation model is competitive with, but not shown superior to, classical GRN inference for this task, and report a prioritised, evidence-tiered BTLA-response TF shortlist rather than a single winner.

## Introduction

BTLA (B- and T-lymphocyte attenuator) is a co-inhibitory receptor that restrains T-cell activation and is an emerging target for immunotherapy. Understanding *how* BTLA engagement reshapes the T-cell transcriptome — and which transcription factors mediate that response — would inform both checkpoint modulation and cell-engineering strategies. The upstream regulators of the BTLA response remain poorly mapped relative to canonical TCR co-stimulation.

Two paradigms dominate regulator inference from single-cell data. The first, exemplified by GENIE3 and the SCENIC ecosystem, learns a gene regulatory network (GRN) from co-expression using tree-based feature importance to weight TF→target edges. It is interpretable and well validated, but it infers association rather than regulation and is sensitive to dataset composition. The second, newer paradigm is the graph-aware single-cell foundation model: GREmLN (Columbia / CZ Biohub NY), trained on ~11M CELLxGENE single-cell profiles, folds GRN topology into the attention mechanism of a transformer via graph signal processing (a diffusion kernel on the graph Laplacian), producing gene embeddings intended to transfer regulatory structure across contexts. Such models promise richer, context-general representations, but their practical value over classical inference is largely untested — and rarely examined against independent perturbation data.

That gap defines this work. We ask a deliberately falsifiable question: does GREmLN prioritise BTLA-response TF regulators more accurately than GENIE3? We answer it not with a single leaderboard number but with a transparency-first, multi-layered benchmark whose keystone is a genome-scale CD4 CRISPRi Perturb-seq screen used as orthogonal perturbational validation — direct TF knockdowns that let us ask whether *predicted* targets actually move when a regulator is silenced. The payoff is twofold: an honest read on where foundation-model embeddings add value, and a prioritised, orthogonally supported shortlist of BTLA-response regulators.

## Methods

**Datasets.** (i) The CZI genome-scale T-cell Perturb-seq resource provided CD4 cells for both GRN construction (non-targeting-control cells) and orthogonal perturbational validation (Rest / Stim8hr / Stim48hr arms). (ii) BTLA differential-expression panels were derived from bulk RNA-seq of anti-BTLA cross-linked primary T cells. The primary panel is **BTLA+TCR vs TCR at 4 h** (249 DEGs; padj < 0.05, |log₂FC| ≥ 0.585), isolating the incremental BTLA effect on a TCR-activated background, with direction retained (BTLA-up / BTLA-down). BTLA source transcriptomic and multi-omic data are subject to publication-permission review (see Code and data availability).

**The definitive model comparison.** The head-to-head compares two models built on the *same* regulatory graph. Stated precisely:

* GENIE3 was inferred from **donor-adjusted (ComBat) CD4 expression** (50k NTC cells, 4k highly variable genes, CPM/log1p), persisted as a 50,000-edge masked graph (334 regulators, 3,557 targets).
* GREmLN used the **corresponding pre-ComBat log1p-CPM expression** of the identical 50k cells, because its tokenizer requires the rank/zero structure of raw counts (ComBat output contains negative values and few true zeros).
* GREmLN used the **masked GENIE3 graph as its regulatory prior**.
* The comparison therefore estimates **the value of the GREmLN model and embedding procedure when supplied with the GENIE3 prior; it is not a comparison of two independent methods on identical numerical expression inputs.**

GREmLN (checkpoint `model.ckpt`, MD5 6e57…; used zero-shot, not fine-tuned) embedded 3,586 genes (342 TFs) in 512 dimensions. Full provenance — checkpoint/expression/graph/TF-list hashes, versions, parameters, seed and device — is recorded in a model-registry table.

**TF ranking and scoring.** GREmLN TF scores use CSLS (Cross-domain Similarity Local Scaling, k = 10; the authors' default), a mutual-neighbour metric correcting the hubness of cosine similarity: each TF is scored by the number of BTLA-DEG seeds among its top-100 embedding neighbours, with dense tie-aware ranking. GENIE3 TFs are scored by **summed outgoing TF→BTLA-DEG edge weight** (from the TF/regulator node's perspective, the summed weight of its edges into DEG targets). All comparisons use the **common 334-TF universe** scorable by both methods (GREmLN-embedded TFs ∩ masked-GENIE3 regulators); to keep the pools symmetric, common-universe TFs with no outgoing GENIE3 edge into the panel are assigned a score of zero and tie at the bottom of the GENIE3 ranking, exactly as GREmLN scores every TF (including zero seed-neighbour counts). The **primary regulator-nomination comparison is seed-excluded** (BTLA panel genes removed, so the ranked TFs are candidate regulators, not response genes); a **seed-inclusive ranking is reported separately as biological context** and its panel genes are flagged, never described as model discoveries.

**Orthogonal CD4 CRISPRi validation.** For each model's top-25 TFs we tested whether its predicted BTLA-DEG targets (GENIE3: masked-graph targets ∩ seeds; GREmLN: top-100 CSLS neighbours ∩ seeds) are themselves differentially expressed when that TF is knocked down in the CD4 CRISPRi **Stim8hr** arm, reporting the functional target-response fraction (confirmed / predicted) and a hypergeometric enrichment p-value; Stim48hr and the Stim8hr∪Stim48hr union are sensitivity arms. The **custom target-response threshold** is adjusted p < 0.05 on the knockdown DE test; the **dataset-provided on-target threshold** is the CZI Perturb-seq `ontarget_significant` flag. TFs failing on-target QC or absent from the screen are reported as such, never as zero support. A signed-agreement statistic records whether confirmed targets move consistently with (or opposite to) the BTLA programme. Because both scoring schemes are **unsigned**, target movement after knockdown does not by itself establish regulatory direction. A superiority decision rule was specified for the definitive comparison: a model is declared superior only if it wins CRISPRi by a meaningful margin (mean functional target-response fraction higher by ≥ 0.02 and at least as many `usable_supportive` TFs) *and* carries at least as much corroborating Paperclip + multi-omics support.

**Corroborating evidence (annotation, not prediction).** Literature support is compiled with the Paperclip full-text retrieval tool; BTLA multi-omics evidence is preserved as layer-specific long-form annotations. Neither is used as a predictive input. Multi-omics layers are reported in two explicitly separated groups (below).

## Results

**A moderately convergent primary comparison.** Over the common 334-TF universe, GREmLN's CSLS ranking and GENIE3's edge-weight ranking were moderately correlated (Spearman ρ = 0.50, seed-excluded primary; ρ = 0.52 seed-inclusive), indicating the embeddings reorganise rather than overturn the co-expression prioritisation. Of the 25 seed-excluded (primary) TFs per method, **7 were shared** (CREM, IRF8, JUNB, PLAGL2, PRDM1, STAT1, ZBTB32) and 18 were unique to each (Figures 1–2). GREmLN uniquely elevated AHR, NFIL3, PRDM1 (also shared), RBPJ, STAT4 and EZH2; GENIE3 uniquely elevated activation-associated EGR2, BHLHE40, HIVEP3, MYC, NR4A2 and ZEB2.

<!-- Figures are regenerated by notebooks/03_benchmark_btla_tf_prioritisation.ipynb into
     results/figures/ and are gitignored until BTLA data-publication clearance. -->
![Figure 1](../results/figures/fig1_top25_seed_excluded.png)

**Figure 1 (primary). Seed-excluded top-25 candidate TFs, GREmLN (CSLS, masked prior) vs GENIE3 (masked raw).** Bars are TF scores over the common universe with BTLA panel genes removed, so entries are candidate regulators. GENIE3 uses summed outgoing TF→BTLA-DEG edge weight; GREmLN uses the CSLS seed-neighbour score.

![Figure 2](../results/figures/fig2_top25_overlap_seed_excluded.png)

**Figure 2 (primary). Seed-excluded top-25 overlap.** Seven candidate TFs are shared; eighteen are model-specific.

**Orthogonal CRISPRi perturbational validation.** The keystone analysis asks whether predicted targets respond to TF knockdown. On the seed-excluded primary lists (Figure 3, Table 1), functional target-response support was modest and overlapping: mean confirmed-target fraction **0.049 for GREmLN** (median 0.00; bootstrap 95% CI 0.00–0.13; 21/25 present in screen, 20 passing on-target QC) and **0.054 for GENIE3** (median 0.00; 95% CI 0.01–0.11; 22/25 present, 19 passing QC). The 0.005 difference is well below the 0.02 decision margin. GREmLN had 1 `usable_supportive`, 18 `usable_unsupported`, 1 `usable_contradictory`, 1 `failed_ontarget_qc` and 4 `not_in_screen` TFs; GENIE3 had 3, 15, 1, 3 and 3 respectively. Predicted-target denominators differed substantially (median 9 targets/TF for GREmLN vs 18 for GENIE3), reflecting GENIE3's larger per-TF target sets. **Neither method met the superiority decision rule specified for the definitive comparison.** We do not claim the methods are equivalent: no equivalence test was performed, and the confirmation distributions simply overlap.

![Figure 3](../results/figures/fig3_crispri_target_response.png)

**Figure 3 (primary). CD4 CRISPRi functional target-response by TF and model (seed-excluded top-25, Stim8hr).** Fraction of each TF's predicted BTLA-DEG targets that are differentially expressed after that TF's knockdown (adj p < 0.05). Support is modest and method-comparable; unsigned, so direction is not established here.

**Table 1. Seed-excluded top-25 CRISPRi summary (Stim8hr, primary).**

| Metric | GREmLN | GENIE3 |
|---|---|---|
| Top-25 TFs | 25 | 25 |
| Present in screen | 21 | 22 |
| Passing on-target QC | 20 | 19 |
| Included in mean | 21 | 22 |
| Mean functional target-response fraction | 0.049 | 0.054 |
| Median | 0.00 | 0.00 |
| Bootstrap 95% CI (mean) | 0.00–0.13 | 0.01–0.11 |
| `usable_supportive` | 1 | 3 |
| `usable_unsupported` | 18 | 15 |
| `usable_contradictory` | 1 | 1 |
| `failed_ontarget_qc` | 1 | 3 |
| `not_in_screen` | 4 | 3 |
| Predicted-target denominators (median) | 9 | 18 |

Thresholds: custom target-response = adj p < 0.05; on-target = CZI Perturb-seq `ontarget_significant` (dataset-provided).

**Candidate shortlist audit (Table 2).** We audited eleven regulators spanning both methods (Table 2). This resolves an earlier over-claim: although the report highlights five GREmLN-specific candidates (AHR, NFIL3, PRDM1, RBPJ, STAT4), **only IRF8 (GREmLN side) and EGR2, BHLHE40, JUNB, REL and TBX21 (GENIE3 side / seeds) reach orthogonal CRISPRi functional support.** The five GREmLN-specific candidates are supported by literature and derived multi-omics but **do not have orthogonal CRISPRi functional support**: AHR shows a high target-response fraction (0.78) but in the *contradictory* direction (net opposing the BTLA programme); NFIL3 and PRDM1 are absent from the screen; RBPJ and STAT4 have zero confirmed targets. They are therefore hypothesis-generating leads requiring targeted perturbation, not perturbation-supported regulators. REL and TBX21 are BTLA panel genes (seeds), retained for context and not counted as candidate discoveries.

**Table 2. Candidate-level evidence audit.** Ranks are seed-excluded primary unless the TF is a seed (—); "func." = CRISPRi functional target-response.

| TF | Seed | GR rank (SE) | G3 rank (SE) | In screen | QC | pred | conf | frac | CRISPRi status | Indep. multi-omics | Derived multi-omics | Category |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| EGR2 | no | 6 | 1 | ✓ | ✓ | 31 | 12 | 0.39 | usable_supportive | none | transcript, TF-activity | functional + corroborating |
| IRF8 | no | 3 | 5 | ✓ | ✓ | 32 | 6 | 0.19 | usable_contradictory* | none | transcript, TF-activity | derived/lit only; no func. support |
| BHLHE40 | no | 5 | 7 | ✓ | ✓ | 29 | 13 | 0.45 | usable_supportive | none | transcript, TF-activity | functional + corroborating |
| JUNB | no | 3 | 2 | ✓ | ✓ | 38 | 2 | 0.05 | usable_supportive | none | transcript, TF-activity | functional + corroborating |
| REL | yes | — | — | ✓ | ✓ | 33 | 2 | 0.06 | usable_supportive | protein | transcript, TF-activity, BIONIC | functional + corroborating (seed) |
| TBX21 | yes | — | — | ✓ | ✓ | 14 | 4 | 0.29 | usable_supportive | none | transcript, TF-activity, BIONIC | functional + corroborating (seed) |
| AHR | no | 2 | 43 | ✓ | ✓ | 9 | 7 | 0.78 | usable_contradictory | none | TF-activity | derived/lit only; contradictory direction |
| NFIL3 | no | 2 | 60 | ✗ | — | 9 | — | — | not_in_screen | none | transcript, TF-activity | annotation/derived only; no func. support |
| PRDM1 | no | 3 | 25 | ✗ | — | 8 | — | — | not_in_screen | none | transcript, TF-activity | annotation/derived only; no func. support |
| RBPJ | no | 2 | 37 | ✓ | ✓ | 9 | 0 | 0.00 | usable_unsupported | protein | transcript, TF-activity | orthogonal (protein) + lit; no func. support |
| STAT4 | no | 2 | 111 | ✓ | ✓ | 9 | 0 | 0.00 | usable_unsupported | none | transcript, TF-activity | annotation/derived only; no func. support |

\*IRF8 status is model-specific: it reaches `usable_supportive` under GREmLN's predicted-target set but `usable_contradictory` under GENIE3's; the table shows the GENIE3-nominated evaluation. This model-dependence of CRISPRi status is itself a caveat.

**Literature (Paperclip).** A complete, template-matched Paperclip review now covers the **full union of both top-25 lists** (all reviewed under one identical query template, with a TF-focused second pass for canonical immune TFs; `results/paperclip/paperclip_union_top25_review.csv`). Of the seed-excluded top-25 candidates, strong-or-moderate T-cell literature exists for **5/7 shared** TFs (CREM, IRF8, JUNB, PRDM1, STAT1), **11/18 GENIE3-only** TFs (BHLHE40, EGR2, FOXM1, FOXP1, ID2, MAF, MDM2, MYC, NR4A2, STAT5B, ZEB2) and **6/18 GREmLN-only** TFs (AHR, BCL3, EZH2, NFIL3, RBPJ, STAT4). GENIE3's method-specific nominations are thus better represented in the prior literature than GREmLN's, and a substantial share of GREmLN-only picks are genes with little or no T-cell transcription-factor literature — several (KIF22, SRP9, SNAPC4, STAU2, HMGN3, HIRIP3, MLX, TFDP1, ETV3L, ZNF121, ZNF706, MSC) are housekeeping or non-transcription-factor genes, a coverage/quality caveat for the GREmLN-only set. Consistent with the protocol, literature is annotation, not a predictive input, and this count does not enter the CRISPRi-anchored verdict.

**Multi-omics support — independent versus derived.** We separate two categories rather than combining them into one count. **Independent orthogonal support** — protein abundance, phosphosite evidence, co-IP — is sparse among the audited candidates (e.g. protein-level support for REL and RBPJ) and is the only multi-omics evidence genuinely independent of the expression data driving both models. **Derived / contextual evidence** — transcript DE, inferred TF activity, inferred kinase activity, BIONIC/GNN modules — is broadly available but is derived from, or correlated with, the same transcriptomic signal, so it corroborates rather than independently validates. Collapsing the two would overstate independent support; we keep them distinct (Figure 4, Table 2).

![Figure 4](../results/figures/fig4_multiomics_heatmap.png)

**Figure 4. BTLA multi-omics support across the union of both top-25 lists.** Per-TF support by layer; ★ marks seeds. Independent orthogonal layers (protein, phospho, co-IP) are read separately from derived/contextual layers (transcript, TF/kinase activity, BIONIC).

**Supporting analyses (not part of the primary verdict).** Ranking TFs across three CD4 GENIE3 **construction/preprocessing variants** (NTC, ComBat-masked, per-batch integrated) — which are alternative preprocessing choices on largely overlapping cells, not independent GRNs — recovered a stable core of 16 TFs (EGR2, IRF8, BHLHE40, REL, HIVEP3, BACH2, TBX21, ZEB2, ZBTB32, NR4A3, CREM, MYC and others), a preprocessing-robustness check consistent with the GENIE3-side functional-support tier (Figure 5).

<!-- Figure 5 (3-GRN construction/preprocessing-variant stability) is a development-archive
     artifact and is not regenerated in this publication repo. See the archive tag
     development-archive-2026-07-12. -->
_Figure 5 (supporting): stability across three CD4 GENIE3 construction/preprocessing variants — available in the development archive (`development-archive-2026-07-12`)._

**Figure 5 (supporting). Stability across three CD4 GENIE3 construction/preprocessing variants.** Green = recovered in all three top-25 lists. These are preprocessing variants of one GENIE3 pipeline, not independent methods.

**Coverage.** GREmLN embeds 3,586 genes / 342 TFs; GENIE3 spans 3,602 genes / 335 regulator-TFs; the common scorable TF universe is 334. Of the 249 BTLA panel genes, 103 are available to GREmLN and 16 are TFs in the common universe. Coverage, not algorithm, is a binding constraint on GREmLN's reach here (Figure 6).

![Figure 6](../results/figures/fig5_coverage.png)

**Figure 6. Gene / TF universe coverage** for GENIE3 and GREmLN, and their overlap with the BTLA panel.

## Conclusions

Three conclusions follow. **On method:** a graph-aware single-cell foundation model, used zero-shot on a shared GENIE3 prior, is competitive with classical co-expression inference for BTLA-response TF prioritisation, but **neither method met the superiority decision rule specified for the definitive comparison** on orthogonal CD4 CRISPRi validation (mean functional target-response 0.049 vs 0.054, within the 0.02 margin). Because the comparison supplies GENIE3's graph to GREmLN and feeds the two models different numerical expression inputs, it measures the value of the GREmLN embedding procedure given that prior, not two independent methods. The likely payoff for foundation models here is expanding the inference universe and context-specific training, not the base checkpoint alone. **On biology:** the regulators with orthogonal CRISPRi functional support are the GENIE3-side activation TFs EGR2, BHLHE40 and JUNB and the panel TFs REL and TBX21; GREmLN-specific candidates (AHR, NFIL3, PRDM1, RBPJ, STAT4) are prioritised by embeddings and supported by literature/derived multi-omics but lack orthogonal perturbational support (AHR contradictory-direction; NFIL3/PRDM1 not in screen), and are hypothesis-generating leads for targeted perturbation. **On methodology:** separating independent orthogonal evidence from derived evidence, excluding response genes from the primary nomination, and reporting CRISPRi with QC/denominator/direction detail prevented several tempting over-claims — most notably that GREmLN-specific candidates were "causally validated".

The limitations are deliberate and stated plainly. The GENIE3 graph is upstream of GREmLN, so the comparison is of embedding value-add, not two independent methods on identical inputs; the common TF universe (334) limits reach; CRISPRi functional target-response fractions are modest, degree-sensitive and unsigned; and derived multi-omics/literature layers are annotation, not prediction. The held-out seed-TF recovery result (in which GREmLN reaches recall@25 = 1.00) is held in the supplement pending an audit for information leakage and inter-method asymmetry, and is not used in the abstract, main verdict or primary figures. SCENIC pruning of the persisted top-50k GENIE3 graph was retained as an exploratory motif-support sensitivity analysis and was not used as the GREmLN prior. The template-matched Paperclip review across the full candidate union is now complete (above); the natural next step is a multi-donor, CD4-trained GREmLN with a genome-scale universe, re-benchmarked against the same CRISPRi validation.

## References

1. Huynh-Thu VA, Irrthum A, Wehenkel L, Geurts P. *Inferring regulatory networks from expression data using tree-based methods* (GENIE3). PLoS ONE, 2010.
2. Van de Sande B, et al. *A scalable SCENIC workflow for single-cell gene regulatory network analysis* (pySCENIC). Nature Protocols, 2020.
3. Aibar S, et al. *SCENIC: single-cell regulatory network inference and clustering.* Nature Methods, 2017.
4. Zhang M, Swamy V, Cassius R, Dupire L, Paull E, AlQuraishi M, Karaletsos T, Califano A. *GREmLN: A Cellular Graph Structure Aware Transcriptomics Foundation Model.* bioRxiv 2025.07.03.663009 (2025). doi:10.1101/2025.07.03.663009.
5. Lambert SA, et al. *The Human Transcription Factors.* Cell, 2018.
6. Chan Zuckerberg Initiative. *Genome-scale T-cell Perturb-seq* dataset. https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq.
7. Conneau A, Lample G, Ranzato M, Denoyer L, Jégou H. *Word translation without parallel data* (introduces the CSLS metric; authors' default k = 10). ICLR, 2018.
8. Watanabe N, et al. *BTLA is a lymphocyte inhibitory receptor with similarities to CTLA-4 and PD-1.* Nature Immunology, 2003.

## Code and data availability

A dedicated, publication-facing repository ([`Lifework-Health/btla-gremln-benchmark`](https://github.com/Lifework-Health/btla-gremln-benchmark)) hosts the three executed notebooks (GENIE3 construction, GREmLN inference, benchmark), the model registry, small publication-permitted derived tables, the complete Paperclip literature review, and figures; GREmLN is pinned as a git submodule. The content is being introduced as a **draft pull request** and is **not yet a public release**: BTLA transcriptomic and multi-omic **source data have access restrictions pending clearance**, so no BTLA source or restricted derived files are committed (they regenerate locally from a documented, access-gated data root), whereas the CD4 Perturb-seq data are public (CZI, MIT license). A tagged release will follow once data-publication permission is confirmed. The development archive is retained separately under the tag `development-archive-2026-07-12`.
