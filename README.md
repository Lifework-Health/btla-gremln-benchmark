# btla-gremln-benchmark

**Benchmarking GREmLN against GENIE3 for BTLA-response transcription-factor prioritisation using orthogonal CRISPRi validation.**

> Status: **report draft (`v1.0.0-report-draft`).** Code is public; **BTLA data are access-restricted
> and not redistributed here** (see [`data/README.md`](data/README.md)). This is a clean,
> publication-facing extract of a larger development archive.

## Scientific question

Does the graph-aware single-cell foundation model **GREmLN** prioritise the transcription-factor
(TF) regulators of the **BTLA cross-linking response** in human CD4⁺ T cells more accurately than
classical co-expression network inference (**GENIE3**), when judged against *orthogonal*
perturbation data?

## What is compared

Two models built on the **same** regulatory graph:

1. **GENIE3 (masked, raw)** — a 50,000-edge CD4 co-expression GRN inferred from donor-adjusted
   (ComBat) expression. TFs scored by **summed outgoing TF→BTLA-DEG edge weight**.
2. **GREmLN (masked prior)** — the GREmLN foundation model run zero-shot on the *pre-ComBat*
   log1p-CPM expression of the same cells, using the masked GENIE3 graph as its **regulatory
   prior**. TFs scored by **CSLS** (k=10) BTLA-seed neighbour count.

This measures the value of the GREmLN embedding procedure *given the GENIE3 prior* — **not** two
independent methods on identical numerical inputs.

## Main conclusion

Over the common 334-TF universe the two rankings are moderately correlated (Spearman ρ ≈ 0.46) and
share 7/25 primary (seed-excluded) candidate TFs. On the keystone **CD4 CRISPRi Perturb-seq**
validation, functional target-response support is modest and comparable (mean confirmed-target
fraction 0.049 GREmLN vs 0.054 GENIE3); **neither method meets the superiority decision rule
specified for the comparison.** We report a prioritised, evidence-tiered BTLA-response TF shortlist
rather than a single winner. Held-out seed-TF recovery is **supplementary only** (pending a
leakage/asymmetry audit) and does not contribute to the verdict.

## Repository structure

```
btla-gremln-benchmark/
├── README.md                     # this file
├── LICENSE                       # MIT (code only; see file for data/submodule scope)
├── CITATION.cff
├── environment.yml               # conda env (fresh build)
├── requirements-lock.txt         # exact pinned versions (reproduce released outputs)
├── .gitmodules                   # pins third_party/GREmLN at a fixed commit
├── config/benchmark_config.yaml  # all parameters + relative paths (no absolute paths)
├── data/
│   ├── README.md                 # data provenance, access, publication gate
│   ├── manifests/                # inputs_manifest.csv (checksums, dims, source, access)
│   └── derived/                  # small permitted derived tables (post-clearance)
├── resources/human_tfs_pySCENIC.txt
├── notebooks/
│   ├── 01_build_genie3_cd4_masked.ipynb
│   ├── 02_build_gremln_from_masked_genie3_prior.ipynb
│   └── 03_benchmark_btla_tf_prioritisation.ipynb
├── scripts/bench_utils.py        # small IO/plot/linear-algebra helpers only
├── results/{model_registry,tables,figures}/   # mostly regenerated (gitignored until clearance)
├── report/btla_gremln_benchmark.md
└── third_party/GREmLN/           # submodule, pinned commit
```

## Environment

```bash
git clone --recurse-submodules <repo-url>
cd btla-gremln-benchmark
conda env create -f environment.yml
conda activate btla-gremln-benchmark
# GPU PyTorch (notebook 02 only) — install the build matching your CUDA, e.g.:
pip install torch==2.4.1 --index-url https://download.pytorch.org/whl/cu121
```

If you already cloned without submodules: `git submodule update --init --recursive`.

## Data access

See [`data/README.md`](data/README.md). Briefly: the CD4 Perturb-seq data are public (CZI, MIT);
the GREmLN checkpoint is public (large, download separately); the **BTLA transcriptomic/multi-omic
data are access-restricted** and gated behind a publication-permission review. No BTLA source or
restricted derived files are committed.

## Notebook execution order, runtime, outputs

| # | Notebook | Does | Needs | Runtime |
|---|---|---|---|---|
| 01 | `01_build_genie3_cd4_masked.ipynb` | CD4 NTC selection, donor/batch audit, normalisation, ComBat masking, TF universe, GENIE3 → 50k-edge graph, model registry, SCENIC top-50k sensitivity | CD4 Perturb-seq NTC cells | Hours (GENIE3) — reuse cell loads the persisted graph |
| 02 | `02_build_gremln_from_masked_genie3_prior.ipynb` | Checkpoint provenance, expression-input rationale, graph-prior + vocab/tokenizer, GREmLN gene-embedding inference, CSLS, coverage audit, saved embeddings/registry | GREmLN checkpoint, GPU | ~minutes on GPU — reuse cell loads saved embeddings |
| 03 | `03_benchmark_btla_tf_prioritisation.ipynb` | BTLA panel (249), common TF universe, **seed-excluded primary** + seed-inclusive context rankings, top-25 comparison, overlap/rank correlation, CD4 CRISPRi, Paperclip, independent-vs-derived multi-omics, SCENIC summary, all report figures/tables, restrained verdict | outputs of 01 + 02 | CPU, minutes |

Execution is controlled by a `DATA_ROOT` (env var / config) pointing at the obtained data and the
regenerated artifacts of the previous notebook. Heavy steps (GENIE3 inference, GPU embedding) are
behind a `REGENERATE` flag; by default the notebooks **reuse validated artifacts** and run quickly.

## Expected outputs

`results/tables/` (TF rankings, top-25 comparison, rank correlation, CRISPRi validation, integrated
evidence, candidate audit), `results/figures/` (top-25 comparison, overlap, CRISPRi target-response,
multi-omics heatmap, coverage; held-out recovery is supplementary), `results/model_registry/`
(provenance). Most are regenerated locally and remain untracked until BTLA data clearance.

## Reproduction status

* **Fully reproducible with public data:** submodule pin, environment, GENIE3 graph construction,
  GREmLN inference (given the public checkpoint + CD4 Perturb-seq download).
* **Requires restricted data:** the BTLA panel and multi-omics joins, hence the final benchmark
  tables/figures. These regenerate once BTLA data access is granted.
* **Not part of the verdict:** held-out seed-TF recovery (supplementary; pending audit).

## License and citations

Code: MIT (see `LICENSE`). `third_party/GREmLN` and the datasets keep their own licenses. Please
cite this repository (`CITATION.cff`) and the GREmLN, GENIE3, CSLS and CZI Perturb-seq sources.

## BTLA data access restrictions

BTLA transcriptomic and multi-omic data are **not** included and **not** redistributable here. The
benchmark's BTLA-dependent outputs can only be regenerated by users with authorised access.
