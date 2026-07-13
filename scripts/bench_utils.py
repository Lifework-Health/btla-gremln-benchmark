"""Small IO / plotting / provenance helpers for the BTLA GREmLN-vs-GENIE3 benchmark notebooks.

Scientific logic (ranking formulae, CSLS, CRISPRi validation, evidence integration, the verdict)
lives IN the notebooks so it is readable. This module only holds path resolution, checksums,
data loaders, and figure rendering shared across notebooks 01/02/03.

Paths are resolved from a single DATA_ROOT so the notebooks contain no absolute local paths:
  * env var  BTLA_BENCH_DATA_ROOT   (preferred), or
  * config/benchmark_config.yaml -> data_root, or
  * the repo root (fallback).
Obtain/regenerate the referenced artifacts per data/README.md.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "notebooks").is_dir() and (c / "config").is_dir():
            return c
    return here.parent.parent


def data_root() -> Path:
    env = os.environ.get("BTLA_BENCH_DATA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    cfg = repo_root() / "config" / "benchmark_config.yaml"
    if cfg.exists():
        try:
            import yaml
            d = yaml.safe_load(cfg.read_text()) or {}
            if d.get("data_root"):
                return Path(str(d["data_root"])).expanduser().resolve()
        except Exception:
            pass
    return repo_root()


# Logical artifact -> path template relative to DATA_ROOT. The development-archive layout is the
# default; override any of these with BTLA_BENCH_<NAME> env vars if your local layout differs.
_LAYOUT = {
    "gremln_emb": "runs/gremln_cd4_masked_raw_prior/2026-07-11_masked_raw_prior_v01/gremln_gene_emb.h5ad",
    "gremln_inputs": "runs/gremln_cd4_masked_raw_prior/2026-07-11_masked_raw_prior_v01/gremln_inputs",
    "gremln_run": "runs/gremln_cd4_masked_raw_prior/2026-07-11_masked_raw_prior_v01",
    "genie3_edges": "runs/genie3_cd4/2026-07-07_cd4_perturb_targeting_masked_v01/grn_edges.tsv",
    "cd4_ntc_counts": "data/raw/cd4_perturb/cd4_ntc_counts.h5ad",
    "gremln_checkpoint": "third_party/GREmLN/checkpoints/model.ckpt",
    "tf_list": "resources/human_tfs_pySCENIC.txt",
    "deg_table": "results/btla_csls_multiomics/btla_4h_deg_gene_level_table.csv",
    "de_stats": "data/raw/cd4_perturb/GWCD4i.DE_stats.h5ad",
    "multiomics_long": "results/btla_csls_multiomics/btla_candidate_multiomics_long_evidence.csv",
    "multiomics_summary": "results/btla_csls_multiomics/btla_candidate_multiomics_summary.csv",
    "paperclip_gremln": "results/btla_csls_multiomics/btla_tf_gremln_top25_paperclip.csv",
    "paperclip_genie3": "results/btla_csls_multiomics/btla_tf_genie3_top25_paperclip.csv",
}


def artifact(name: str) -> Path:
    env = os.environ.get(f"BTLA_BENCH_{name.upper()}")
    if env:
        return Path(env).expanduser().resolve()
    return data_root() / _LAYOUT[name]


def md5(path: Path, limit: int = 2_000_000_000) -> str:
    path = Path(path)
    if not path.exists():
        return "missing"
    if path.stat().st_size > limit:
        return f"skipped_large_{path.stat().st_size}"
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def redact(p) -> str:
    """Render a path relative to DATA_ROOT so committed provenance carries no absolute paths."""
    try:
        return str(Path(p).resolve().relative_to(data_root()))
    except Exception:
        return Path(p).name


# --------------------------------------------------------------------------- loaders
def load_tfs(path: Path = None) -> set:
    # The TF list is a small public file committed to the repo; prefer that copy so both a
    # clean public clone and a reuse-mode run (BTLA_BENCH_DATA_ROOT set) resolve it correctly.
    if path is None:
        env = os.environ.get("BTLA_BENCH_TF_LIST")
        if env:
            path = Path(env).expanduser()
        else:
            repo_copy = repo_root() / "resources" / "human_tfs_pySCENIC.txt"
            path = repo_copy if repo_copy.exists() else artifact("tf_list")
    return {ln.strip() for ln in Path(path).read_text().splitlines()
            if ln.strip() and not ln.startswith("#")}


def load_edges(path: Path = None) -> pd.DataFrame:
    path = path or artifact("genie3_edges")
    e = pd.read_csv(path, sep="\t")
    e["regulator"] = e["regulator"].astype(str)
    e["target"] = e["target"].astype(str)
    if "weight" not in e.columns:
        e["weight"] = 1.0
    e["weight"] = e["weight"].astype(float)
    return e


def load_gremln_embeddings(path: Path = None) -> Tuple[List[str], np.ndarray]:
    import anndata as ad
    path = path or artifact("gremln_emb")
    g = ad.read_h5ad(path)
    X = np.asarray(g.X, dtype=np.float64)
    hugo = g.obs["hugo"].astype(str).values
    seen, genes, rows = set(), [], []
    for i, sym in enumerate(hugo):
        if sym and sym not in ("nan", "None") and sym not in seen:
            seen.add(sym); genes.append(sym); rows.append(i)
    return genes, X[rows]


def load_seeds(path: Path = None) -> Tuple[set, Dict[str, str], pd.DataFrame]:
    path = path or artifact("deg_table")
    deg = pd.read_csv(path)
    s = deg[deg["contrast_id"] == "BTLA_TCR_vs_TCR"].copy()
    s["gene"] = s["gene"].astype(str)
    seeds = set(s["gene"])
    direction = dict(zip(s["gene"], s["direction"]))
    return seeds, direction, s


# --------------------------------------------------------------------------- figures
def _save(fig, path):
    fig.tight_layout(); fig.savefig(path, dpi=200, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)


def fig_top25_comparison(gr_top25, g3_top25, out_path, seeds=None, mode="seed-excluded"):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    seeds = set(seeds) if seeds is not None else set()
    fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), constrained_layout=True)
    panels = (
        (axes[0], gr_top25, "gremln_csls_score", "GREmLN — CSLS (masked prior)",
         "#c0392b", "#f1948a", "CSLS seed-neighbour score"),
        (axes[1], g3_top25, "genie3_score", "GENIE3 — masked raw",
         "#1f618d", "#7fb3d5", "Summed outgoing TF\u2192BTLA-DEG edge weight"),
    )
    for ax, df, scol, title, cand_c, seed_c, xlab in panels:
        d = df.head(25).iloc[::-1]
        y = np.arange(len(d))
        is_seed = d["gene"].isin(seeds).values
        colors = [(seed_c if s else cand_c) for s in is_seed]
        vmax = float(d[scol].max()) if len(d) else 1.0
        ax.barh(y, d[scol].values, color=colors, edgecolor="white", linewidth=0.5, zorder=3)
        for yi, val in enumerate(d[scol].values):
            ax.text(val + vmax * 0.01, yi, f"{val:g}", va="center", ha="left", fontsize=6.5,
                    color="#333", zorder=4)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{g} \u2605" if s else g for g, s in zip(d["gene"].values, is_seed)],
                           fontsize=8)
        ax.set_title(title, fontsize=12, fontweight="bold", pad=8)
        ax.set_xlabel(xlab, fontsize=9); ax.set_xlim(0, vmax * 1.14); ax.margins(y=0.01)
        ax.grid(axis="x", color="#e6e6e6", linewidth=0.8, zorder=0)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=0)
    handles = [Patch(facecolor="#555", label="candidate (non-seed, darker bar)"),
               Patch(facecolor="#bdbdbd", label="\u2605 BTLA_vs_TCR seed (lighter bar)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"{mode.capitalize()} top-25 transcription factors (common TF universe)",
                 fontsize=13, fontweight="bold")
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_overlap(gr_set, g3_set, out_path, seeds=None, title="Top-25 overlap"):
    import matplotlib.pyplot as plt
    seeds = set(seeds) if seeds is not None else set()
    mark = lambda g: f"\u2605{g}" if g in seeds else g
    shared = sorted(set(gr_set) & set(g3_set))
    gr_only = sorted(set(gr_set) - set(g3_set))
    g3_only = sorted(set(g3_set) - set(gr_set))
    n_seed_shared = len(set(shared) & seeds)
    try:
        from matplotlib_venn import venn2
        fig, ax = plt.subplots(figsize=(7, 6))
        venn2([set(gr_set), set(g3_set)], set_labels=("GREmLN CSLS", "GENIE3 masked"), ax=ax)
        ax.set_title(f"{title}\n\u2605 seeds among shared: {n_seed_shared}/{len(shared)}")
    except Exception:
        fig, ax = plt.subplots(figsize=(8, 6)); ax.axis("off")
        txt = (f"Shared ({len(shared)}; \u2605 seeds={n_seed_shared}):\n  " + ", ".join(mark(g) for g in shared) +
               f"\n\nGREmLN-only ({len(gr_only)}):\n  " + ", ".join(mark(g) for g in gr_only) +
               f"\n\nGENIE3-only ({len(g3_only)}):\n  " + ", ".join(mark(g) for g in g3_only))
        ax.text(0.02, 0.98, txt, va="top", ha="left", fontsize=10, family="monospace")
        ax.set_title(f"{title}  (\u2605 = BTLA_vs_TCR seed)")
    _save(fig, out_path)


def fig_crispri(crispri, out_path):
    import matplotlib.pyplot as plt
    d = crispri[crispri["crispri_arm_present_8hr"] == True].copy()
    if d.empty:
        fig, ax = plt.subplots(figsize=(6, 4)); ax.text(0.5, 0.5, "no CRISPRi arms", ha="center")
        _save(fig, out_path); return
    seed_tfs = set(d.loc[d.get("is_BTLA_vs_TCR_seed", False) == True, "TF"]) \
        if "is_BTLA_vs_TCR_seed" in d.columns else set()
    piv = d.pivot_table(index="TF", columns="model", values="frac_confirmed_8hr", aggfunc="max").fillna(0)
    piv = piv.sort_values(list(piv.columns)[0], ascending=True)
    fig, ax = plt.subplots(figsize=(10, max(6, len(piv) * 0.32)))
    y = np.arange(len(piv)); colors = {"GREmLN": "#d6604d", "GENIE3": "#4393c3"}; ncol = len(piv.columns)
    for i, model in enumerate(piv.columns):
        ax.barh(y + (i - (ncol - 1) / 2) * 0.4, piv[model].values, 0.4, label=model,
                color=colors.get(model, "#888"))
    ax.set_yticks(y)
    ax.set_yticklabels([f"\u2605 {t}" if t in seed_tfs else t for t in piv.index], fontsize=8)
    ax.set_xlabel("CD4 CRISPRi functional target-response fraction (Stim8hr)")
    ax.set_title("CD4 CRISPRi functional target-response by TF and model (\u2605 = seed)")
    ax.legend(); _save(fig, out_path)


def fig_evidence_heatmap(integrated, out_path):
    import matplotlib.pyplot as plt
    layers = ["transcript_support", "protein_support", "phosphosite_support", "coip_support",
              "tf_activity_support", "kinase_activity_support", "bionic_support",
              "early_synapse_trafficking_support"]
    layers = [c for c in layers if c in integrated.columns]
    if not layers:
        fig, ax = plt.subplots(figsize=(6, 4)); ax.text(0.5, 0.5, "no multi-omics", ha="center")
        _save(fig, out_path); return
    d = integrated.set_index("TF")

    def to_num(v):
        s = str(v).lower()
        if s in ("nan", "none", ""): return 0.0
        if s in ("strong", "high", "yes", "true", "significant"): return 2.0
        if s in ("moderate", "medium", "weak", "partial"): return 1.0
        try: return float(v)
        except Exception: return 1.0 if s not in ("no", "false", "none") else 0.0

    seed_tfs = set(d.index[d["is_BTLA_vs_TCR_seed"] == True]) if "is_BTLA_vs_TCR_seed" in d.columns else set()
    mat = np.array([[to_num(d.loc[tf, l]) for l in layers] for tf in d.index])
    fig, ax = plt.subplots(figsize=(max(8, len(layers) * 1.1), max(6, len(d) * 0.3)))
    im = ax.imshow(mat, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(layers)))
    ax.set_xticklabels([l.replace("_support", "") for l in layers], rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([f"\u2605 {t}" if t in seed_tfs else t for t in d.index], fontsize=7)
    ax.set_title("Integrated BTLA multi-omics support (union top-25; \u2605 = seed)")
    fig.colorbar(im, ax=ax, shrink=0.6, label="support (0/1/2)")
    _save(fig, out_path)


def fig_coverage(n_expr, n_genie3, n_gremln, gremln_tfs, genie3_tfs, seeds_common, out_path):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    v1 = [n_expr, n_genie3, n_gremln]
    axes[0].bar(["masked expr", "GENIE3 graph", "GREmLN embedded"], v1,
                color=["#999", "#4393c3", "#d6604d"])
    for i, v in enumerate(v1):
        axes[0].text(i, v, str(v), ha="center", va="bottom", fontsize=9)
    axes[0].set_title("Gene-universe coverage"); axes[0].set_ylabel("n genes")
    common = len(set(gremln_tfs) & set(genie3_tfs))
    v2 = [len(genie3_tfs), len(gremln_tfs), common, seeds_common]
    axes[1].bar(["GENIE3 TFs", "GREmLN TFs", "common TFs", "seeds\u2229common"], v2,
                color=["#4393c3", "#d6604d", "#7fbf7b", "#b2abd2"])
    for i, v in enumerate(v2):
        axes[1].text(i, v, str(v), ha="center", va="bottom", fontsize=9)
    axes[1].set_title("TF universe coverage"); axes[1].tick_params(axis="x", labelrotation=20)
    _save(fig, out_path)
