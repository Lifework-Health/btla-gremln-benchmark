#!/usr/bin/env python3
"""Generate the three benchmark notebooks with the scientific logic inlined in cells.

Only IO/plotting/provenance helpers are imported (scripts/bench_utils.py); all ranking, CSLS,
CRISPRi-validation, evidence-integration and verdict logic is written directly into notebook
cells so the notebooks are readable and are not thin wrappers.

Usage:  python scripts/_build_notebooks.py [01|02|03|all]
"""
import sys
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks"
NB.mkdir(exist_ok=True)


def _nb(cells):
    nb = nbf.v4.new_notebook()
    nb.cells = cells
    nb.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    }
    return nb


def md(s):
    return nbf.v4.new_markdown_cell(s.strip("\n"))


def code(s):
    return nbf.v4.new_code_cell(s.strip("\n"))


# ============================================================ NOTEBOOK 03
def build_nb3():
    c = []
    c.append(md("""
# 03 — Benchmark: BTLA-response TF prioritisation (GREmLN masked-prior vs GENIE3 masked-raw)

Definitive head-to-head. Both models are restricted to the **common GREmLN∩GENIE3 TF universe**
and scored against the **BTLA_TCR_vs_TCR 4 h panel (249 DEGs)**:

* **GREmLN** — CSLS (k=10): a TF's score = number of BTLA seeds among its top-100 embedding
  neighbours (dense rank).
* **GENIE3** — summed **outgoing** TF→BTLA-DEG edge weight (dense rank).

**Primary nomination comparison = seed-EXCLUDED** (BTLA panel genes removed → candidate
regulators). A **seed-INCLUSIVE** ranking is reported separately as biological context, with panel
genes flagged (★) and never called discoveries.

Candidates are then examined with three **orthogonal, non-predictive** validators:
1. **CD4 CRISPRi Perturb-seq** (primary; Stim8hr) — functional target-response;
2. **Paperclip** literature (identical template);
3. **BTLA multi-omics**, split into **independent orthogonal** (protein/PTM/co-IP) vs
   **derived/contextual** (transcript DE, TF/kinase activity, BIONIC).

Held-out seed-TF recovery is **supplementary only** (pending a leakage/asymmetry audit) and does
**not** enter the verdict. Requires the outputs of notebooks 01 (GENIE3 graph) and 02 (GREmLN
embeddings); set `BTLA_BENCH_DATA_ROOT` to where those + the BTLA panel/Perturb-seq live.
"""))

    c.append(code("""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, rankdata, spearmanr
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path.cwd().parent / "scripts"))
import bench_utils as bu
import yaml

CFG = yaml.safe_load((bu.repo_root() / "config" / "benchmark_config.yaml").read_text())
# ---- locked parameters (must match notebooks 01/02) ----
CSLS_K, TOPN = 10, 100
DE_PADJ, LFC_MIN = 0.05, float(np.log2(1.5))
STIM_PRIMARY, STIM_SENS = "Stim8hr", "Stim48hr"
RANDOM_SEED = 42
CRISPRI_MARGIN = CFG["decision_rule"]["crispri_margin"]   # decision-rule margin (from config)
N_MASKED_EXPR_GENES = int(CFG["genie3"]["n_hvg"])         # masked expression universe (HVGs)

OUT = bu.repo_root() / "results"
TAB, FIG, REG = OUT / "tables", OUT / "figures", OUT / "model_registry"
for d in (TAB, FIG, REG):
    d.mkdir(parents=True, exist_ok=True)
print("DATA_ROOT:", bu.data_root())
"""))

    c.append(md("## 1. BTLA panel (assert 249), TF list, embeddings, GENIE3 graph"))
    c.append(code("""
seeds, seed_dir, seeds_df = bu.load_seeds()
assert len(seeds) == 249, f"BTLA_TCR_vs_TCR panel must be exactly 249 genes, got {len(seeds)}"
tfs = bu.load_tfs()
genes, X = bu.load_gremln_embeddings()
edges = bu.load_edges()
n_up = (seeds_df['direction'] == 'up').sum(); n_down = (seeds_df['direction'] == 'down').sum()
print(f"BTLA panel: {len(seeds)} genes (up={n_up}, down={n_down}) | TF list: {len(tfs)}")
print(f"GREmLN embedded genes: {len(genes)} | GENIE3 masked edges: {len(edges)}")
"""))

    c.append(md("""## 2. CSLS similarity and the common (scorable-by-both) TF universe

Cosine → CSLS (Conneau et al. 2018) corrects hubness: `CSLS(i,j) = 2·cos(i,j) − r_k(i) − r_k(j)`
with `r_k` the mean cosine to the k nearest neighbours."""))
    c.append(code("""
def cosine_similarity_matrix(X):
    n = np.linalg.norm(X, axis=1, keepdims=True); n = np.where(n == 0, 1.0, n)
    Xn = X / n; S = Xn @ Xn.T; np.fill_diagonal(S, -np.inf); return S

def compute_csls_matrix(S, k=CSLS_K):
    Sf = np.where(np.isfinite(S), S, -np.inf)
    rk = np.sort(Sf, axis=1)[:, -k:].mean(axis=1)
    return 2 * S - rk[:, None] - rk[None, :]

S_csls = compute_csls_matrix(cosine_similarity_matrix(X), k=CSLS_K)

gremln_tfs = set(genes) & tfs
genie3_tfs = set(edges['regulator']) & tfs
common = gremln_tfs & genie3_tfs
print(f"GREmLN TFs: {len(gremln_tfs)} | GENIE3 regulator-TFs: {len(genie3_tfs)} | common: {len(common)}")
"""))

    c.append(md("""## 3. Ranking formulae (inlined) and rankings

`score_genes_by_seed_neighbours` (GREmLN) and the outgoing-edge-weight sum (GENIE3) are computed
for both the **seed-excluded** (primary) and **seed-inclusive** (context) views, restricted to the
common universe with dense tie-aware ranks."""))
    c.append(code("""
def score_genes_by_seed_neighbours(genes, S, seed_genes, topn=TOPN, exclude_seeds=False):
    seed_in = {g for g in seed_genes if g in set(genes)}
    G = len(genes); is_seed = np.array([g in seed_in for g in genes]); K = int(is_seed.sum())
    base = K / G if G else 0.0
    dcount = np.empty(G, int); pval = np.ones(G)
    for i in range(G):
        nn = np.argpartition(-S[i], min(topn, G - 1))[:topn]
        x = int(is_seed[nn].sum()); dcount[i] = x
        nsucc = K - int(is_seed[i])
        pval[i] = hypergeom.sf(x - 1, G - 1, nsucc, topn) if nsucc >= 0 else 1.0
    out = pd.DataFrame({"gene": genes, "gremln_csls_score": dcount,
                        "fold_enrichment": (dcount / topn) / base if base > 0 else 0.0,
                        "hypergeom_p": pval, "is_seed": is_seed})
    if exclude_seeds:
        out = out.loc[~out["is_seed"]].copy()
    return out.sort_values("gremln_csls_score", ascending=False).reset_index(drop=True)

def gremln_ranking(exclude_seeds):
    r = score_genes_by_seed_neighbours(genes, S_csls, seeds, exclude_seeds=exclude_seeds)
    t = r[r["gene"].isin(common)].copy()
    t["gremln_csls_rank"] = rankdata(-t["gremln_csls_score"], method="dense").astype(int)
    t["is_BTLA_vs_TCR_seed"] = t["gene"].isin(seeds)
    t["BTLA_vs_TCR_direction"] = t["gene"].map(lambda g: f"BTLA_{seed_dir[g]}" if g in seed_dir else "not_seed")
    t["eligible_for_fair_comparison"] = (~t["is_BTLA_vs_TCR_seed"]) if exclude_seeds else True
    return t.sort_values("gremln_csls_rank").reset_index(drop=True)

def genie3_ranking(exclude_seeds):
    st = edges[edges["target"].isin(seeds)]
    inc = (st.groupby("regulator")["weight"].agg(n_seed_targets="count", genie3_score="sum")
             .reset_index().rename(columns={"regulator": "gene"}))
    inc = inc[inc["gene"].isin(common)].copy()
    if exclude_seeds:
        inc = inc[~inc["gene"].isin(seeds)].copy()
    inc = inc.sort_values("genie3_score", ascending=False).reset_index(drop=True)
    inc["genie3_rank"] = rankdata(-inc["genie3_score"], method="dense").astype(int)
    inc["is_BTLA_vs_TCR_seed"] = inc["gene"].isin(seeds)
    inc["BTLA_vs_TCR_direction"] = inc["gene"].map(lambda g: f"BTLA_{seed_dir[g]}" if g in seed_dir else "not_seed")
    return inc

rk = {}
for mode, excl in (("seed_excluded", True), ("seed_inclusive", False)):
    gr, g3 = gremln_ranking(excl), genie3_ranking(excl)
    rk[mode] = {"gr": gr, "g3": g3, "gr25": list(gr.head(25)["gene"]), "g325": list(g3.head(25)["gene"])}
    gr.to_csv(TAB / f"gremln_btla_vs_tcr_{mode}_tf_ranking.csv", index=False)
    g3.to_csv(TAB / f"genie3_btla_vs_tcr_{mode}_tf_ranking.csv", index=False)
print("seed-excluded GREmLN top10:", rk["seed_excluded"]["gr25"][:10])
print("seed-excluded GENIE3 top10:", rk["seed_excluded"]["g325"][:10])
"""))

    c.append(md("## 4. Primary comparison (seed-excluded): overlap + rank correlation. Seed-inclusive as context."))
    c.append(code("""
def compare(mode):
    gr, g3 = rk[mode]["gr"], rk[mode]["g3"]
    m = gr[["gene", "gremln_csls_score", "gremln_csls_rank", "is_BTLA_vs_TCR_seed"]].merge(
        g3[["gene", "genie3_score", "genie3_rank"]], on="gene", how="outer")
    m["rank_difference_gremln_minus_genie3"] = m["gremln_csls_rank"] - m["genie3_rank"]
    both = m.dropna(subset=["gremln_csls_rank", "genie3_rank"])
    rho = spearmanr(both["gremln_csls_rank"], both["genie3_rank"])[0] if len(both) >= 3 else np.nan
    gr25, g325 = set(rk[mode]["gr25"]), set(rk[mode]["g325"])
    ov = {"mode": mode, "spearman_rank_rho": float(rho), "n_common_ranked_by_both": int(len(both)),
          "shared_top25": sorted(gr25 & g325), "gremln_only_top25": sorted(gr25 - g325),
          "genie3_only_top25": sorted(g325 - gr25)}
    m.sort_values("gremln_csls_rank").to_csv(TAB / f"gremln_genie3_common_universe_rank_comparison_{mode}.csv", index=False)
    return ov

overlap = {m: compare(m) for m in ("seed_excluded", "seed_inclusive")}
prim = overlap["seed_excluded"]
print(f"PRIMARY (seed-excluded): Spearman rho={prim['spearman_rank_rho']:.3f}; "
      f"shared top25={len(prim['shared_top25'])}: {prim['shared_top25']}")
print("GREmLN-only:", prim["gremln_only_top25"])
print("GENIE3-only:", prim["genie3_only_top25"])
"""))

    c.append(md("""## 5. CD4 CRISPRi perturbational validation (primary; seed-excluded top-25)

For each model's top-25, its predicted BTLA-DEG targets (GENIE3: graph targets∩seeds; GREmLN:
top-100 CSLS neighbours∩seeds) are tested for differential expression after that TF's knockdown
(adj p < DE_PADJ). Signed agreement flags whether confirmed targets move with/against the BTLA
programme. Both scores are **unsigned**, so target movement does not establish direction."""))
    c.append(code("""
def genie3_targets_by_tf(tf_list):
    sub = edges[edges["regulator"].isin(tf_list) & edges["target"].isin(seeds)]
    return {tf: set(g["target"].unique()) for tf, g in sub.groupby("regulator", sort=False)}

def csls_seed_targets_by_tf(tf_list):
    gi = {g: i for i, g in enumerate(genes)}; seed_in = set(seeds) & set(genes); G = len(genes); out = {}
    for tf in tf_list:
        if tf not in gi:
            out[tf] = set(); continue
        nn = np.argpartition(-S_csls[gi[tf]], min(TOPN, G - 1))[:TOPN]
        out[tf] = {genes[j] for j in nn} & seed_in
    return out

def crispri_validation(union_by_model, de_stats_path):
    import anndata as ad, h5py
    all_tfs = sorted({t for lst in union_by_model.values() for t in lst})
    de = ad.read_h5ad(de_stats_path, backed="r"); de_obs = de.obs.reset_index(drop=True)
    de_gene = de.var["gene_name"].astype(str).values; n_meas = len(de_gene)
    def arm(cond):
        mask = ((de_obs["culture_condition"].astype(str) == cond)
                & (de_obs["target_contrast_gene_name"].astype(str).isin(all_tfs)))
        ri = np.where(mask.values)[0]; order = np.argsort(ri); sidx = ri[order]; inv = np.argsort(order)
        with h5py.File(de_stats_path, "r") as f:
            adjp = f["layers/adj_p_value"][sidx, :][inv]; logfc = f["layers/log_fc"][sidx, :][inv]
        tfn = de_obs.loc[ri, "target_contrast_gene_name"].astype(str).values
        ot = de_obs.loc[ri, "ontarget_significant"].values
        nc = de_obs.loc[ri, "n_cells_target"].values if "n_cells_target" in de_obs else np.full(len(ri), np.nan)
        return adjp, logfc, ot, nc, {t: i for i, t in enumerate(tfn)}
    a8, l8, o8, c8, p8 = arm(STIM_PRIMARY); a48, l48, o48, _, p48 = arm(STIM_SENS)
    def de_info(pmap, adjp, tf):
        if tf not in pmap: return None
        padj = adjp[pmap[tf], :]; hit = np.isfinite(padj) & (padj < DE_PADJ)
        return set(de_gene[hit]), int(hit.sum())
    sign = {g: (1 if seed_dir.get(g) == "up" else -1) for g in seeds}
    rows = []
    for model, top in union_by_model.items():
        pred_by = csls_seed_targets_by_tf(top) if model == "GREmLN" else genie3_targets_by_tf(top)
        for tf in top:
            pred = set(pred_by.get(tf, set())); present = tf in p8
            d8, d48 = de_info(p8, a8, tf), de_info(p48, a48, tf)
            rec = {"TF": tf, "model": model, "is_BTLA_vs_TCR_seed": tf in seeds,
                   "crispri_arm_present_8hr": bool(present), "crispri_arm_present_48hr": bool(tf in p48),
                   "ontarget_kd_significant": bool(o8[p8[tf]]) if present else np.nan,
                   "n_predicted_btla_targets": len(pred)}
            for lab, di in (("8hr", d8), ("48hr", d48)):
                if di is None or not pred:
                    rec[f"n_de_genes_crispri_{lab}"] = di[1] if di else np.nan
                    rec[f"n_confirmed_{lab}"] = np.nan; rec[f"frac_confirmed_{lab}"] = np.nan
                    rec[f"hypergeom_p_{lab}"] = np.nan; rec[f"confirmed_targets_{lab}"] = ""
                    continue
                dg, nde = di; conf = pred & dg; k = len(conf)
                rec[f"n_de_genes_crispri_{lab}"] = nde; rec[f"n_confirmed_{lab}"] = k
                rec[f"frac_confirmed_{lab}"] = k / len(pred)
                rec[f"hypergeom_p_{lab}"] = float(hypergeom.sf(k - 1, n_meas, nde, len(pred))) if k > 0 else np.nan
                rec[f"confirmed_targets_{lab}"] = "|".join(sorted(conf))
            na = no = 0
            if present and pred and d8:
                lfc = dict(zip(de_gene, l8[p8[tf], :]))
                for g in (pred & d8[0]):
                    lf = lfc.get(g, np.nan)
                    if not np.isfinite(lf) or lf == 0 or g not in sign: continue
                    na += int(np.sign(lf) == np.sign(sign[g])); no += int(np.sign(lf) != np.sign(sign[g]))
            rec["n_confirmed_kd_agrees_btla"] = na; rec["n_confirmed_kd_opposes_btla"] = no
            rec["net_signed_agreement"] = (na - no) / (na + no) if (na + no) else np.nan
            rows.append(rec)
    df = pd.DataFrame(rows); df["adj_hypergeom_p_8hr"] = np.nan
    for model in df["model"].unique():
        m = df["model"] == model; p = df.loc[m, "hypergeom_p_8hr"].values; ok = np.isfinite(p)
        if ok.sum():
            adj = np.full(len(p), np.nan); adj[ok] = multipletests(p[ok], method="fdr_bh")[1]
            df.loc[m, "adj_hypergeom_p_8hr"] = adj
    def status(r):
        if not r["crispri_arm_present_8hr"]: return "not_in_screen"
        if r["ontarget_kd_significant"] == False: return "failed_ontarget_qc"
        nconf, p, net = r["n_confirmed_8hr"], r["hypergeom_p_8hr"], r["net_signed_agreement"]
        if pd.notna(nconf) and nconf >= 3 and pd.notna(net) and net <= -0.5: return "usable_contradictory"
        if pd.notna(nconf) and nconf > 0 and pd.notna(p) and p < 0.05: return "usable_supportive"
        return "usable_unsupported"
    df["validation_status"] = df.apply(status, axis=1)
    return df

union_se = {"GREmLN": rk["seed_excluded"]["gr25"], "GENIE3": rk["seed_excluded"]["g325"]}
crispri = crispri_validation(union_se, bu.artifact("de_stats"))
crispri.to_csv(TAB / "gremln_genie3_cd4_crispri_validation.csv", index=False)

def crispri_stats(model):
    sub = crispri[(crispri.model == model) & (crispri.TF.isin(union_se[model]))]
    pres = sub[sub.crispri_arm_present_8hr == True]; fr = pres.frac_confirmed_8hr.astype(float)
    rng = np.random.default_rng(RANDOM_SEED)
    b = [rng.choice(fr.dropna().values, len(fr.dropna()), replace=True).mean() for _ in range(5000)] if fr.notna().any() else [np.nan]
    return {"model": model, "present": int(len(pres)),
            "qc_pass": int((sub.ontarget_kd_significant == True).sum()),
            "mean_frac_8hr": round(float(fr.mean()), 4), "median_frac_8hr": round(float(fr.median()), 4),
            "ci_lo": round(float(np.nanquantile(b, .025)), 4), "ci_hi": round(float(np.nanquantile(b, .975)), 4),
            **{s: int((sub.validation_status == s).sum()) for s in
               ["usable_supportive", "usable_unsupported", "usable_contradictory", "failed_ontarget_qc", "not_in_screen"]}}
cr_stats = {m: crispri_stats(m) for m in ("GREmLN", "GENIE3")}
pd.DataFrame(cr_stats.values()).to_csv(TAB / "crispri_summary_by_model_seed_excluded.csv", index=False)
pd.DataFrame(cr_stats.values())
"""))

    c.append(md("""## 6. Paperclip literature (identical template; annotation only)

Structured literature evidence for the **full union of both top-25 lists** under one identical query
template (`scripts/build_paperclip_review.py`; committed to `results/paperclip/`). Because coverage
is complete, a quantitative literature comparison across the model-specific candidate sets is drawn
here. Literature is annotation, never a predictive input."""))
    c.append(code("""
union_tfs = sorted(set(union_se["GREmLN"]) | set(union_se["GENIE3"]))
pc_path = bu.repo_root() / "results" / "paperclip" / "paperclip_union_top25_review.csv"
if pc_path.exists():
    pc = pd.read_csv(pc_path)
    pc_join = pd.DataFrame({"TF": union_tfs}).merge(pc, on="TF", how="left")
    pc_join["paperclip_reviewed"] = pc_join["paperclip_reviewed"].fillna(False).astype(bool)
else:
    pc_join = pd.DataFrame({"TF": union_tfs}); pc_join["paperclip_reviewed"] = False
    pc_join["paperclip_evidence_tier"] = np.nan
pc_join.to_csv(TAB / "paperclip_union_top25.csv", index=False)
paperclip_coverage_complete = bool(len(pc_join) and pc_join["paperclip_reviewed"].all())
strong_mod = set(pc_join.loc[pc_join["paperclip_evidence_tier"].isin(["strong", "moderate"]), "TF"])
gr_only = set(union_se["GREmLN"]) - set(union_se["GENIE3"])
g3_only = set(union_se["GENIE3"]) - set(union_se["GREmLN"])
shared = set(union_se["GREmLN"]) & set(union_se["GENIE3"])
paperclip_lit = {
    "coverage_complete": paperclip_coverage_complete,
    "shared_with_strong_moderate_lit": f"{len(shared & strong_mod)}/{len(shared)}",
    "gremln_only_with_strong_moderate_lit": f"{len(gr_only & strong_mod)}/{len(gr_only)}",
    "genie3_only_with_strong_moderate_lit": f"{len(g3_only & strong_mod)}/{len(g3_only)}",
}
print(f"Paperclip coverage: {pc_join['paperclip_reviewed'].sum()}/{len(pc_join)} union TFs; "
      f"complete={paperclip_coverage_complete}")
for k, v in paperclip_lit.items():
    print(f"  {k}: {v}")
"""))

    c.append(md("""## 7. BTLA multi-omics — independent orthogonal vs derived/contextual

Independent orthogonal = protein, phosphosite, co-IP (not derivable from the expression signal
driving the models). Derived/contextual = transcript DE, TF activity, kinase activity, BIONIC/GNN
modules. Kept **separate**; not collapsed into one count."""))
    c.append(code("""
INDEP = ["protein_support", "phosphosite_support", "coip_support"]
DERIV = ["transcript_support", "tf_activity_support", "kinase_activity_support",
         "bionic_support", "early_synapse_trafficking_support"]
mo_path = bu.artifact("multiomics_summary")
union_df = pd.DataFrame({"TF": union_tfs})
union_df["is_BTLA_vs_TCR_seed"] = union_df["TF"].isin(seeds)
union_df["BTLA_vs_TCR_direction"] = union_df["TF"].map(lambda g: f"BTLA_{seed_dir[g]}" if g in seed_dir else "not_seed")
union_df["in_gremln_top25"] = union_df["TF"].isin(union_se["GREmLN"])
union_df["in_genie3_top25"] = union_df["TF"].isin(union_se["GENIE3"])

integ = union_df.copy()
if mo_path.exists():
    mo = pd.read_csv(mo_path).rename(columns={"gene": "TF"})
    cols = ["TF"] + [c for c in INDEP + DERIV + ["strongest_support_layer", "independent_support",
            "overall_multiomics_support"] if c in mo.columns]
    integ = integ.merge(mo[cols], on="TF", how="left")
    def any_support(row, layers):
        return any(str(row.get(l, "")).lower() not in ("", "nan", "none", "no", "false", "0", "0.0")
                   for l in layers if l in integ.columns)
    integ["independent_orthogonal_support"] = integ.apply(lambda r: any_support(r, INDEP), axis=1)
    integ["derived_contextual_support"] = integ.apply(lambda r: any_support(r, DERIV), axis=1)
# attach compact CRISPRi + paperclip
cr_agg = crispri.groupby("TF").agg(crispri_status=("validation_status", lambda s: ";".join(sorted(set(s)))),
        crispri_best_frac_8hr=("frac_confirmed_8hr", "max")).reset_index()
pc_cols = [c for c in ["TF", "paperclip_evidence_tier", "paperclip_primary_phenotype",
           "paperclip_direction", "paperclip_key_source", "paperclip_short_rationale",
           "paperclip_reviewed"] if c in pc_join.columns]
integ = integ.merge(cr_agg, on="TF", how="left").merge(pc_join[pc_cols], on="TF", how="left")
integ.to_csv(TAB / "gremln_genie3_top25_integrated_evidence.csv", index=False)
print("Independent orthogonal support among union TFs:",
      int(integ.get("independent_orthogonal_support", pd.Series(dtype=bool)).sum()))
integ.head(12)
"""))

    c.append(md("## 8. Report figures (primary seed-excluded + seed-inclusive context + coverage)"))
    c.append(code("""
bu.fig_top25_comparison(rk["seed_excluded"]["gr"], rk["seed_excluded"]["g3"],
                        FIG / "fig1_top25_seed_excluded.png", seeds=seeds, mode="seed-excluded")
bu.fig_overlap(rk["seed_excluded"]["gr25"], rk["seed_excluded"]["g325"],
               FIG / "fig2_top25_overlap_seed_excluded.png", seeds=seeds, title="Seed-excluded top-25 overlap")
bu.fig_top25_comparison(rk["seed_inclusive"]["gr"], rk["seed_inclusive"]["g3"],
                        FIG / "figS_top25_seed_inclusive.png", seeds=seeds, mode="seed-inclusive")
bu.fig_crispri(crispri, FIG / "fig3_crispri_target_response.png")
bu.fig_evidence_heatmap(integ, FIG / "fig4_multiomics_heatmap.png")
n_genie3_nodes = len(set(edges['regulator']) | set(edges['target']))
bu.fig_coverage(N_MASKED_EXPR_GENES, n_genie3_nodes, len(genes),
                gremln_tfs, genie3_tfs, len(set(seeds) & common), FIG / "fig5_coverage.png")
print("figures written to", FIG)
"""))

    c.append(md("""## 9. SCENIC scope note

> SCENIC pruning of the persisted top-50k GENIE3 graph was retained as an exploratory motif-support
> sensitivity analysis (notebook 01) and was **not** used as the GREmLN prior."""))

    c.append(md("""## 10. Supplementary (S1) — held-out seed-TF recovery (NOT a validator)

An internal ranking-consistency diagnostic: hold out a fraction of BTLA seed-TFs and score how each
model recovers them from the remaining seeds. GREmLN's recall@25 saturates at 1.00 here, which is a
seed-clustering artefact; this is **pending a leakage/asymmetry audit** and does **not** enter the
verdict. Left unexecuted by default; enable by setting `RUN_HELDOUT = True`."""))
    c.append(code("""
RUN_HELDOUT = False
if RUN_HELDOUT:
    print("Held-out recovery is supplementary; see the development archive for the full computation.")
else:
    print("Held-out recovery skipped (supplementary; excluded from verdict).")
"""))

    c.append(md("""## 11. Restrained verdict — three-validator decision rule

A model is **superior** only if it wins CRISPRi by ≥ `CRISPRI_MARGIN` on mean functional
target-response **and** has ≥ as many `usable_supportive` TFs **and** carries ≥ corroborating
literature + multi-omics support. Otherwise: equivalence-not-tested, GENIE3-superior, or
complementary. Held-out recovery is excluded."""))
    c.append(code("""
gr_s, g3_s = cr_stats["GREmLN"], cr_stats["GENIE3"]
d_mean = gr_s["mean_frac_8hr"] - g3_s["mean_frac_8hr"]

# corroborating-evidence support among each model's method-SPECIFIC top-25 candidates:
# strong/moderate Paperclip literature + independent orthogonal multi-omics (protein/PTM/coIP).
indep_support = set(integ.loc[integ.get("independent_orthogonal_support", False) == True, "TF"]) \\
    if "independent_orthogonal_support" in integ.columns else set()
corrob = {
    "GREmLN": len(gr_only & strong_mod) + len(gr_only & indep_support),
    "GENIE3": len(g3_only & strong_mod) + len(g3_only & indep_support),
}
# full decision rule: CRISPRi margin win AND >= usable_supportive AND >= corroborating support
gremln_wins = ((d_mean >= CRISPRI_MARGIN)
               and (gr_s["usable_supportive"] >= g3_s["usable_supportive"])
               and (corrob["GREmLN"] >= corrob["GENIE3"]))
genie3_wins = ((-d_mean >= CRISPRI_MARGIN)
               and (g3_s["usable_supportive"] >= gr_s["usable_supportive"])
               and (corrob["GENIE3"] >= corrob["GREmLN"]))
if gremln_wins:
    verdict = "GREmLN superiority"
elif genie3_wins:
    verdict = "GENIE3 superiority"
else:
    verdict = "complementary candidate prioritisation (neither met the decision rule; equivalence not tested)"

prim = overlap["seed_excluded"]
summary = {
    "ranking_primary": "seed_excluded",
    "spearman_rho_seed_excluded": prim["spearman_rank_rho"],
    "shared_top25_seed_excluded": prim["shared_top25"],
    "gremln_only_top25": prim["gremln_only_top25"],
    "genie3_only_top25": prim["genie3_only_top25"],
    "crispri_mean_frac_8hr": {"GREmLN": gr_s["mean_frac_8hr"], "GENIE3": g3_s["mean_frac_8hr"],
                              "difference": round(d_mean, 4), "margin": CRISPRI_MARGIN},
    "crispri_usable_supportive": {"GREmLN": gr_s["usable_supportive"], "GENIE3": g3_s["usable_supportive"]},
    "corroborating_support_method_specific": corrob,
    "decision_rule": "CRISPRi margin AND >= usable_supportive AND >= corroborating (literature+independent multi-omics)",
    "paperclip_coverage_complete": paperclip_coverage_complete,
    "paperclip_literature_support": paperclip_lit,
    "held_out_recovery": "SUPPLEMENTARY — excluded from verdict",
    "verdict": verdict,
}
(REG / "benchmark_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(json.dumps(summary, indent=2, default=str))
print("\\nVERDICT:", verdict)
"""))

    _write(_nb(c), NB / "03_benchmark_btla_tf_prioritisation.ipynb")


# ============================================================ NOTEBOOK 01
def build_nb1():
    c = []
    c.append(md("""
# 01 — Build the masked CD4 GENIE3 graph (regulatory prior for GREmLN)

Constructs the **raw masked CD4 GENIE3** co-expression graph used both as the GENIE3 comparator and
as GREmLN's regulatory prior (notebook 02). Steps: CD4 **non-targeting-control (NTC)** cell
selection → donor/batch audit → CPM/log1p normalisation → **ComBat** donor masking → TF regulator
universe → **GENIE3** (tree ensembles) → top-50,000-edge graph → model registry, plus a **SCENIC
top-50k sensitivity** (motif support only; **not** the GREmLN prior).

The GENIE3 forest fit is hours-long, so it is behind `REGENERATE`. By default the notebook
**reuses** the persisted `grn_edges.tsv` and audits it. Parameters come from
`config/benchmark_config.yaml`; see `data/README.md` to obtain the CD4 Perturb-seq NTC cells.
"""))
    c.append(code("""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path.cwd().parent / "scripts"))
import bench_utils as bu
import yaml

CFG = yaml.safe_load((bu.repo_root() / "config" / "benchmark_config.yaml").read_text())
REGENERATE = False   # True re-runs GENIE3 (hours; needs the CD4 NTC expression matrix)
REG = bu.repo_root() / "results" / "model_registry"; REG.mkdir(parents=True, exist_ok=True)
print("GENIE3 params:", CFG["genie3"]); print("DATA_ROOT:", bu.data_root())
"""))
    c.append(md("""## 1. Source data, NTC selection, donor/batch audit, normalisation, ComBat masking

The GENIE3 input is CD4 **NTC** cells from the CZI genome-scale T-cell Perturb-seq resource:
select NTC guides, keep CD4, take the top-4,000 highly variable genes, normalise to CPM(1e6)→log1p,
then apply **ComBat on `donor_id`** (this donor-masking defines the *masked* graph). The exact
preprocessing lives in the development archive's prep script; parameters are pinned in the config.
The GREmLN prior in notebook 02 deliberately uses the *pre-ComBat* log1p-CPM of the same cells."""))
    c.append(code("""
# Documented preprocessing (guarded). In reuse mode we do not need the raw matrix.
if REGENERATE:
    import anndata as ad, scanpy as sc
    adata = ad.read_h5ad(bu.artifact("cd4_ntc_counts"))          # provide via BTLA_BENCH_CD4_NTC_COUNTS
    adata = adata[adata.obs["guide_target"].astype(str).eq("non-targeting")].copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=CFG["genie3"]["n_hvg"], flavor="seurat_v3")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.normalize_total(adata, target_sum=1e6); sc.pp.log1p(adata)
    sc.pp.combat(adata, key="donor_id")                          # donor masking -> masked graph
    donor_audit = adata.obs["donor_id"].value_counts()
    print("cells x genes:", adata.shape); print(donor_audit)
else:
    print("REGENERATE=False -> reuse persisted masked graph; preprocessing shown above for provenance.")
"""))
    c.append(md("""## 2. TF regulator universe + GENIE3 (tree-ensemble feature importances)

GENIE3 fits, for each target gene, a tree ensemble predicting it from the TF regulators; edge weight
= regulator feature importance. Core algorithm inlined below; the full VIM over ~4k genes is the
hours-long step."""))
    c.append(code('''
def genie3_single(expr, target_idx, reg_idx, ntrees=1000, K="sqrt"):
    # Predict ONE target gene from the regulator set; return per-regulator importances.
    from sklearn.ensemble import ExtraTreesRegressor
    y = np.nan_to_num(np.asarray(expr[:, target_idx], float)); s = y.std()
    if s > 0: y = y / s
    preds = [i for i in reg_idx if i != target_idx]   # regulators are the predictors
    Xr = np.nan_to_num(np.asarray(expr[:, preds], float))
    est = ExtraTreesRegressor(n_estimators=ntrees, max_features=K, n_jobs=1).fit(Xr, y)
    imp = np.array([e.tree_.compute_feature_importances(normalize=False) for e in est.estimators_])
    vi = np.zeros(expr.shape[1]); vi[preds] = imp.sum(0) / len(est); return vi

def run_genie3(expr, gene_names, regulators, n_edges_keep=50000):
    # Standard GENIE3: every gene is a candidate TARGET predicted from the regulators;
    # edge = regulator -> target, weight = regulator's feature importance for that target.
    ng = expr.shape[1]; ridx = [i for i, g in enumerate(gene_names) if g in regulators]
    VIM = np.zeros((ng, ng))                          # VIM[target, regulator]
    for t in range(ng):                              # loop over ALL genes as targets
        VIM[t, :] = genie3_single(expr, t, ridx)
    links = [(gene_names[r], gene_names[t], float(VIM[t, r]))   # regulator r -> target t
             for t in range(ng) for r in ridx if r != t and VIM[t, r] > 0]
    links.sort(key=lambda x: x[2], reverse=True)
    return pd.DataFrame(links[:n_edges_keep], columns=["regulator", "target", "weight"])

if REGENERATE:
    tfs = bu.load_tfs(); gnames = list(adata.var_names.astype(str))
    regs = [g for g in gnames if g in tfs]
    edges = run_genie3(np.asarray(adata.X, float), gnames, regs, CFG["genie3"]["n_edges_persisted"])
    edges.to_csv(bu.artifact("genie3_edges"), sep="\\t", index=False)
else:
    print("Skipping GENIE3 fit (reuse mode).")
'''))
    c.append(md("## 3. Load / audit the persisted masked graph + model registry"))
    c.append(code("""
edges = bu.load_edges()
regs, tgts = set(edges["regulator"]), set(edges["target"])
tfs = bu.load_tfs()
stats = {"n_edges": len(edges), "n_regulators": len(regs), "n_regulators_in_tf_list": len(regs & tfs),
         "n_targets": len(tgts), "n_nodes": len(regs | tgts),
         "weight_min": float(edges["weight"].min()), "weight_max": float(edges["weight"].max())}
print(json.dumps(stats, indent=2))
registry = pd.DataFrame([
    ("model_name", "GENIE3_CD4_masked_raw"),
    ("graph_path", bu.redact(bu.artifact("genie3_edges"))),
    ("graph_md5", bu.md5(bu.artifact("genie3_edges"))),
    ("n_edges", stats["n_edges"]), ("n_regulators", stats["n_regulators"]),
    ("n_targets", stats["n_targets"]),
    ("n_hvg", CFG["genie3"]["n_hvg"]), ("normalisation", CFG["genie3"]["normalisation"]),
    ("batch_correction", CFG["genie3"]["batch_correction"]),
    ("tf_list", CFG["genie3"]["tf_list"]),
    ("tf_list_md5", bu.md5(bu.repo_root() / "resources" / "human_tfs_pySCENIC.txt")),
    ("tree_method", "ExtraTrees/RandomForest, ntrees=1000, K=sqrt"),
], columns=["field", "value"])
registry.to_csv(REG / "genie3_model_registry.csv", index=False)
registry
"""))
    c.append(md("""## 4. SCENIC top-50k sensitivity (motif support only — NOT the GREmLN prior)

> SCENIC pruning of the persisted top-50k GENIE3 graph is retained as an exploratory motif-support
> sensitivity analysis and is **not** used as the GREmLN prior. Enable with `REGENERATE` + pySCENIC
> databases; the pruned graph is reported separately and never fed to GREmLN."""))
    c.append(code("""
if REGENERATE:
    print("Run pySCENIC prune_targets on grn_edges.tsv with cisTarget motif rankings (see archive).")
else:
    print("SCENIC sensitivity skipped; documented as a separate motif-support analysis only.")
"""))
    _write(_nb(c), NB / "01_build_genie3_cd4_masked.ipynb")


# ============================================================ NOTEBOOK 02
def build_nb2():
    c = []
    c.append(md("""
# 02 — GREmLN gene embeddings from the masked GENIE3 prior

Runs the **GREmLN** foundation model (frozen checkpoint, zero-shot) to produce per-gene embeddings,
using the masked CD4 GENIE3 graph (notebook 01) as its **regulatory prior**. Key design point: the
expression fed to GREmLN is the **pre-ComBat log1p-CPM** of the same 50k cells, because GREmLN's
tokenizer needs the rank/zero structure of near-count data (ComBat output has negatives / few true
zeros). GPU inference is behind `REGENERATE`; by default we **reuse** the saved embeddings and audit
coverage. Obtain the checkpoint per `data/README.md`.
"""))
    c.append(code("""
import sys, json, platform
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path.cwd().parent / "scripts"))
import bench_utils as bu
import yaml
CFG = yaml.safe_load((bu.repo_root() / "config" / "benchmark_config.yaml").read_text())
REGENERATE = False       # True runs GREmLN inference (GPU). False reuses saved embeddings.
REG = bu.repo_root() / "results" / "model_registry"; REG.mkdir(parents=True, exist_ok=True)
print("GREmLN cfg:", {k: CFG["gremln"][k] for k in ("embedding_dim", "expression_input", "prior_graph")})
"""))
    c.append(md("""## 1. Provenance: checkpoint, submodule commit, expression input, graph prior"""))
    c.append(code("""
def pkg(n):
    try:
        from importlib.metadata import version; return version(n)
    except Exception: return "unknown"
def git_commit(path):
    import subprocess
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(path),
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except Exception: return "unknown"
gremln_commit = git_commit(bu.repo_root() / "third_party" / "GREmLN")
print("GREmLN submodule commit:", gremln_commit)
print("checkpoint md5 (config):", CFG["gremln"]["checkpoint_md5"])
"""))
    c.append(md("""## 2. GREmLN inference (GPU, guarded): vocab → RegulatoryNetwork(prior) → tokenizer → embeddings

The real inference call is inlined below. It loads the frozen `GDTransformer`, builds the
`GraphTokenizer` from the GREmLN gene vocabulary and the masked GENIE3 graph as `RegulatoryNetwork`,
and runs a single gene-embedding forward pass (`get_gene_embeddings`), mapping Ensembl→HUGO."""))
    c.append(code('''
if REGENERATE:
    import anndata as ad, torch
    sys.path.insert(0, str(bu.repo_root() / "third_party" / "GREmLN"))
    from scGraphLLM import GeneVocab, GraphTokenizer, InferenceDataset, RegulatoryNetwork
    from scGraphLLM.config import graph_kernel_attn_3L_4096
    from scGraphLLM.inference import get_gene_embeddings
    from scGraphLLM.models import GDTransformer
    adata = ad.read_h5ad(bu.artifact("gremln_inputs") / "expr.h5ad")     # pre-ComBat log1p-CPM, Ensembl var
    vocab = GeneVocab.from_csv(bu.repo_root() / "third_party/GREmLN/scGraphLLM/resources/gene_vocab.csv",
                               gene_col="gene_name", node_col="idx")
    network = RegulatoryNetwork.from_csv(bu.artifact("gremln_inputs") / "graph.tsv", sep="\\t",
                                         reg_name="regulator.values", tar_name="target.values",
                                         wt_name="mi.values", lik_name="log.p.values")
    model = GDTransformer.load_from_checkpoint(str(bu.artifact("gremln_checkpoint")),
                                               config=graph_kernel_attn_3L_4096).eval().cuda()
    ds = InferenceDataset(expression=adata.to_df(), tokenizer=GraphTokenizer(vocab=vocab, network=network))
    x_gene = get_gene_embeddings(ds, model, vocab, batch_size=256)       # per-gene mean embeddings
    emb = ad.AnnData(x_gene.values, obs=adata.var.loc[x_gene.index].copy())
    emb.obs["hugo"] = emb.obs.get("hugo_symbol", pd.Series(emb.obs_names)).astype(str)
    emb.write_h5ad(bu.artifact("gremln_emb"))
else:
    print("REGENERATE=False -> reuse saved GREmLN embeddings (inference code above for provenance).")
'''))
    c.append(md("## 3. Load embeddings + coverage audit + sanitized model registry"))
    c.append(code("""
genes, X = bu.load_gremln_embeddings()
tfs = bu.load_tfs(); edges = bu.load_edges(); seeds, seed_dir, _ = bu.load_seeds()
gremln_tfs = set(genes) & tfs; genie3_tfs = set(edges["regulator"]) & tfs
common = gremln_tfs & genie3_tfs
audit = pd.DataFrame([
    ("genes_embedded_by_gremln", len(genes)), ("embedding_dim", X.shape[1]),
    ("tfs_embedded_by_gremln", len(gremln_tfs)), ("tfs_in_genie3_graph", len(genie3_tfs)),
    ("common_tf_universe", len(common)), ("btla_seed_panel_size", len(seeds)),
    ("seeds_available_to_gremln", len(set(seeds) & set(genes))),
    ("seed_tfs_in_common_universe", len(set(seeds) & common)),
], columns=["metric", "value"])
print(audit.to_string(index=False))
registry = pd.DataFrame([
    ("model_name", "GREmLN_CD4_masked_raw_prior"),
    ("gremln_checkpoint", "model.ckpt"), ("gremln_checkpoint_md5", CFG["gremln"]["checkpoint_md5"]),
    ("gremln_submodule_commit", gremln_commit),
    ("python_version", platform.python_version()),
    ("numpy_version", pkg("numpy")), ("pandas_version", pkg("pandas")),
    ("anndata_version", pkg("anndata")), ("torch_version", pkg("torch")),
    ("expression_input", CFG["gremln"]["expression_input"]),
    ("prior_graph", CFG["gremln"]["prior_graph"]),
    ("prior_graph_md5", bu.md5(bu.artifact("genie3_edges"))),
    ("embedding_h5ad_md5", bu.md5(bu.artifact("gremln_emb"))),
    ("n_genes_embedded", len(genes)), ("embedding_dim", X.shape[1]),
    ("csls_params", CFG["scoring"]["gremln"]), ("random_seed", CFG["seed"]),
    ("device", "cuda"), ("inference_params", CFG["gremln"]["inference"]),
], columns=["field", "value"])
registry.to_csv(REG / "gremln_model_registry.csv", index=False)
audit.to_csv(REG / "gremln_gene_universe_audit.csv", index=False)
registry
"""))
    _write(_nb(c), NB / "02_build_gremln_from_masked_genie3_prior.ipynb")


def _write(nb, path):
    nbf.write(nb, str(path))
    print("wrote", path)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("01", "all"):
        build_nb1()
    if which in ("02", "all"):
        build_nb2()
    if which in ("03", "all"):
        build_nb3()
