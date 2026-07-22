#!/usr/bin/env python3
"""Sections 4-5: protein/phospho threshold sensitivity + source-version reconciliation.

Evidence checks only; no manuscript / verdict changes. Restricted numeric values
are never written to committed outputs — only statuses, counts, hashes and
schema.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

PADJ_MAX = 0.05
PROT_LOG2 = math.log2(1.3)      # 0.379 source-native protein/phospho
HARM_LOG2 = math.log2(1.5)      # 0.585 harmonised sensitivity
MS_TPS = ["2mn", "20mn", "4h", "24h"]
ACUTE = {"2mn", "20mn"}

GM = Path("/mnt/R0/Projects/POIAZ/gmorgan/gremln-tcells/data/multiomics")
VT = Path("/mnt/R0/Projects/POIAZ/Viet/BulkFormer/BTLA_CSLS_multiomics/data/multiomics")

GM_FILES = {
    "transcriptomics": GM / "Transcriptomics/@Transcriptomics/Tables/DEG_wide_statistics.xlsx",
    "proteomics": GM / "Proteomics/@Proteomics/Tables/WP_DEP_wide_statistics.xlsx",
    "phosphoproteomics": GM / "Phosphoproteomics/@Phosphoproteomics/@Tables/DEPPs_All_statistics.xlsx",
}
VT_FILES = {
    "transcriptomics": VT / "DEG_wide_statistics.xlsx",
    "proteomics": VT / "WP_DEP_wide_statistics.xlsx",
    "phosphoproteomics": VT / "DEPPs_All_statistics.xlsx",
}


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for c in [here.parent, *here.parents]:
        if (c / "scripts").is_dir() and (c / "results").is_dir():
            return c
    return here.parent.parent


def sha256_file(p: Path) -> str:
    if not p.exists():
        return "missing"
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_union(repo):
    u = pd.read_csv(repo / "results/publication_data/top25_union_primary.csv")

    def grp(r):
        if r["in_gremln_top25"] and r["in_genie3_top25"]:
            return "shared"
        return "GREmLN_specific" if r["in_gremln_top25"] else "GENIE3_specific"
    u["candidate_group"] = u.apply(grp, axis=1)
    return u


# --------------------------------------------------------------------------- #
def protein_qualifying(prot_path, cand, log2_min):
    """Return dict tf -> list of qualifying timepoints (Btla contrast)."""
    prot = pd.read_excel(prot_path)
    gene_col = "Gene" if "Gene" in prot.columns else "gene"
    by = {str(g): sub for g, sub in prot.groupby(gene_col)}
    res = {}
    universe = set(prot[gene_col].astype(str))
    for tf in cand:
        if tf not in by:
            continue
        r = by[tf].iloc[0]
        tps = []
        for tp in MS_TPS:
            eff, padj = r.get(f"Btla_{tp}_logFC"), r.get(f"Btla_{tp}_adj.P.Val")
            if pd.notna(eff) and pd.notna(padj) and padj < PADJ_MAX and abs(eff) > log2_min:
                tps.append(tp)
        if tps:
            res[tf] = tps
    return res, universe


def phospho_qualifying(phos_path, cand, log2_min):
    """Return dict tf -> list of (site,timepoint) qualifying (Btla sheets)."""
    xls = pd.ExcelFile(phos_path)
    frames = {s.replace("Btla_", ""): pd.read_excel(phos_path, sheet_name=s)
              for s in xls.sheet_names if s.startswith("Btla_")}
    universe = set()
    for df in frames.values():
        universe |= set(df["Gene"].astype(str))
    res = {}
    for tf in cand:
        hits = []
        for tp, df in frames.items():
            sub = df[df["Gene"].astype(str) == tf]
            for _, s in sub.iterrows():
                eff, padj = s.get("logFC"), s.get("adj.P.Val")
                if pd.notna(eff) and pd.notna(padj) and padj < PADJ_MAX and abs(eff) > log2_min:
                    hits.append((str(s.get("Gene_STY_canonical", "")), tp))
        if hits:
            res[tf] = hits
    return res, universe


def section4(repo, union):
    cand = list(union["TF"])
    grp = dict(zip(union["TF"], union["candidate_group"]))
    prot_p = GM_FILES["proteomics"]
    phos_p = GM_FILES["phosphoproteomics"]

    prot_prim, _ = protein_qualifying(prot_p, cand, PROT_LOG2)
    prot_harm, _ = protein_qualifying(prot_p, cand, HARM_LOG2)
    phos_prim, _ = phospho_qualifying(phos_p, cand, PROT_LOG2)
    phos_harm, _ = phospho_qualifying(phos_p, cand, HARM_LOG2)

    def assoc(protd, phosd):
        return set(protd) | set(phosd)

    prim = assoc(prot_prim, phos_prim)
    harm = assoc(prot_harm, phos_harm)

    rows = []
    for tf in sorted(prim | harm):
        in_prim = tf in prim
        in_harm = tf in harm
        # layer+timepoint detail (site residues NOT emitted; count + timepoint only)
        prot_tps_p = ";".join(prot_prim.get(tf, []))
        prot_tps_h = ";".join(prot_harm.get(tf, []))
        phos_tps_p = ";".join(sorted(set(tp for _, tp in phos_prim.get(tf, []))))
        phos_tps_h = ";".join(sorted(set(tp for _, tp in phos_harm.get(tf, []))))
        rows.append({
            "TF": tf, "candidate_group": grp.get(tf, ""),
            "primary_association": in_prim, "harmonised_association": in_harm,
            "change": ("retained" if in_prim and in_harm else
                       "lost_under_harmonised" if in_prim and not in_harm else
                       "gained_under_harmonised"),
            "protein_qual_tps_primary": prot_tps_p,
            "protein_qual_tps_harmonised": prot_tps_h,
            "phospho_qual_tps_primary": phos_tps_p,
            "phospho_qual_tps_harmonised": phos_tps_h,
            "n_phospho_sites_primary": len(phos_prim.get(tf, [])),
            "n_phospho_sites_harmonised": len(phos_harm.get(tf, [])),
        })
    sens = pd.DataFrame(rows)

    def group_counts(s):
        return {g: sorted(union[(union["candidate_group"] == g) &
                                (union["TF"].isin(s))]["TF"].tolist())
                for g in ["GREmLN_specific", "shared", "GENIE3_specific"]}

    summary = {
        "primary_threshold": "protein/phospho padj<0.05 & |log2FC|>log2(1.3)=%.3f" % PROT_LOG2,
        "harmonised_threshold": "padj<0.05 & |log2FC|>=log2(1.5)=%.3f" % HARM_LOG2,
        "primary_association_candidates": sorted(prim),
        "harmonised_association_candidates": sorted(harm),
        "retained": sorted(prim & harm),
        "lost_under_harmonised": sorted(prim - harm),
        "gained_under_harmonised": sorted(harm - prim),
        "primary_group_counts": {k: len(v) for k, v in group_counts(prim).items()},
        "harmonised_group_counts": {k: len(v) for k, v in group_counts(harm).items()},
        "primary_groups": group_counts(prim),
        "harmonised_groups": group_counts(harm),
    }
    return sens, summary


# --------------------------------------------------------------------------- #
def schema_of(path, layer):
    if not path.exists():
        return {"exists": False}
    xls = pd.ExcelFile(path)
    sheets = xls.sheet_names
    first = pd.read_excel(path, sheet_name=sheets[0], nrows=5)
    info = {
        "exists": True, "sha256": sha256_file(path), "n_sheets": len(sheets),
        "sheets": sheets[:12], "n_cols_sheet0": int(first.shape[1]),
        "columns_sheet0": list(map(str, first.columns)),
    }
    return info


def section5(repo, union):
    cand = list(union["TF"])
    out = {}
    for layer in ["transcriptomics", "proteomics", "phosphoproteomics"]:
        gm, vt = GM_FILES[layer], VT_FILES[layer]
        gm_s, vt_s = schema_of(gm, layer), schema_of(vt, layer)
        entry = {
            "gmorgan": {"filename": gm.name, **gm_s},
            "viet": {"filename": vt.name, **vt_s},
            "sha256_match": gm_s.get("sha256") == vt_s.get("sha256"),
            "schema_match": gm_s.get("columns_sheet0") == vt_s.get("columns_sheet0")
                            and gm_s.get("sheets") == vt_s.get("sheets"),
        }
        # candidate coverage + qualifying diff
        if layer == "proteomics" and vt_s.get("exists"):
            gm_q, gm_u = protein_qualifying(gm, cand, PROT_LOG2)
            vt_q, vt_u = protein_qualifying(vt, cand, PROT_LOG2)
            entry["coverage"] = {"gmorgan_in_universe": len(set(cand) & gm_u),
                                 "viet_in_universe": len(set(cand) & vt_u)}
            entry["qualifying_diff"] = {"gmorgan_only": sorted(set(gm_q) - set(vt_q)),
                                        "viet_only": sorted(set(vt_q) - set(gm_q)),
                                        "both": sorted(set(gm_q) & set(vt_q))}
        elif layer == "phosphoproteomics" and vt_s.get("exists"):
            gm_q, gm_u = phospho_qualifying(gm, cand, PROT_LOG2)
            vt_q, vt_u = phospho_qualifying(vt, cand, PROT_LOG2)
            entry["coverage"] = {"gmorgan_in_universe": len(set(cand) & gm_u),
                                 "viet_in_universe": len(set(cand) & vt_u)}
            entry["qualifying_diff"] = {"gmorgan_only": sorted(set(gm_q) - set(vt_q)),
                                        "viet_only": sorted(set(vt_q) - set(gm_q)),
                                        "both": sorted(set(gm_q) & set(vt_q))}
        else:  # transcriptomics: coverage only
            if vt_s.get("exists"):
                gm_df = pd.read_excel(gm); vt_df = pd.read_excel(vt)
                gm_u = set(gm_df["gene"].astype(str)); vt_u = set(vt_df["gene"].astype(str))
                entry["coverage"] = {"gmorgan_in_universe": len(set(cand) & gm_u),
                                     "viet_in_universe": len(set(cand) & vt_u),
                                     "gmorgan_universe_size": len(gm_u),
                                     "viet_universe_size": len(vt_u)}
        out[layer] = entry
    return out


def main():
    repo = repo_root()
    outdir = repo / "results/multiomics/audit_v2"
    (outdir / "local_restricted").mkdir(parents=True, exist_ok=True)
    union = load_union(repo)

    sens, s4 = section4(repo, union)
    sens.to_csv(outdir / "protein_phospho_threshold_sensitivity.csv", index=False)
    s5 = section5(repo, union)

    report = {"section4_threshold_sensitivity": s4, "section5_source_version": s5}
    (outdir / "local_restricted" / "threshold_source_checks.json").write_text(
        json.dumps(report, indent=2, default=str))

    print("=== SECTION 4: threshold sensitivity (protein+phospho) ===")
    print("primary association:", s4["primary_association_candidates"])
    print("harmonised association:", s4["harmonised_association_candidates"])
    print("retained:", s4["retained"])
    print("LOST under harmonised:", s4["lost_under_harmonised"])
    print("gained under harmonised:", s4["gained_under_harmonised"])
    print("primary group counts:", s4["primary_group_counts"])
    print("harmonised group counts:", s4["harmonised_group_counts"])
    print("\n=== SECTION 5: source version ===")
    for layer, e in s5.items():
        print(f"\n[{layer}] sha256_match={e['sha256_match']} schema_match={e['schema_match']}")
        print("  gmorgan sha:", e["gmorgan"].get("sha256", "")[:16],
              "| viet sha:", e["viet"].get("sha256", "")[:16], "| viet exists:", e["viet"]["exists"])
        print("  coverage:", e.get("coverage"))
        if "qualifying_diff" in e:
            print("  qualifying diff:", e["qualifying_diff"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
