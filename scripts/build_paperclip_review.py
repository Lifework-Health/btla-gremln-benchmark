#!/usr/bin/env python3
"""Structured Paperclip literature review of the 50-TF union (both models x seed-excluded/inclusive
top-25), using one identical query template per TF plus a TF-focused second pass for canonical
immune TFs whose BTLA-scoped query returned only generic hits. Evidence tiers are assigned from the
retrieved full-text papers (public PMC). Output: results/paperclip/paperclip_union_top25_review.csv.

This is literature annotation (public papers), NOT a predictive input and NOT BTLA source data.
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import bench_utils as bu  # noqa

# tier, phenotype, tcr/checkpoint?, activation/exhaustion?, direction, key PMC, summary
R = {
 "AHR":     ("moderate","activation_differentiation","no","yes","context_dependent","PMC10089284","AhR regulates ILC/T-cell transcriptional programmes; context-dependent effector/regulatory roles."),
 "BACH2":   ("strong","exhaustion_differentiation","yes","yes","restrains_exhaustion","PMC7906956","BACH2 enforces stem-like CD8 fate and prevents terminal exhaustion."),
 "BCL3":    ("moderate","activation_function","no","no","context_dependent","PMC10552310","Bcl-3 regulates T-cell energy metabolism and proliferation via mTOR."),
 "BHLHE40": ("strong","checkpoint_effector","yes","yes","promotes_effector","PMC9164498","BHLHE40 is required for T-cell effector function and immune-checkpoint-therapy efficacy."),
 "CREM":    ("moderate","activation_differentiation","no","yes","restrains","PMC4770719","CREM is a negative regulator of Th2 responses; deficiency worsens allergic inflammation."),
 "CYCS":    ("none","none","no","no","na","","Cytochrome c; not a transcription factor / no T-cell TF role."),
 "EGR2":    ("strong","exhaustion","yes","yes","promotes_exhaustion","PMC8119420","EGR2 is induced by chronic antigen and required for exhausted CD8 stability/maintenance."),
 "ETV3L":   ("none","none","no","no","na","","No TF-specific T-cell literature retrieved."),
 "EZH2":    ("strong","differentiation_epigenetic","yes","yes","context_dependent","PMC7574680","EZH2 epigenetically regulates CD8 T-cell fate and function."),
 "FOXM1":   ("moderate","proliferation","no","no","promotes_proliferation","PMC2821927","FoxM1 is a master cell-cycle regulator for mature T cells."),
 "FOXP1":   ("strong","quiescence_exhaustion","yes","yes","context_dependent","PMC9576946","FOXP1 regulates T-cell quiescence, differentiation and exhaustion."),
 "GAR1":    ("none","none","no","no","na","","H/ACA ribonucleoprotein; not a T-cell TF."),
 "GTF2A2":  ("none","none","no","no","na","","General transcription factor IIA; housekeeping, no T-cell-specific role."),
 "GTF3C2":  ("none","none","no","no","na","","General transcription factor IIIC; housekeeping."),
 "HIRIP3":  ("none","none","no","no","na","","Histone-related; no T-cell TF literature retrieved."),
 "HIVEP3":  ("weak","general_tcell","no","no","context_dependent","PMC4255756","Appears only in general T-helper lineage reviews; no TF-specific T-cell paper retrieved."),
 "HMGN3":   ("none","none","no","no","na","","No HMGN3-specific T-cell literature (HMGB2 is the related exhaustion factor)."),
 "ID2":     ("strong","differentiation_exhaustion","yes","yes","promotes","PMC10902300","Id2 controls CD8 exhaustion (Tcf3-LSD1) and reinforces Th1 differentiation."),
 "IRF8":    ("moderate","differentiation","no","yes","restrains_Th17","PMC3112536","IRF8 directs a silencing programme restraining Th17 differentiation."),
 "JUNB":    ("strong","activation_treg","yes","yes","promotes_effector","PMC11884676","JunB is required for CD8 effector responses and effector-Treg homeostasis (BATF/IRF4)."),
 "KIF22":   ("none","none","no","no","na","","Kinesin motor protein; not a transcription factor."),
 "MAF":     ("strong","checkpoint_tolerance","yes","yes","context_dependent","PMC7033575","c-Maf is a multifaceted checkpoint/tolerance TF across T-cell subsets."),
 "MDM2":    ("moderate","treg","no","no","context_dependent","PMC7318079","MDM2 stabilises FOXP3 and modulates human regulatory T-cell function."),
 "MLX":     ("none","none","no","no","na","","Max-like protein X (metabolic bHLH); no T-cell TF role retrieved."),
 "MSC":     ("none","none","no","no","na","","No Musculin(MSC) T-cell TF literature retrieved (hit was mesenchymal stromal cells)."),
 "MYC":     ("moderate","activation_metabolism","no","yes","promotes_activation","PMC6938709","c-Myc is essential for Treg homeostasis and transitional activation (metabolism)."),
 "NFIL3":   ("strong","treg_differentiation","yes","yes","context_dependent","PMC6802641","NFIL3 controls Treg function/stability by downregulating Foxp3."),
 "NR4A2":   ("strong","tcr_induced_differentiation","yes","yes","promotes","PMC7883379","NR4A family is TCR-induced and shapes T-cell differentiation, exhaustion and Treg fate."),
 "NR4A3":   ("strong","tcr_induced_differentiation","yes","yes","promotes","PMC7883379","NR4A family is TCR-induced and shapes T-cell differentiation, exhaustion and Treg fate."),
 "PLAGL2":  ("none","none","no","no","na","","No T-cell TF literature retrieved."),
 "POU2AF1": ("moderate","differentiation","no","yes","restrains","PMC12069319","OCA-B/Pou2af1 restricts Th2 differentiation (via GATA3)."),
 "POU2F2":  ("weak","general_tcell","no","no","context_dependent","PMC7387970","Appears only in general Tfh differentiation reviews; no TF-specific paper retrieved."),
 "PRDM1":   ("strong","hyporesponsiveness_exhaustion","yes","yes","promotes_exhaustion","PMC9097451","PRDM1/Blimp-1 drives human T-cell hyporesponsiveness (transcriptome/epigenome)."),
 "RBPJ":    ("moderate","notch_signaling","no","yes","context_dependent","PMC9371899","RBPJ is the Notch effector TF (co-activator/repressor); Notch shapes T-cell fate."),
 "REL":     ("strong","activation_differentiation","yes","yes","promotes_effector","PMC3310234","c-Rel (NF-kB) is central to T-lymphocyte differentiation and effector functions."),
 "SMAP2":   ("none","none","no","no","na","","ArfGAP (SMAP2); not a transcription factor."),
 "SNAPC4":  ("none","none","no","no","na","","snRNA-activating complex; housekeeping, no T-cell TF role."),
 "SRP9":    ("none","none","no","no","na","","Signal recognition particle 9; not a transcription factor."),
 "STAT1":   ("strong","differentiation_th1","yes","yes","promotes_Th1","PMC8250105","STAT1 is crucial for Th1 differentiation and anti-tumour function."),
 "STAT4":   ("strong","differentiation_th1","yes","yes","promotes","PMC7097918","STAT4 is an immunoregulator central to Th1 and autoimmune/inflammatory disease."),
 "STAT5B":  ("moderate","function","no","yes","promotes","PMC3907443","STAT5B target genes in human T cells enrich for immune function/proliferation."),
 "STAU2":   ("none","none","no","no","na","","Staufen-2 RNA-binding protein; not a transcription factor."),
 "TBX21":   ("strong","differentiation_th1","yes","yes","promotes_Th1_effector","PMC4036734","T-bet (TBX21) is the master Th1/effector TF regulating cytokine production."),
 "TFDP1":   ("none","none","no","no","na","","DP-1 (E2F partner) cell-cycle factor; no immune-specific role retrieved."),
 "TGIF1":   ("none","none","no","no","na","","TGFB-induced homeobox 1; no T-cell TF literature retrieved."),
 "TRAF4":   ("none","none","no","no","na","","TNF-receptor-associated adaptor; not a transcription factor."),
 "ZBTB32":  ("weak","lymphocyte","no","no","context_dependent","PMC5869932","ZBTB family regulates lymphocyte responses (evidence mainly B-cell)."),
 "ZEB2":    ("strong","differentiation_effector","yes","yes","promotes_effector","PMC4647262","ZEB2 (downstream of T-bet) drives terminal CD8 effector/memory differentiation."),
 "ZNF121":  ("none","none","no","no","na","","No T-cell TF literature retrieved."),
 "ZNF706":  ("none","none","no","no","na","","No T-cell TF literature retrieved (Zbtb20 is the related memory-CD8 factor)."),
}

seeds, seed_dir, _ = bu.load_seeds()
files = {
 ("seed_excluded","GREmLN"):"results/tables/gremln_btla_vs_tcr_seed_excluded_tf_ranking.csv",
 ("seed_excluded","GENIE3"):"results/tables/genie3_btla_vs_tcr_seed_excluded_tf_ranking.csv",
 ("seed_inclusive","GREmLN"):"results/tables/gremln_btla_vs_tcr_seed_inclusive_tf_ranking.csv",
 ("seed_inclusive","GENIE3"):"results/tables/genie3_btla_vs_tcr_seed_inclusive_tf_ranking.csv",
}
rankcol = {"GREmLN":"gremln_csls_rank","GENIE3":"genie3_rank"}
top = {k: set(pd.read_csv(ROOT/v).sort_values(rankcol[k[1]]).head(25)["gene"]) for k,v in files.items()}

rows = []
for tf,(tier,phen,tcr,ae,direction,src,summ) in R.items():
    rows.append({
        "TF": tf, "is_BTLA_vs_TCR_seed": tf in seeds,
        "paperclip_evidence_tier": tier, "paperclip_primary_phenotype": phen,
        "tcr_checkpoint_evidence": tcr, "activation_exhaustion_evidence": ae,
        "btla_specific_evidence": "no",  # no paper directly implicates the TF in BTLA signalling
        "paperclip_direction": direction, "paperclip_key_source": src,
        "paperclip_short_rationale": summ,
        "in_gremln_top25_seed_excluded": tf in top[("seed_excluded","GREmLN")],
        "in_genie3_top25_seed_excluded": tf in top[("seed_excluded","GENIE3")],
        "in_gremln_top25_seed_inclusive": tf in top[("seed_inclusive","GREmLN")],
        "in_genie3_top25_seed_inclusive": tf in top[("seed_inclusive","GENIE3")],
        "paperclip_reviewed": True,
        "query_template": "'<TF> transcription factor CD4 T cell activation exhaustion checkpoint BTLA regulation' (+ TF-focused second pass for canonical TFs)",
    })
out = pd.DataFrame(rows).sort_values("TF")
dst = ROOT / "results" / "paperclip" / "paperclip_union_top25_review.csv"
out.to_csv(dst, index=False)

# coverage + tier comparison for the report (seed-excluded primary)
strong_mod = out[out["paperclip_evidence_tier"].isin(["strong","moderate"])]["TF"]
gr_only = top[("seed_excluded","GREmLN")] - top[("seed_excluded","GENIE3")]
g3_only = top[("seed_excluded","GENIE3")] - top[("seed_excluded","GREmLN")]
shared = top[("seed_excluded","GREmLN")] & top[("seed_excluded","GENIE3")]
def lit(s): return sorted(set(s) & set(strong_mod))
print("coverage:", len(out), "TFs, all reviewed =", bool(out["paperclip_reviewed"].all()))
print("tier counts:", out["paperclip_evidence_tier"].value_counts().to_dict())
print(f"seed-excluded shared with strong/moderate lit: {len(lit(shared))}/{len(shared)} -> {lit(shared)}")
print(f"GREmLN-only with strong/moderate lit: {len(lit(gr_only))}/{len(gr_only)} -> {lit(gr_only)}")
print(f"GENIE3-only with strong/moderate lit: {len(lit(g3_only))}/{len(g3_only)} -> {lit(g3_only)}")
print("wrote", dst)
