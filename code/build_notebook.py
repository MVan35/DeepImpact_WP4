# -*- coding: utf-8 -*-
"""Builds Paper_figures.ipynb (figures and Python-side tables, organised by the manuscript). python build_notebook.py"""
import json, uuid
cells = []
def md(t): cells.append({"cell_type": "markdown", "metadata": {}, "id": uuid.uuid4().hex[:8], "source": t.strip("\n")})
def code(t): cells.append({"cell_type": "code", "metadata": {}, "id": uuid.uuid4().hex[:8], "execution_count": None, "outputs": [], "source": t.strip("\n")})

md("""
# Neighbouring weed species influence the root microbiota and growth of wheat — figures and Python-side tables

Sections follow the manuscript. Run `Paper_analyses.qmd` first: this notebook reads its outputs from `results/`.
Section 3.1.2 also calls R (`pooled_check.R`, `vs_alone.R`; DESeq2, edgeR, metagenomeSeq). The last section copies every display item to
`results/display_items/` under its manuscript name.
""")
md("## Setup")
code('''
import os, sys, shutil, subprocess, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT    = os.environ.get("PAPER_PROJECT", "/Users/mathieuvanhove/Desktop/Work/DeepImpact/Weeds_wheat")
SUBMISSION = os.path.join(PROJECT, "submission")
DATA_DIR   = os.environ.get("PAPER_DATA", os.path.join(SUBMISSION, "data"))   # or the unpacked Zenodo data deposit
RES_DIR    = os.path.join(DATA_DIR, "results_v4")
MATRIX_DIR = os.path.join(DATA_DIR, "matrice_phase1", "df")
RESULTS    = os.environ.get("PAPER_RESULTS", os.path.join(SUBMISSION, "results"))   # written by Paper_analyses.qmd
OUT_DIR    = os.environ.get("PAPER_OUT", RESULTS)                                    # written here
FIG, TAB   = os.path.join(OUT_DIR, "figures"), os.path.join(OUT_DIR, "tables")
SUPP       = os.path.join(RESULTS, "supporting")
REV        = os.path.join(SUPP, "revision")            # revision_analyses.R (qmd 3.1.1)
TRAITS     = os.path.join(SUPP, "trait_models")        # trait_models_v3.R (qmd 3.3)
for d in (FIG, TAB): os.makedirs(d, exist_ok=True)
PAPER1     = os.path.join(PROJECT, "paper_1", "data")  # fungal DA exports and NST results, until staged
# edgeR results depend on the edgeR version (see README). The exports used for the paper are deposited in
# data/da_exports_paper and are used first when present; set PAPER_EDGER=regenerated to use the qmd re-export instead.
PAPER_EDGER = os.path.join(DATA_DIR, "da_exports_paper")
DA_ROOTS   = ([PAPER_EDGER] if os.path.isdir(PAPER_EDGER) and os.environ.get("PAPER_EDGER") != "regenerated" else []) + \\
             [os.path.join(RESULTS, "da_exports"), RES_DIR, os.path.join(PAPER1, "results_v4")]
NST_DIR    = os.path.join(DATA_DIR, "NST") if os.path.isdir(os.path.join(DATA_DIR, "NST")) else os.path.join(PAPER1, "NST")
RSCRIPT    = os.environ.get("RSCRIPT", shutil.which("Rscript") or "Rscript")
sys.path.insert(0, os.getcwd())
NEEDED = ["Table1_variance_partition.csv", "Table_S5_biomass_associated_taxa.csv", "Table_S6_biomass_ancova.csv", "da_exports/edgeR_species",
          "supporting/revision/fig4_direct_test.csv", "supporting/trait_models/trait_joint.csv"]
missing = [f for f in NEEDED if not os.path.exists(os.path.join(RESULTS, f))]
assert not missing, f"run Paper_analyses.qmd first; missing: {missing}"
import main_figures as mf, si_figures as sf, da_consensus as dc, biomass_tests as bt, control_size as csz
''')

md("## 3.1.1 Community composition and assembly")
md("### Fig. 1; Fig. S5 (α-diversity); Fig. S6 (Glomeromycota); Fig. S7 (bacterial assembly metrics); turnover share")
code('''
fig1 = mf.fig1(MATRIX_DIR, os.path.join(FIG, "Fig1"), mf.table1_from_csv(os.path.join(RESULTS, "Table1_variance_partition.csv")))
import fig_alpha_diversity as fa
alpha = fa.load_alpha(data_dir=MATRIX_DIR); fa.alpha_table(alpha, out_dir=TAB)
fa.fig_alpha_diversity(alpha, out_dir=FIG, fname="FigS5_alpha_diversity.png"); plt.close("all")
figS6, glom = sf.fig_glomeromycota(MATRIX_DIR, os.path.join(FIG, "FigS6_glomeromycota"))
figS7 = sf.fig_assembly_bacteria(NST_DIR, MATRIX_DIR, os.path.join(FIG, "FigS7_assembly_bacteria"))
turn = sf.turnover_table(MATRIX_DIR); turn.to_csv(os.path.join(TAB, "turnover_share.csv"), index=False)
print(glom.round(2)); turn
''')

md("## 3.1.2 Neighbour-enriched taxa")
md("### Fig. 2, Table S3; Methods S3 method agreement and cross-family convergence")
code('''
DA_OUT = os.path.join(OUT_DIR, "differential_abundance"); os.makedirs(DA_OUT, exist_ok=True)
r = subprocess.run([sys.executable, "da_core_check.py", MATRIX_DIR, DA_OUT] + DA_ROOTS, capture_output=True, text=True)
assert r.returncode == 0, r.stderr[-2000:]
conv = pd.read_csv(os.path.join(DA_OUT, "da_consensus_definitions.csv"))
conv[conv.definition.str.startswith("D6")][["kingdom", "soil", "calls", "ASVs", "ASVs_in_3plus_families", "expected_3plus", "P_convergence"]]
''')
md("### Methods S3, Table S3: control-size check (all weeds pooled against wheat × wheat; each call against the other weeds; each call against wheat alone)")
code('''
CS_OUT = os.path.join(SUPP, "control_size"); os.makedirs(CS_OUT, exist_ok=True)
fig2v = pd.read_csv(os.path.join(DA_OUT, "Fig2_values.csv"))
_spec = os.path.join(CS_OUT, "fig2_neighbour_vs_other_weeds.csv"); _k = ["kingdom", "soil", "neighbour", "ASV"]
_fresh = os.path.exists(os.path.join(CS_OUT, "vs_wheat_alone_all.csv")) and os.path.exists(_spec) and set(map(tuple, pd.read_csv(_spec)[_k].values)) == set(map(tuple, fig2v[_k].values))
if not _fresh:                                                    # c. 5-20 min; rerun whenever the Fig. 2 calls change
    csz.run(MATRIX_DIR, os.path.join(DA_OUT, "Fig2_values.csv"), CS_OUT, RSCRIPT)
pooled = pd.read_csv(os.path.join(CS_OUT, "pooled_weeds_vs_wxw.csv")); spec = pd.read_csv(os.path.join(CS_OUT, "fig2_neighbour_vs_other_weeds.csv"))
alone = pd.read_csv(os.path.join(CS_OUT, "vs_wheat_alone_all.csv"))
calls = csz.classify(fig2v, pooled, spec, alone)
csz.table_s3(os.path.join(DA_OUT, "Table_S7_enriched_ASVs.xlsx"), calls, pooled, os.path.join(TAB, "Table_S3_enriched_ASVs.xlsx"))
csz.summary(calls)
''')
md("### Results 3.2.2: overlap of the Fig. 2 ASVs with the TITAN2 indicators")
code('''
r = subprocess.run([sys.executable, "core_vs_indicators.py", MATRIX_DIR, os.path.join(RESULTS, "Table_S5_biomass_associated_taxa.csv"),
                    os.path.join(DA_OUT, "Fig2_values.csv"), os.path.join(TAB, "fig2_vs_indicators.csv")], capture_output=True, text=True)
assert r.returncode == 0, r.stderr[-1500:]
pd.read_csv(os.path.join(TAB, "fig2_vs_indicators.csv"))[["kingdom", "soil", "z+_shared", "z+_expected", "z+_P", "z+_ASVs"]]
''')

md("## 3.2.1 Wheat biomass")
md("### Fig. 3 (between-soil tests, BH across 30 treatments); log-scale ANOVA; legume test; Table S4")
code('''
META_B = os.path.join(MATRIX_DIR, "metadata_bacteria_1314ASVs.csv")
m = bt.load(META_B)
fig3, soil_dependent = mf.fig3(MATRIX_DIR, os.path.join(FIG, "Fig3"))
aov = bt.anova_tables(m); aov.to_csv(os.path.join(TAB, "biomass_anova.csv"), index=False)
bt.scale_check(m).to_csv(os.path.join(TAB, "biomass_scale_check.csv"), index=False)   # Methods 2.7 and S3: why log scale
between = bt.between_soils(m); between.to_csv(os.path.join(TAB, "biomass_between_soils.csv"), index=False)
leg = bt.legume_test(m); leg.to_csv(os.path.join(TAB, "biomass_legumes.csv"), index=False)
t4 = bt.table_s3(m, os.path.join(TAB, "Table_S4_biomass_contrasts.xlsx"))     # Table S4 in the SI
print("asterisks:", soil_dependent, "| above W×W (q < 0.05):", int((t4.q_vs_wxw < 0.05).sum()), "| above wheat alone:", int((t4.delta_vs_alone_g > 0).sum()))
aov[aov.scale == "log"].round(3)
''')

md("## 3.2.2 Microbiota shifts and wheat growth")
md("### Fig. 4 (species-centroid db-RDA on species-mean biomass; tests from qmd 3.1.1); Table S6 (sensitivity of the Fig. 4 tests); Table S7 (TITAN2 indicators, exploratory)")
code('''
fig4 = mf.fig4(os.path.join(REV, "fig4_scores.csv"), os.path.join(REV, "fig4_direct_test.csv"), os.path.join(FIG, "Fig4"))
s5 = pd.read_csv(os.path.join(RESULTS, "Table_S5_biomass_associated_taxa.csv"))
t6 = s5[s5.TITAN2_direction.isin(["z+", "z-"])][["kingdom", "soil", "ASV", "Genus", "TITAN2_direction", "spearman_rho_biomass"]].copy()
t6["ASV"] = t6.ASV.str.replace("Cluster_", "ASV_")
t6["in_Fig2"] = [(k, a) in set(zip(fig2v.kingdom, fig2v.ASV.str.replace("Cluster_", "ASV_"))) for k, a in zip(t6.kingdom, t6.ASV)]
t6["responds_in_both_soils"] = t6.groupby(["kingdom", "ASV"]).soil.transform("nunique") > 1
t6.to_excel(os.path.join(TAB, "Table_S7_TITAN2_indicators.xlsx"), index=False)
t6s = pd.read_csv(os.path.join(REV, "fig4_direct_test.csv"))[["kingdom", "soil", "R2_pct", "P", "P_given_monocot", "P_given_monocot_legume", "P_given_Btot", "P_without_legumes"]]
t6s.round(3).to_csv(os.path.join(TAB, "Table_S6_fig4_sensitivity.csv"), index=False)
print(t6.groupby(["kingdom", "soil", "TITAN2_direction"]).size())
pd.read_csv(os.path.join(REV, "fig4_direct_test.csv")).round(4)
''')

md("## 3.3 Weed traits and phylogeny")
md("### Table 2")
code('''
comp = pd.read_csv(os.path.join(TRAITS, "trait_comp.csv")); joint = pd.read_csv(os.path.join(TRAITS, "trait_joint.csv"))
mant = pd.read_csv(os.path.join(TRAITS, "trait_mantel.csv"))
star = lambda p: "***" if p <= 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
cells_ = [(k, s) for k in ("bacteria", "fungi") for s in ("Cereal-dominated (AF193)", "Mixed-system (AF190)")]
lab = [("Phylogeny", "Distance to wheat", "Distance"), ("", "Monocot vs eudicot", "Monocot"), ("", "Segetal status", "Segetal"),
       ("Habitat", "Ellenberg N", "Nitrogen"), ("", "Ellenberg L", "Light"), ("Root traits", "SRL", "SRL"), ("", "RDMC", "RDMC"), ("", "Root:shoot", "RootShoot")]
rows = []
for c1, c2, k in lab:
    rows.append([c1, c2] + [f"{r.F:.2f}{star(r.P_species)}" for r in [comp[(comp.kingdom == a) & (comp.soil == b) & (comp.trait == k)].iloc[0] for a, b in cells_]])
J = joint.set_index(["kingdom", "soil"]); MT = mant.set_index(["kingdom", "soil"])
rows.append(["Joint model", "R² of seven traits (%)"] + [f"{100 * J.loc[c].R2_total:.1f} (P = {J.loc[c].P_species:.2f})" for c in cells_])
rows.append(["", "R² of species identity (%)"] + [f"{100 * J.loc[c].R2_species:.1f}" for c in cells_])
rows.append(["Phylogenetic signal", "Mantel r"] + [f"{MT.loc[c].mantel_r:.2f}{star(MT.loc[c].P)}" for c in cells_])
t2 = pd.DataFrame(rows, columns=["Trait category", "Variable", "Bacteria AF193", "Bacteria AF190", "Fungi AF193", "Fungi AF190"])
t2.to_excel(os.path.join(TAB, "Table_2_traits.xlsx"), index=False); t2
''')

md("## Graphical abstract")
code('''
import build_ga
build_ga.build(path=os.path.join(FIG, "graphical_abstract.png")); plt.close("all")
''')

md("""
## Display items

Every item under its manuscript name in `display_items/`. Figs S1–S3 are drawings or images supplied with the manuscript
(design, soil characterisation, rarefaction); Fig. S4 and Tables 1 and S5 are written by `Paper_analyses.qmd`.
""")
code('''
ITEMS = {
 "Fig_1": os.path.join(FIG, "Fig1.png"), "Fig_2": os.path.join(DA_OUT, "Fig2.png"), "Fig_3": os.path.join(FIG, "Fig3.png"), "Fig_4": os.path.join(FIG, "Fig4.png"),
 "Table_1": os.path.join(RESULTS, "Table1_variance_partition.csv"), "Table_2": os.path.join(TAB, "Table_2_traits.xlsx"),
 "Fig_S4": os.path.join(RESULTS, "FigS4_permdisp.png"), "Fig_S5": os.path.join(FIG, "FigS5_alpha_diversity.png"),
 "Fig_S6": os.path.join(FIG, "FigS6_glomeromycota.png"), "Fig_S7": os.path.join(FIG, "FigS7_assembly_bacteria.png"),
 "Table_S3": os.path.join(TAB, "Table_S3_enriched_ASVs.xlsx"), "Table_S4": os.path.join(TAB, "Table_S4_biomass_contrasts.xlsx"),
 "Table_S5": os.path.join(RESULTS, "Table_S6_biomass_ancova.csv"), "Table_S6": os.path.join(TAB, "Table_S6_fig4_sensitivity.csv"), "Table_S7": os.path.join(TAB, "Table_S7_TITAN2_indicators.xlsx"),
 "Graphical_abstract": os.path.join(FIG, "graphical_abstract.png")}
DI = os.path.join(OUT_DIR, "display_items"); os.makedirs(DI, exist_ok=True)
status = {}
for name, src in ITEMS.items():
    stem, ext = os.path.splitext(src)
    for e in ([ext] + [x for x in (".pdf", ".svg") if ext == ".png"]):
        if os.path.exists(stem + e): shutil.copy(stem + e, os.path.join(DI, name + e))
    status[name] = os.path.exists(src)
missing = [k for k, v in status.items() if not v]
print(f"{sum(status.values())} of {len(status)} display items written to {DI}; missing: {missing}")
''')
nb = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                   "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
json.dump(nb, open("Paper_figures.ipynb", "w"), indent=1, ensure_ascii=False)
print("cells", len(cells))
