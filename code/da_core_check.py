# -*- coding: utf-8 -*-
"""Differential abundance for Results 3.1.2: Fig. 2, the enriched-ASV workbook (Table S3) and the Methods S3 tests.
A call needs enrichment under the neighbour (log2 fold change > 2) by at least two of DESeq2, edgeR and metagenomeSeq:
at adjusted P < 0.05 for Fig. 2 (DESeq2, edgeR BH; metagenomeSeq Bonferroni), at nominal method-specific thresholds for the broader screen (Table S3).
Methods S3: method agreement against chance (hypergeometric) and cross-family convergence (999 permutations).
Also draws a genus-level map of the screen (FigS18_genus_screen; not a display item).
Usage (from this folder): python da_core_check.py <matrix_dir> <out_dir> <export_root> [<export_root> ...]"""
import os, sys, io, contextlib
sys.path.insert(0, os.getcwd())
import numpy as np, pandas as pd
from scipy.stats import hypergeom
import da_consensus as dc

MATRIX, OUT, ROOTS = sys.argv[1], sys.argv[2], sys.argv[3:]
os.makedirs(OUT, exist_ok=True)
TT = dc.load_tables(MATRIX)
rng = np.random.default_rng(1)


def calls_full(kg, soil):
    out = {}
    for m, (sub, pat) in dc.FILES.items():
        f = next(os.path.join(r, sub, pat.format(kg=kg, soil=soil)) for r in ROOTS
                 if os.path.exists(os.path.join(r, sub, pat.format(kg=kg, soil=soil))))
        raw = pd.read_csv(f); df = dc.read_export(f, m)
        df["q"] = raw.padj.values if m == "DESeq2" else (raw.FDR.values if m == "edgeR" else df.p)   # metagenomeSeq p is adjusted
        df = df[df.weed != "wheat alone"].copy()
        sign, _, _ = dc.orient_and_check(TT, kg, soil, df, os.path.basename(f)); df["l"] = sign * df.lfc
        out[m] = df[["weed", "Cluster", "p", "q", "l"]].copy()
    return out


CF = {(kg, s): calls_full(kg, s) for kg in ("bacteria", "fungi") for s in ("Resistant", "Permissive")}


def universe(kg, soil, w):
    X = TT[kg]["X"]; md = TT[kg]["md"]
    idx = md[(md.soil == dc.SOILCODE[soil]) & md.species.isin([w, dc.CONTROL])].index.intersection(X.index)
    return set(X.columns[(X.loc[idx] > 0).any()])


UNI = {(kg, s, w): universe(kg, s, w) for (kg, s), d in CF.items() for w in d["metagenomeSeq"].weed.unique()}


def called(d, m, mode, w, mthr=None):
    df = d[m]; df = df[df.weed == w]; _, pthr, lthr = dc.THR[m]
    if m == "metagenomeSeq" and mthr is not None: lthr = mthr
    sel = (df.p < pthr) if mode == "raw" else (df.q < 0.05)
    return set(df.loc[sel & (df.l > lthr), "Cluster"])


rows = []
for mode, mthr in (("raw", None), ("BH", None), ("raw", 2.0), ("BH", 2.0)):
    for (kg, s), d in CF.items():
        for A, B in (("DESeq2", "edgeR"), ("DESeq2", "metagenomeSeq"), ("edgeR", "metagenomeSeq")):
            o = e = nt = ns = 0
            for w in d["metagenomeSeq"].weed.unique():
                U = UNI[(kg, s, w)]; a = called(d, A, mode, w, mthr) & U; b = called(d, B, mode, w, mthr) & U
                if not a or not b: continue
                N, k = len(U), len(a & b); o += k; e += len(a) * len(b) / N; nt += 1; ns += hypergeom(N, len(a), len(b)).sf(k - 1) < 0.05
            rows.append(dict(kingdom=kg, soil=dc.SOILLABEL[s], threshold=mode + ("" if mthr is None else f", metagenomeSeq log2FC>{mthr:g}"),
                             pair=f"{A} ∩ {B}", neighbours_tested=nt, observed=o, expected=round(e, 1), fold=round(o / e, 1) if e else np.nan,
                             neighbours_p05=ns))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "da_method_agreement.csv"), index=False)

DEFS = {"D1 raw, ≥2 of 3 (v1.7-4)": ("raw", None, None, 2), "D2 raw, ≥2 of 3, common log2FC > 2 (adopted screen)": ("raw", 2.0, None, 2),
        "D3 raw, DESeq2 ∩ edgeR": ("raw", None, ("DESeq2", "edgeR"), 2), "D4 raw, all three, common log2FC > 2": ("raw", 2.0, None, 3),
        "D5 FDR, ≥2 of 3": ("BH", None, None, 2), "D6 FDR, ≥2 of 3, common log2FC > 2 (adopted core, Fig. 2)": ("BH", 2.0, None, 2)}
fam = dict(zip(TT["bacteria"]["md"].species, TT["bacteria"]["md"].family))


def pairs_for(dfn):
    mode, mthr, meths, kmin = dfn; out = {}
    for (kg, s), d in CF.items():
        P = set()
        for w in d["metagenomeSeq"].weed.unique():
            U = UNI[(kg, s, w)]; sets = [called(d, m, mode, w, mthr) & U for m in (meths or dc.FILES)]
            P |= {(w, c) for c in set().union(*sets) if sum(c in x for x in sets) >= kmin}
        out[(kg, s)] = P
    return out


def convergence(P, kg, s, nperm=999, kfam=3):
    if not P: return 0, np.nan, np.nan
    byw = pd.Series([w for w, _ in P]).value_counts()
    def count(pp):
        df = pd.DataFrame(pp, columns=["w", "c"]); df["f"] = df.w.map(fam); return int((df.groupby("c").f.nunique() >= kfam).sum())
    obs = count(list(P)); null = []
    for _ in range(nperm):
        pp = []
        for w, k in byw.items():
            U = sorted(UNI[(kg, s, w)]); pp += [(w, c) for c in rng.choice(U, size=k, replace=False)]
        null.append(count(pp))
    null = np.array(null); return obs, float(null.mean()), float((np.sum(null >= obs) + 1) / (nperm + 1))


def prev(kg, s, asvs, sp):
    X = TT[kg]["X"]; md = TT[kg]["md"]; idx = md[(md.soil == dc.SOILCODE[s]) & (md.species == sp)].index.intersection(X.index)
    return (X.loc[idx, list(asvs)] > 0).mean()


rows, PAIRS = [], {}
for name, dfn in DEFS.items():
    PP = pairs_for(dfn); PAIRS[name] = PP
    for (kg, s), P in PP.items():
        asv = {c for _, c in P}; nn = pd.Series([c for _, c in P]).value_counts() if P else pd.Series(dtype=int)
        o, e, p = convergence(P, kg, s)
        pw = prev(kg, s, asv, dc.CONTROL) if asv else pd.Series(dtype=float); pa = prev(kg, s, asv, "wheat alone") if asv else pd.Series(dtype=float)
        rows.append(dict(definition=name, kingdom=kg, soil=dc.SOILLABEL[s], calls=len(P), ASVs=len(asv),
                         genera=len({TT[kg]["tax"].Genus.get(c) for c in asv}),
                         median_neighbours_per_ASV=float(nn.median()) if len(nn) else np.nan,
                         ASVs_in_3plus_families=o, expected_3plus=round(e, 1) if e == e else np.nan, P_convergence=round(p, 3) if p == p else np.nan,
                         pct_absent_from_all_wxw_pots=round(100 * (pw == 0).mean()) if len(pw) else np.nan,
                         pct_absent_wxw_but_present_wheat_alone=round(100 * ((pw == 0) & (pa > 0)).mean()) if len(pw) else np.nan))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "da_consensus_definitions.csv"), index=False)

core = PAIRS["D6 FDR, ≥2 of 3, common log2FC > 2 (adopted core, Fig. 2)"]
lfc = {(kg, s, w, c): l for (kg, s), d in CF.items() for w, c, l in zip(d["DESeq2"].weed, d["DESeq2"].Cluster, d["DESeq2"].l)}
fam28 = {w: fam[w] for w in sorted(set(TT["bacteria"]["md"].species) - {"wheat x wheat", "wheat alone"})}
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import si_figures as sf, bacteria_dotplot_revised as bdp, fungi_dotplot_cell as fdp
import main_figures as mf
mf.fig2(core, lfc, TT, fam28, os.path.join(OUT, "Fig2"))   # Times New Roman
scr = PAIRS["D2 raw, ≥2 of 3, common log2FC > 2 (adopted screen)"]
L = pd.DataFrame([dict(Cluster=c, Plant_Species=w, Soil=s, Kingdom=kg, Method="consensus") for (kg, s), P in scr.items() for w, c in P])
L.to_csv(os.path.join(OUT, "da_consensus_screen_long.csv"), index=False)
fi = dc.fig2_inputs(L, TT, min_combos_bacteria=5)
with contextlib.redirect_stdout(io.StringIO()):
    fa, _ = bdp.create_combined_bacteria_plot(fi["bact"], fi["genus_b"], fi["phylum_b"], fi["family"], dc.FUNCTION_BY_GENUS,
                                              dc.FUNCTION_PALETTE, dc.BACT_PHYLUM_COLOURS)
    fb, _ = fdp.create_combined_fungi_plot(fi["fungi"], fi["genus_f"], fi["phylum_f"], fi["family"])
    dc.match_fig2_panels(fa, fb); dc.save_fig2(fa, fb, os.path.join(OUT, "FigS17_genus_screen"))
plt.close("all")
# ---------------------------------------------------------------- Table S7 workbook (v2.3)
def tier_calls(mode):
    rows = []
    for (kg, s), d in CF.items():
        for w in d["metagenomeSeq"].weed.unique():
            U = UNI[(kg, s, w)]; sets = {m: called(d, m, mode, w, 2.0) & U for m in dc.FILES}
            for c in set().union(*sets.values()):
                ms = [m for m in dc.FILES if c in sets[m]]
                if len(ms) >= 2:
                    rows.append(dict(kingdom=kg, soil=dc.SOILLABEL[s], neighbour=w, ASV=c, methods=", ".join(ms)))
    return pd.DataFrame(rows)


SCR, FDR = tier_calls("raw"), tier_calls("BH")
assert {(k, w, c) for k, w, c in zip(FDR.kingdom, FDR.neighbour, FDR.ASV)} == {(k, w, c) for (k, s), P in core.items() for w, c in P}
fdr_keys = set(zip(FDR.kingdom, FDR.soil, FDR.neighbour, FDR.ASV))
SCR["FDR_supported"] = ["yes" if k in fdr_keys else "no" for k in zip(SCR.kingdom, SCR.soil, SCR.neighbour, SCR.ASV)]


def per_asv(kg):
    g = SCR[SCR.kingdom == kg]; tx = TT[kg]["tax"]; out = []
    for c, h in g.groupby("ASV"):
        hf = h[h.FDR_supported == "yes"]
        out.append(dict(ASV=c.replace("Cluster_", "ASV_"), Phylum=tx.Phylum.get(c), Genus=tx.Genus.get(c), Species=tx.Species.get(c),
                        Soil=", ".join(sorted(h.soil.unique())), n_neighbours=h.neighbour.nunique(),
                        n_families=h.neighbour.map(fam).nunique(),
                        Method_combinations="; ".join(sorted(h.methods.unique())),
                        Neighbours="; ".join(f"{s}: " + ", ".join(sorted(x.neighbour)) for s, x in h.groupby("soil")),
                        In_Fig2="yes" if len(hf) else "no",
                        Neighbours_FDR="; ".join(f"{s}: " + ", ".join(sorted(x.neighbour)) for s, x in hf.groupby("soil"))))
    return pd.DataFrame(out).sort_values(["Phylum", "Genus", "ASV"])


ov = pd.concat([SCR.assign(tier="unadjusted P (screen)"), FDR.assign(tier="FDR q < 0.05 (Fig. 2)")])
ov = ov.groupby(["tier", "kingdom", "soil", "methods"]).size().unstack("methods", fill_value=0).reset_index()
ov["calls"] = ov[[c for c in ov.columns if "," in c]].sum(1)
ag = pd.read_csv(os.path.join(OUT, "da_method_agreement.csv"))
ag = ag[ag.threshold.isin(["raw, metagenomeSeq log2FC>2", "BH, metagenomeSeq log2FC>2"])].copy()
ag["threshold"] = ag.threshold.map({"raw, metagenomeSeq log2FC>2": "unadjusted P (screen)", "BH, metagenomeSeq log2FC>2": "FDR q < 0.05 (Fig. 2)"})
readme = pd.DataFrame({"Table S7": [
    "Neighbour-enriched ASVs in the wheat root endosphere. An ASV × neighbour × soil call requires enrichment under the neighbour relative to "
    "wheat × wheat by at least two of DESeq2, edgeR and metagenomeSeq, each with log2 fold change > 2.",
    "Fig. 2 (adjusted P): DESeq2 padj and edgeR FDR (Benjamini-Hochberg), metagenomeSeq Bonferroni-adjusted P, each < 0.05. Screen (unadjusted): DESeq2 P < 0.05, edgeR P < 0.01, "
    "metagenomeSeq adjusted P < 0.05. Every Fig. 2 call is also a screen call.",
    "The fold-change sign of every export was checked against log2(neighbour/control) computed from the abundance table.",
    "Sheets: bacteria, fungi (one row per ASV; In_Fig2 = passes at FDR); calls (one row per ASV × neighbour × soil); method_overlap "
    "(calls by method combination, both tiers); method_agreement (observed vs chance overlap of each pair of methods; Methods S3)."]})
xl = os.path.join(OUT, "Table_S7_enriched_ASVs.xlsx")
with pd.ExcelWriter(xl) as w:
    readme.to_excel(w, sheet_name="README", index=False)
    per_asv("bacteria").to_excel(w, sheet_name="bacteria", index=False); per_asv("fungi").to_excel(w, sheet_name="fungi", index=False)
    SCR.assign(ASV=SCR.ASV.str.replace("Cluster_", "ASV_")).sort_values(["kingdom", "soil", "neighbour", "ASV"]).to_excel(w, sheet_name="calls", index=False)
    ov.to_excel(w, sheet_name="method_overlap", index=False)
    ag[["threshold", "kingdom", "soil", "pair", "neighbours_tested", "observed", "expected", "fold", "neighbours_p05"]].to_excel(w, sheet_name="method_agreement", index=False)
print("Table S7:", xl)

print("done:", OUT)
