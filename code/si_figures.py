# -*- coding: utf-8 -*-
"""Supporting Information figures (Figs S8, S9, S11, S12, S14, S16) and the ASV panel used by Fig. 2.
Called from Paper_figures.ipynb. Soil colours: cereal-dominated (AF193) blue, mixed-system (AF190) orange."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

SOILCODE = {"resistant": "Cereal-dominated (AF193)", "Permissive": "Mixed-system (AF190)"}
SOILCOL = {"Cereal-dominated (AF193)": "#3B75AF", "Mixed-system (AF190)": "#E08A2C"}
LIGHT = {"Cereal-dominated (AF193)": "#A9C4E2", "Mixed-system (AF190)": "#F3C994"}
CTRL = ["wheat x wheat", "wheat alone"]
RC = {"font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
      "legend.fontsize": 6, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
      "font.family": "sans-serif", "pdf.fonttype": 42, "svg.fonttype": "none"}
FAMCOL = {"Apiaceae": "#d62728", "Asteraceae": "#1f77b4", "Boraginaceae": "#bcbd22", "Brassicaceae": "#2ca02c",
          "Campanulaceae": "#e377c2", "Caryophyllaceae": "#98df8a", "Fabaceae": "#7f7f7f", "Geraniaceae": "#9467bd",
          "Juncaceae": "#8c564b", "Papaveraceae": "#17becf", "Plantaginaceae": "#ff9896", "Poaceae": "#ff7f0e",
          "Polygonaceae": "#ffbb78", "Ranunculaceae": "#aec7e8", "Wheat": "#000000"}
# old labels in the NST exports (same mislabelling as the old DA exports, see da_consensus.SYN)
NST_SYN = {"Centaurea cyanus": "Cyanus segetum", "Bromus sterilis": "Anisantha sterilis", "Festuca arundinacea": "Lolium arundinaceum",
           "Festuca myuros": "Lolium arundinaceum", "Viola arvensis": "Vulpia myuros"}


def _save(fig, out):
    for ext in ("png", "pdf", "svg"): fig.savefig(f"{out}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)
    return out + ".png"


def _meta(DF, kg="bacteria"):
    n = 1314 if kg == "bacteria" else 647
    return pd.read_csv(os.path.join(DF, f"metadata_{kg}_{n}ASVs.csv"))


def species_order(meta):
    f = meta[~meta.species.isin(CTRL)].groupby("species").family.first()
    return list(f.reset_index().sort_values(["family", "species"]).species), f.to_dict()


def beta_pair(PA):
    A = PA.astype(np.int32); a = A @ A.T; s = A.sum(1); b = s[:, None] - a; c = s[None, :] - a
    mn = np.minimum(b, c)
    with np.errstate(invalid="ignore", divide="ignore"):
        sor = (b + c) / (2 * a + b + c); sim = mn / (a + mn)
    iu = np.triu_indices(len(A), 1)
    return np.nanmean(sim[iu]), np.nanmean(sor[iu] - sim[iu]), np.nanmean(sor[iu])


def turnover_table(DF):
    rows = []
    for kg, n in (("bacteria", 1314), ("fungi", 647)):
        X = pd.read_csv(os.path.join(DF, f"otu_{kg}_table_R_{n}ASVs.csv"), sep="\t", index_col=0).T.astype(float)
        md = _meta(DF, kg).set_index("sample"); md = md.loc[md.index.intersection(X.index)]
        for code, sl in SOILCODE.items():
            m = md[md.soil == code]; PA = (X.loc[m.index] > 0).values
            t, ne, tot = beta_pair(PA); rows.append(dict(kingdom=kg, soil=sl, group="All pots", n=len(m), turnover=t, nestedness=ne, total=tot))
            grp = m.species.map({"wheat x wheat": "Wheat × wheat", "wheat alone": "Wheat alone"}).fillna(m.family)
            for g in grp.unique():
                idx = (grp == g).values
                if idx.sum() < 3: continue
                t, ne, tot = beta_pair(PA[idx]); rows.append(dict(kingdom=kg, soil=sl, group=g, n=int(idx.sum()), turnover=t, nestedness=ne, total=tot))
    d = pd.DataFrame(rows); d["pct_turnover"] = 100 * d.turnover / d.total
    return d


def fig_turnover(DF, out):
    d = turnover_table(DF); d.to_csv(out + "_values.csv", index=False)
    fams = sorted(set(d.group) - {"All pots", "Wheat × wheat", "Wheat alone"})
    order = ["All pots"] + fams + ["Wheat alone", "Wheat × wheat"]
    with plt.rc_context(RC):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.4), sharey=True, gridspec_kw=dict(wspace=0.08))
        for ax, kg in zip(axes, ("bacteria", "fungi")):
            g = d[d.kingdom == kg].set_index(["group", "soil"])
            for i, grp in enumerate(order):
                for k, s in enumerate(SOILCOL):
                    if (grp, s) not in g.index: continue
                    r = g.loc[(grp, s)]; y = i + (k - 0.5) * 0.38
                    ax.barh(y, r.turnover, height=0.36, color=SOILCOL[s], lw=0)
                    ax.barh(y, r.nestedness, left=r.turnover, height=0.36, color=LIGHT[s], lw=0)
                    ax.text(r.total + 0.01, y, f"{r.pct_turnover:.0f}%", va="center", fontsize=5.5, color="#333333")
            ax.axhline(0.5, color="#999999", lw=0.5); ax.axhline(len(order) - 2.5, color="#999999", lw=0.5)
            ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
            ax.set_xlim(0, 0.95); ax.set_xlabel("Mean pairwise Sørensen dissimilarity")
            ax.set_title(("Bacteria (rpoB)" if kg == "bacteria" else "Fungi (ITS1)"), loc="left", fontweight="bold")
        axes[0].invert_yaxis()
        h = [Patch(color=SOILCOL[s], label=f"{s}: turnover") for s in SOILCOL] + \
            [Patch(facecolor=LIGHT[s], label=f"{s}: nestedness") for s in SOILCOL]
        fig.legend(handles=h, frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.09))
        fig.text(0.99, 0.005, "labels: turnover as % of total dissimilarity", ha="right", fontsize=5.5, color="#555555")
        return _save(fig, out), d


def nst_table(NST_DIR):
    t = []
    for kg, f1, f2 in (("bacteria", "NST_bacteria.csv", "pNST_bacteria.csv"), ("fungi", "NST_fungi.csv", "pnst_fungi.csv")):
        a = pd.read_csv(os.path.join(NST_DIR, f1)); b = pd.read_csv(os.path.join(NST_DIR, f2))
        b["Species"] = b.Species.str.replace("_", " ").str.replace("bursa pastoris", "bursa-pastoris").str.replace("speculum veneris", "speculum-veneris").str.replace("foenum graecum", "foenum-graecum")
        m = a[["Soil", "Species", "NST"]].merge(b[["Soil", "Species", "pNST"]], on=["Soil", "Species"])
        m["Species"] = m.Species.replace(NST_SYN); m["kingdom"] = kg; m["soil"] = m.Soil.map(SOILCODE); t.append(m)
    return pd.concat(t, ignore_index=True)


def fig_nst(NST_DIR, DF, out):
    d = nst_table(NST_DIR); sp, fam = species_order(_meta(DF))
    dropped = sorted(set(d.Species) - set(sp) - set(CTRL)); d = d[d.Species.isin(sp + CTRL)]
    d.to_csv(out + "_values.csv", index=False)
    order = sp + CTRL; lab = [s if s not in CTRL else ("wheat × wheat" if s == "wheat x wheat" else "wheat alone") for s in order]
    with plt.rc_context(RC):
        fig, axes = plt.subplots(1, 4, figsize=(7.2, 4.8), sharey=True, gridspec_kw=dict(wspace=0.14))
        for ax, (kg, col, ttl) in zip(axes, (("bacteria", "NST", "Bacteria, taxonomic NST"), ("bacteria", "pNST", "Bacteria, phylogenetic NST"),
                                             ("fungi", "NST", "Fungi, taxonomic NST"), ("fungi", "pNST", "Fungi, phylogenetic NST"))):
            g = d[d.kingdom == kg]
            for k, s in enumerate(SOILCOL):
                gs = g[g.soil == s].set_index("Species")[col].reindex(order)
                ax.scatter(gs.values, np.arange(len(order)) + (k - 0.5) * 0.3, s=9, color=SOILCOL[s], lw=0, label=s, zorder=3)
            ax.axvline(0.5, color="black", lw=0.6, ls="--")
            ax.set_xlim(-0.02, 1.02); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"]); ax.set_title(ttl, loc="left"); ax.set_xlabel("NST")
            for j in range(1, len(order)):
                fa = fam.get(order[j], "Wheat"); fb = fam.get(order[j - 1], "Wheat")
                if fa != fb: ax.axhline(j - 0.5, color="#e5e5e5", lw=0.5, zorder=0)
            ax.text(0.25, -1.2, "deterministic", ha="center", fontsize=5.5, color="#555555"); ax.text(0.75, -1.2, "stochastic", ha="center", fontsize=5.5, color="#555555")
        axes[0].set_yticks(range(len(order))); axes[0].set_yticklabels(lab, style="italic"); axes[0].invert_yaxis()
        for t_, s in zip(axes[0].get_yticklabels(), order):
            t_.set_color(FAMCOL.get(fam.get(s, "Wheat"), "black")) if False else None
        fig.legend(handles=[Line2D([], [], marker="o", ls="", color=c, label=s, markersize=4) for s, c in SOILCOL.items()],
                   frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.05))
        return _save(fig, out), dropped


def fig_rci_family(DF, out):
    meta = _meta(DF).dropna(subset=["Dry_wheat_mass"]).copy(); meta["soilL"] = meta.soil.map(SOILCODE)
    alone = meta[meta.species == "wheat alone"].groupby("soilL").Dry_wheat_mass.mean()
    meta["RCI"] = (meta.soilL.map(alone) - meta.Dry_wheat_mass) / meta.soilL.map(alone)
    sp, fam = species_order(meta); order = sp + ["wheat x wheat"]
    meta[meta.species != "wheat alone"][["soilL", "species", "RCI"]].to_csv(out + "_values.csv", index=False)
    with plt.rc_context(RC):
        fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.4), sharex=True, sharey=True, gridspec_kw=dict(hspace=0.12))
        for ax, s in zip(axes, SOILCOL):
            mu = meta[(meta.soilL == s) & ~meta.species.isin(CTRL)].RCI.mean()
            for i, x in enumerate(order):
                v = meta.loc[(meta.soilL == s) & (meta.species == x), "RCI"].dropna().values
                c = FAMCOL[fam.get(x, "Wheat")]
                bp = ax.boxplot(v, positions=[i], widths=0.65, patch_artist=True, showfliers=True,
                                flierprops=dict(marker="o", markersize=1.8, markerfacecolor="black", markeredgecolor="none"),
                                medianprops=dict(color="black", lw=0.8), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
                bp["boxes"][0].set(facecolor=c, edgecolor="black", lw=0.5)
            ax.axhline(mu, color="black", lw=0.7, ls="--"); ax.axhline(0, color="#777777", lw=0.5, ls=":")
            ax.set_ylabel("Relative Competition Index (RCI)"); ax.set_title(s, loc="left", fontweight="bold", color=SOILCOL[s])
            ax.set_title(f"dashed line: mean RCI of neighbour pots ({mu:.2f}); dotted: RCI = 0 (wheat alone)", loc="right", fontsize=5.5)
        axes[1].set_xticks(range(len(order)))
        axes[1].set_xticklabels([("wheat × wheat" if o == "wheat x wheat" else o) for o in order], rotation=80, style="italic")
        hs = [Patch(facecolor=FAMCOL[f], edgecolor="black", lw=0.4, label=("Wheat (wheat × wheat)" if f == "Wheat" else f))
              for f in list(dict.fromkeys([fam[s] for s in sp])) + ["Wheat"]]
        fig.legend(handles=hs, title="Family", frameon=False, loc="center left", bbox_to_anchor=(0.91, 0.5))
        return _save(fig, out)


def fig_fdr_asv(P, lfc, T, fam, out):
    """P: dict (kingdom, soil) -> set of (neighbour, ASV); lfc: dict (kingdom, soil, neighbour, ASV) -> DESeq2 log2FC."""
    rows = [dict(kingdom=kg, soil=SOILCODE["resistant" if s == "Resistant" else "Permissive"], neighbour=w, ASV=c,
                 l=lfc.get((kg, s, w, c), np.nan)) for (kg, s), pp in P.items() for w, c in pp]
    d = pd.DataFrame(rows)
    d["genus"] = [T[k]["tax"].Genus.get(a) for k, a in zip(d.kingdom, d.ASV)]
    d["phylum"] = [T[k]["tax"].Phylum.get(a) for k, a in zip(d.kingdom, d.ASV)]
    d.to_csv(out + "_values.csv", index=False)
    nb = sorted(fam, key=lambda s: (fam[s], s))
    asv = d.groupby(["kingdom", "phylum", "genus", "ASV"]).neighbour.nunique().reset_index().sort_values(["kingdom", "phylum", "genus", "neighbour"], ascending=[True, True, True, False])
    lab = [f"{g} {a.replace('Cluster_', 'ASV_')}" for g, a in zip(asv.genus, asv.ASV)]
    ypos = {(k, a): i for i, (k, a) in enumerate(zip(asv.kingdom, asv.ASV))}
    with plt.rc_context(RC):
        fig, ax = plt.subplots(figsize=(7.2, 0.16 * len(asv) + 2.0))
        for k, s in enumerate(SOILCOL):
            g = d[d.soil == s]
            ax.scatter([nb.index(w) + (k - 0.5) * 0.32 for w in g.neighbour], [ypos[(a, b)] for a, b in zip(g.kingdom, g.ASV)],
                       s=14, color=SOILCOL[s], lw=0, label=s, zorder=3)
        for j in range(1, len(nb)):
            if fam[nb[j]] != fam[nb[j - 1]]: ax.axvline(j - 0.5, color="#dddddd", lw=0.5, zorder=0)
        kb = (asv.kingdom == "bacteria").sum(); ax.axhline(kb - 0.5, color="black", lw=0.6)
        ax.set_yticks(range(len(asv))); ax.set_yticklabels(lab, style="italic"); ax.invert_yaxis()
        ax.set_xticks(range(len(nb))); ax.set_xticklabels(nb, rotation=80, style="italic")
        prev = None
        for j, w in enumerate(nb):
            if fam[w] != prev:
                k_ = sum(1 for x in nb if fam[x] == fam[w]); ax.text(j - 0.4 + (k_ - 1) / 2 + 0.4, -1.2, fam[w][:4] + ".", ha="center", fontsize=5, color="#555555", rotation=45)
                prev = fam[w]
        ax.set_xlim(-0.6, len(nb) - 0.4)
        ax.text(len(nb) - 0.3, kb / 2, "Bacteria", rotation=-90, va="center", fontsize=7, fontweight="bold")
        ax.text(len(nb) - 0.3, kb + (len(asv) - kb) / 2, "Fungi", rotation=-90, va="center", fontsize=7, fontweight="bold")
        h = [Line2D([], [], marker="o", ls="", color=c, label=s, markersize=4) for s, c in SOILCOL.items()]
        ax.legend(handles=h, frameon=False, loc="upper left", bbox_to_anchor=(1.03, 1.0))
        return _save(fig, out), d


# ---------------- figures carried over (Figs S8, S14, S16; soil texture)
SYN = {'Centaurea cyanus': 'Cyanus segetum', 'Bromus sterilis': 'Anisantha sterilis', 'Festuca arundinacea': 'Lolium arundinaceum'}

def _family_order(meta):
    f = meta[~meta.species.isin(CTRL)].groupby("species").family.first()
    return f.sort_index().sort_values(kind="stable")



def fig_glomeromycota(DF, out):
    otu = pd.read_csv(os.path.join(DF, "otu_fungi_table_R_647ASVs.csv"), sep="\t", index_col=0).T.astype(float)
    tax = pd.read_csv(os.path.join(DF, "tax_fungi_R_647ASVs.csv"), sep=None, engine="python", index_col=0)
    meta = pd.read_csv(os.path.join(DF, "metadata_fungi_647ASVs.csv")).set_index("sample")
    meta = meta.loc[meta.index.intersection(otu.index)].copy()
    glom = [c for c in otu.columns if str(tax.Phylum.get(c, "")).startswith("Glomeromycota")]
    X = otu.loc[meta.index]
    meta["glom"] = 100 * X[glom].sum(1) / X.sum(1)
    meta["soilL"] = meta.soil.map(SOILCODE)
    meta["fam"] = np.where(meta.species.isin(CTRL), "Controls (wheat)", meta.family)
    fams = sorted(meta.loc[~meta.species.isin(CTRL), "fam"].unique()) + ["Controls (wheat)"]
    with plt.rc_context(RC):
        fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw=dict(width_ratios=[1, 5.2], wspace=0.12), sharey=True)
        rng = np.random.default_rng(1)
        for i, s in enumerate(SOILCOL):
            v = meta.loc[meta.soilL == s, "glom"]
            a.boxplot(v, positions=[i], widths=0.55, showfliers=False, medianprops=dict(color="black"),
                      boxprops=dict(color=SOILCOL[s]), whiskerprops=dict(color=SOILCOL[s]), capprops=dict(color=SOILCOL[s]))
            a.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=3, color=SOILCOL[s], alpha=0.5, lw=0)
        a.set_xticks([0, 1]); a.set_xticklabels(["Cereal-\ndominated", "Mixed-\nsystem"])
        a.set_ylabel("Glomeromycota (% of fungal reads)"); a.set_title("Soil", loc="left")
        for j, f in enumerate(fams):
            for k, s in enumerate(SOILCOL):
                v = meta.loc[(meta.fam == f) & (meta.soilL == s), "glom"]
                if len(v) == 0: continue
                x = j + (k - 0.5) * 0.36
                b.boxplot(v, positions=[x], widths=0.3, showfliers=False, medianprops=dict(color="black", lw=0.8),
                          boxprops=dict(color=SOILCOL[s], lw=0.7), whiskerprops=dict(color=SOILCOL[s], lw=0.7), capprops=dict(color=SOILCOL[s], lw=0.7))
                b.scatter(x + rng.uniform(-0.08, 0.08, len(v)), v, s=2, color=SOILCOL[s], alpha=0.5, lw=0)
        b.set_xticks(range(len(fams))); b.set_xticklabels(fams, rotation=45, ha="right")
        b.set_title("Neighbour family", loc="left")
        for s, c in SOILCOL.items(): b.plot([], [], color=c, lw=4, label=s)
        b.legend(frameon=False, loc="upper right")
        a.margins(y=0.03)
        return _save(fig, out), meta.groupby("soilL").glom.agg(["mean", "median", "size"]).round(2)



def fig_trait_F(trait_comp_csv, out):
    comp = pd.read_csv(trait_comp_csv)
    lab = {"Distance": "Phylogenetic distance to wheat", "Segetal": "Segetal status", "Nitrogen": "Ellenberg N", "Light": "Ellenberg L",
           "Btot": "Neighbour biomass (Btot)", "SRL": "Specific root length", "RDMC": "Root dry matter content", "RootShoot": "Root:shoot ratio"}
    traits = list(lab)[::-1]
    cells = [("bacteria", "Cereal-dominated (AF193)"), ("bacteria", "Mixed-system (AF190)"), ("fungi", "Cereal-dominated (AF193)"), ("fungi", "Mixed-system (AF190)")]
    with plt.rc_context(RC):
        fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.6), sharey=True, sharex=True, gridspec_kw=dict(wspace=0.06))
        for ax, (kg, s) in zip(axes, cells):
            g = comp[(comp.kingdom == kg) & (comp.soil == s)].set_index("trait")
            for y, t in enumerate(traits):
                r = g.loc[t]; P = r.P_pot if t == "Btot" else r.P_species
                if P < 0.05:   mk = dict(marker="o", s=26, color=SOILCOL[s], edgecolor="black", lw=0.5)
                elif r.P_pot < 0.05: mk = dict(marker="o", s=26, facecolor="white", edgecolor=SOILCOL[s], lw=1.0)
                else:          mk = dict(marker="o", s=10, color="#bbbbbb", lw=0)
                ax.scatter(r.F, y, zorder=3, **mk)
                ax.plot([1, r.F], [y, y], color="#cccccc", lw=0.6, zorder=1)
            ax.axvline(1, color="black", lw=0.5, ls=":")
            ax.set_title(f"{kg.capitalize()}\n{s}", loc="left"); ax.set_xlabel("db-RDA pseudo-F")
        axes[0].set_yticks(range(len(traits))); axes[0].set_yticklabels([lab[t] for t in traits])
        axes[3].scatter([], [], marker="o", s=26, color="#777777", edgecolor="black", lw=0.5, label="P < 0.05, species-level permutation")
        axes[3].scatter([], [], marker="o", s=26, facecolor="white", edgecolor="#777777", lw=1.0, label="P < 0.05 only if pots permuted")
        axes[3].scatter([], [], marker="o", s=10, color="#bbbbbb", label="P ≥ 0.05")
        fig.legend(frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12))
        for ax in axes: ax.margins(x=0.08)
        return _save(fig, out)



def fig_meqo_performance(s7_csv, out):
    t = pd.read_csv(s7_csv)
    cols = [("R2_insample", "In-\nsample", "#8B2E2E"), ("null_R2_mean", "Shuffled-\nbiomass\nsearch", "#9aa3ad"),
            ("R2_cv", "Cross-\nvalidated", "#4a78b8"), ("R2_oos", "Out-of-\nsample", "#3e6b3e")]
    with plt.rc_context(RC):
        fig, axes = plt.subplots(2, 2, figsize=(5.6, 4.6), sharey=True, gridspec_kw=dict(wspace=0.08, hspace=0.55))
        axes = axes.ravel()
        for ax, (_, r) in zip(axes, t.iterrows()):
            for i, (c, l, col) in enumerate(cols):
                v = max(r[c], -0.29)
                ax.bar(i, v, color=col, width=0.7)
                if r[c] < -0.29:
                    ax.text(i, -0.28, f"{r[c]:.2f}", ha="center", va="bottom", fontsize=6, color="white")
                else:
                    ax.text(i, v + (0.012 if v >= 0 else -0.012), f"{r[c]:.2f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=6)
            ax.plot([0.62, 1.38], [r.null_R2_p95] * 2, color="black", lw=0.8)
            ax.text(1.0, r.null_R2_p95 + 0.01, "95th pct", fontsize=5.5, ha="center", va="bottom")
            ax.axhline(0, color="black", lw=0.5)
            ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l, _ in cols], fontsize=5.5)
            ax.set_title(f"{r.kingdom.capitalize()}, {r.soil}\nP = {r.P_null:.3f} vs shuffled biomass ({int(r.n_perm)} runs)", loc="left")
        for ax in axes[::2]: ax.set_ylabel("Wheat biomass variance\nexplained (R²)")
        axes[0].set_ylim(-0.3, 0.36)
        return _save(fig, out)




def fig_assembly_bacteria(NST_DIR, DF, out):
    """Fig. S7: bacterial assembly metrics per neighbour treatment and soil: (a) mean betaNTI, (b) taxonomic and phylogenetic NST.
    Fungal values are not shown: the ITS1 tree cannot be aligned reliably across fungal phyla (Methods S2)."""
    d = nst_table(NST_DIR); d = d[d.kingdom == "bacteria"]
    b = pd.read_csv(os.path.join(NST_DIR, "bnti_bacteria.csv"))
    b["Species"] = b.Species.str.replace("_", " ").str.replace("bursa pastoris", "bursa-pastoris").str.replace("speculum veneris", "speculum-veneris").str.replace("foenum graecum", "foenum-graecum")
    b["Species"] = b.Species.replace(NST_SYN); b["soil"] = b.Soil.map(SOILCODE)
    sp, fam = species_order(_meta(DF)); order = sp + CTRL
    lab = [s if s not in CTRL else ("wheat × wheat" if s == "wheat x wheat" else "wheat alone") for s in order]
    vals = d.merge(b[["Species", "soil", "mean_bNTI"]], on=["Species", "soil"], how="outer")
    vals = vals[vals.Species.isin(order)]; vals.to_csv(out + "_values.csv", index=False)
    with plt.rc_context(RC):
        fig, axes = plt.subplots(1, 3, figsize=(7.2, 4.8), sharey=True, gridspec_kw=dict(wspace=0.14))
        for ax, (col, ttl, xl) in zip(axes, (("mean_bNTI", "(a) Mean βNTI", "βNTI"), ("NST", "(b) Taxonomic NST", "NST"), ("pNST", "(c) Phylogenetic NST", "NST"))):
            for k, s in enumerate(SOILCOL):
                gs = vals[vals.soil == s].set_index("Species")[col].reindex(order)
                ax.scatter(gs.values, np.arange(len(order)) + (k - 0.5) * 0.3, s=9, color=SOILCOL[s], lw=0, zorder=3)
            if col == "mean_bNTI":
                for x in (-2, 2): ax.axvline(x, color="black", lw=0.6, ls="--")
                ax.set_xlim(-3, 4)
            else:
                ax.axvline(0.5, color="black", lw=0.6, ls="--"); ax.set_xlim(-0.02, 1.02); ax.set_xticks([0, 0.5, 1])
            ax.set_title(ttl, loc="left"); ax.set_xlabel(xl)
            for j in range(1, len(order)):
                if fam.get(order[j], "Wheat") != fam.get(order[j - 1], "Wheat"): ax.axhline(j - 0.5, color="#e5e5e5", lw=0.5, zorder=0)
        axes[0].set_yticks(range(len(order))); axes[0].set_yticklabels(lab, style="italic"); axes[0].invert_yaxis()
        fig.legend(handles=[Line2D([], [], marker="o", ls="", color=c, label=s, markersize=4) for s, c in SOILCOL.items()],
                   frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.05))
        return _save(fig, out)
