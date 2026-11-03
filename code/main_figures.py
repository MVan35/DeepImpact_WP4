# -*- coding: utf-8 -*-
"""Main-text figures (Times New Roman). Called from Paper_figures.ipynb.
  fig1(DF, out)                    PCoA, variance partition (sequential PERMANOVA), phylum composition
  fig2(core_csv, DF, out)          ASV-level FDR core 
  fig3(DF, out)                    wheat biomass by neighbour and soil
  fig4(cf, assoc, module_dir, res_dir, out_dir)   convergence_figures.fig4_ranking_and_neighbour in Times
All functions write PNG (400 dpi), PDF and SVG (text kept as text). Soil colours as everywhere in the paper."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Patch
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu, chi2

SOILCODE = {"resistant": "Cereal-dominated (AF193)", "Permissive": "Mixed-system (AF190)"}
SOILCOL = {"Cereal-dominated (AF193)": "#3B75AF", "Mixed-system (AF190)": "#E08A2C"}
KCOL = {"bacteria": "#3B75AF", "fungi": "#E08A2C"}
CTRL = ["wheat x wheat", "wheat alone"]
from si_figures import FAMCOL as _F21          # one family palette for Fig. 1 and Figs S8, S12
FAMCOL = dict(_F21)
TNR = {"font.family": "serif", "font.serif": ["Times New Roman", "Times"], "mathtext.fontset": "custom",
       "mathtext.rm": "Times New Roman", "mathtext.it": "Times New Roman:italic", "font.size": 8, "axes.labelsize": 8,
       "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False,
       "axes.spines.right": False, "axes.linewidth": 0.6, "pdf.fonttype": 42, "svg.fonttype": "none"}


def _save(fig, out):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{out}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)
    return out + ".png"


def _load(DF, kg):
    n = 1314 if kg == "bacteria" else 647
    X = pd.read_csv(os.path.join(DF, f"otu_{kg}_table_R_{n}ASVs.csv"), sep="\t", index_col=0).T.astype(float)
    md = pd.read_csv(os.path.join(DF, f"metadata_{kg}_{n}ASVs.csv")).set_index("sample")
    md = md.loc[md.index.intersection(X.index)]
    tax = pd.read_csv(os.path.join(DF, f"tax_{kg}_R_{n}ASVs.csv"), sep=None, engine="python", index_col=0)
    return X.loc[md.index], md, tax


def pcoa(D):
    n = D.shape[0]; J = np.eye(n) - 1 / n; B = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(B); o = np.argsort(w)[::-1]; w, V = w[o], V[:, o]
    pos = w > 0
    return V[:, :2] * np.sqrt(w[:2]), 100 * w[:2] / w[pos].sum()


def _ellipse(ax, xy, col, q=0.95):
    mu = xy.mean(0); C = np.cov(xy.T); vals, vecs = np.linalg.eigh(C)
    k = np.sqrt(chi2.ppf(q, 2)); ang = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))
    ax.add_patch(Ellipse(mu, 2 * k * np.sqrt(vals[1]), 2 * k * np.sqrt(vals[0]), angle=ang, fc=col, ec=col, alpha=0.10, lw=0))
    ax.add_patch(Ellipse(mu, 2 * k * np.sqrt(vals[1]), 2 * k * np.sqrt(vals[0]), angle=ang, fc="none", ec=col, lw=1.0))


def table1_from_docx(path):
    """Read Table 1 (sequential PERMANOVA) from the manuscript: {kingdom: {soil|species|inter: (R2 %, F)}}."""
    from docx import Document
    tb = next(t for t in Document(path).tables if t.rows[0].cells[0].text.strip() == "Model")
    out, kg = {}, None
    key = {"Soil": "soil", "Species": "species", "Soil × Species": "inter"}
    for r in tb.rows[1:]:
        c = [x.text.strip() for x in r.cells]
        if c[0]: kg = "bacteria" if c[0].startswith("Bacterial") else "fungi"
        if c[1] in key: out.setdefault(kg, {})[key[c[1]]] = (float(c[3]), float(c[4]))
    return out


def table1_from_csv(path):
    """Read Table 1 as written by Paper_analyses.qmd (Table1_variance_partition.csv): {kingdom: {soil|species|inter: (R2 %, F)}}."""
    t = pd.read_csv(path)
    def key(term):                      # robust to how R wrote the multiplication sign in a non-UTF-8 locale
        term = str(term)
        if term.startswith("Soil of origin"): return "soil"
        if term.startswith("Neighbour identity"): return "species"
        if term.startswith("Soil") and "neighbour" in term: return "inter"
        return None
    out = {}
    for _, r in t.iterrows():
        k = key(r.Term)
        if k is None: continue
        kg = "bacteria" if r.Kingdom.startswith("Bacteria") else "fungi"
        out.setdefault(kg, {})[k] = (float(r.R2_pct), float(r.F))
    return out


def fig1(DF, out, table1):
    """table1: {kingdom: {'soil': (R2%, F), 'species': (R2%, F), 'inter': (R2%, F)}}, e.g. table1_from_docx(manuscript)."""
    with plt.rc_context(TNR):
        fig = plt.figure(figsize=(7.2, 6.4))
        gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.95], height_ratios=[1.15, 1], wspace=0.34, hspace=0.42)
        famlab = {}
        for j, kg in enumerate(("bacteria", "fungi")):
            X, md, _ = _load(DF, kg)
            D = squareform(pdist(X.values, "braycurtis")); Y, ev = pcoa(D)
            ax = fig.add_subplot(gs[0, j]); md = md.assign(PC1=Y[:, 0], PC2=Y[:, 1], soilL=md.soil.map(SOILCODE))
            for s, c in SOILCOL.items():
                g = md[md.soilL == s]; ax.scatter(g.PC1, g.PC2, s=2.5, color=c, alpha=0.35, lw=0, zorder=1)
                _ellipse(ax, g[["PC1", "PC2"]].values, c)
            md["fam"] = np.where(md.species.isin(CTRL), "Wheat", md.family)
            cen = md.groupby(["soilL", "species", "fam"])[["PC1", "PC2"]].mean().reset_index()
            for (s, mk) in zip(SOILCOL, ("o", "^")):
                g = cen[cen.soilL == s]
                ax.scatter(g.PC1, g.PC2, s=16, marker=mk, c=[FAMCOL[f] for f in g.fam], edgecolor="black", lw=0.3, zorder=3)
            famlab.update({f: FAMCOL[f] for f in cen.fam.unique()})
            ax.set_xlabel(f"PCoA axis 1 ({ev[0]:.1f}%)"); ax.set_ylabel(f"PCoA axis 2 ({ev[1]:.1f}%)")
            ax.set_title(("Bacteria (rpoB)" if kg == "bacteria" else "Fungi (ITS1)"), loc="left", fontweight="bold")
            t = table1[kg]
            ax.text(0.02, 0.99, f"Soil R² = {t['soil'][0]:.1f}%\nSpecies R² = {t['species'][0]:.1f}%\nSoil × species R² = {t['inter'][0]:.1f}%",
                    transform=ax.transAxes, fontsize=6, va="top", style="italic", linespacing=1.2)
            ax.margins(y=0.12)
        axl = fig.add_subplot(gs[0, 2]); axl.axis("off")
        h = [Patch(fc=c, ec=c, alpha=0.4, label=s.replace(" (", "\n(")) for s, c in SOILCOL.items()]
        h += [Line2D([], [], ls="", marker="o", mfc="white", mec="black", label="centroid, cereal-dominated"),
              Line2D([], [], ls="", marker="^", mfc="white", mec="black", label="centroid, mixed-system")]
        h += [Line2D([], [], ls="", marker="s", mfc=famlab[f], mec="black", mew=0.3, label=f) for f in sorted(famlab, key=lambda x: (x == "Wheat", x))]
        h[-1].set_label("Wheat (controls)")
        axl.legend(handles=h, frameon=False, loc="upper left", fontsize=6.5, title="Soil and neighbour family", title_fontsize=7)
        # (b) variance partition
        axb = fig.add_subplot(gs[1, 0:2])
        terms = [("soil", "Soil"), ("species", "Species"), ("inter", "Soil × species")]
        for i, (k, lab) in enumerate(terms):
            for dy, kg in ((-0.18, "bacteria"), (0.18, "fungi")):
                r2, F = table1[kg][k]
                axb.barh(i + dy, r2, height=0.34, color=KCOL[kg], lw=0)
                axb.text(r2 + 0.15, i + dy, f"{r2:.1f}%  (F = {F:.2f})", va="center", fontsize=6.5)
        axb.set_yticks(range(3)); axb.set_yticklabels([l for _, l in terms]); axb.invert_yaxis()
        axb.set_xlim(0, 14.5); axb.set_xlabel("Variance explained (R², %)")
        axb.legend(handles=[Patch(color=KCOL["bacteria"], label="Bacteria"), Patch(color=KCOL["fungi"], label="Fungi")], frameon=False, loc="upper right")
        axb.set_title("Sequential PERMANOVA (Bray–Curtis, 999 permutations)", loc="left")
        # (c, d) phylum composition, all samples pooled
        sub = gs[1, 2].subgridspec(2, 1, hspace=0.95)
        for r_, kg in enumerate(("fungi", "bacteria")):
            X, md, tax = _load(DF, kg)
            ph = tax.Phylum.reindex(X.columns).fillna("unidentified").astype(str)
            n_asv = ph.value_counts(); ra = (X.div(X.sum(1), axis=0).mean()).groupby(ph).sum()
            top = list(ra.sort_values(ascending=False).index[:7])
            n_asv = n_asv.groupby(lambda p: p if p in top else "other").sum(); ra = ra.groupby(lambda p: p if p in top else "other").sum()
            order = [p for p in ra.sort_values(ascending=False).index if p != "other"] + (["other"] if "other" in ra.index else [])
            cmap = plt.get_cmap("tab10"); colp = {p: (cmap(i) if p != "other" else "#cccccc") for i, p in enumerate(order)}
            axc = fig.add_subplot(sub[r_]); b1 = b2 = 0
            for p in order:
                v1 = n_asv.get(p, 0) / n_asv.sum(); v2 = ra.get(p, 0) / ra.sum()
                axc.bar(0, v1, bottom=b1, color=colp[p], width=0.7, lw=0); axc.bar(1, v2, bottom=b2, color=colp[p], width=0.7, lw=0, label=p)
                b1 += v1; b2 += v2
            axc.set_xticks([0, 1]); axc.set_xticklabels([f"ASVs\n(n = {n_asv.sum():,})", "Relative\nabundance"], fontsize=6.5)
            axc.set_ylim(0, 1); axc.set_yticks([0, 0.5, 1]); axc.set_ylabel("Proportion"); axc.set_title(("Fungi" if kg == "fungi" else "Bacteria"), loc="left", fontweight="bold", pad=3)
            axc.text(-0.34, 1.10, "(c)" if kg == "fungi" else "(d)", transform=axc.transAxes, fontsize=10, fontweight="bold", va="bottom")
            axc.legend(frameon=False, fontsize=5.8, loc="center left", bbox_to_anchor=(1.0, 0.5), handlelength=1)
        fig.axes[0].text(-0.28, 1.06, "(a)", transform=fig.axes[0].transAxes, fontsize=10, fontweight="bold", va="bottom")
        axb.text(-0.14, 1.06, "(b)", transform=axb.transAxes, fontsize=10, fontweight="bold", va="bottom")
        return _save(fig, out)


def fig2(P, lfc, T, fam, out):
    """ASV x neighbour calls passing >= 2 of 3 methods at FDR-adjusted q < 0.05 (log2FC > 2).
    P: {(kingdom, soil): set of (neighbour, ASV)}; soil keys 'Resistant'/'Permissive' as in da_consensus;
    lfc: {(kingdom, soil, neighbour, ASV): DESeq2 log2FC}; T[kingdom]['tax']: taxonomy; fam: {neighbour: family}."""
    rows = [dict(kingdom=kg, soil=SOILCODE["resistant" if s == "Resistant" else "Permissive"], neighbour=w, ASV=c,
                 l=lfc.get((kg, s, w, c), np.nan)) for (kg, s), pp in P.items() for w, c in pp]
    d = pd.DataFrame(rows)
    d["genus"] = [T[k]["tax"].Genus.get(a) for k, a in zip(d.kingdom, d.ASV)]
    d["phylum"] = [T[k]["tax"].Phylum.get(a) for k, a in zip(d.kingdom, d.ASV)]
    d.sort_values(["kingdom", "soil", "neighbour", "ASV"]).to_csv(out + "_values.csv", index=False)
    nb = sorted(fam, key=lambda s: (fam[s], s))
    asv = d.groupby(["kingdom", "phylum", "genus", "ASV"]).neighbour.nunique().reset_index() \
           .sort_values(["kingdom", "phylum", "genus", "neighbour"], ascending=[True, True, True, False])
    lab = [f"{g} {a.replace('Cluster_', 'ASV_')}" for g, a in zip(asv.genus, asv.ASV)]
    ypos = {(k, a): i for i, (k, a) in enumerate(zip(asv.kingdom, asv.ASV))}
    with plt.rc_context(TNR):
        fig, ax = plt.subplots(figsize=(7.4, 0.205 * len(asv) + 2.6))
        for k, s in enumerate(SOILCOL):
            g = d[d.soil == s]
            ax.scatter([nb.index(w) for w in g.neighbour], [ypos[(a, b)] for a, b in zip(g.kingdom, g.ASV)],
                       s=30, color=SOILCOL[s], lw=0.3, edgecolor="white", label=s, zorder=3)
        for i in range(len(asv)):
            ax.axhline(i, color="#f0f0f0", lw=0.4, zorder=0)
        bounds, prev = [], None
        for j, w in enumerate(nb):
            if fam[w] != prev: bounds.append(j); prev = fam[w]
        bounds.append(len(nb))
        for j in bounds[1:-1]: ax.axvline(j - 0.5, color="#bbbbbb", lw=0.6, zorder=1)
        kb = (asv.kingdom == "bacteria").sum(); ax.axhline(kb - 0.5, color="black", lw=0.8)
        ax.set_yticks(range(len(asv))); ax.set_yticklabels(lab, style="italic", fontsize=8); ax.invert_yaxis()
        ax.set_xticks(range(len(nb))); ax.set_xticklabels(nb, rotation=75, style="italic", fontsize=8, ha="right", rotation_mode="anchor")
        ax.set_xlim(-0.6, len(nb) - 0.4); ax.set_ylim(len(asv) - 0.4, -0.6)
        for a0, a1 in zip(bounds[:-1], bounds[1:]):
            ax.text((a0 + a1 - 1) / 2, -0.9, fam[nb[a0]], ha="left", va="bottom", rotation=50, fontsize=8, rotation_mode="anchor")
        ax.text(len(nb) - 0.1, (kb - 1) / 2, "Bacteria", rotation=-90, va="center", fontsize=9, fontweight="bold")
        ax.text(len(nb) - 0.1, kb + (len(asv) - kb - 1) / 2, "Fungi", rotation=-90, va="center", fontsize=9, fontweight="bold")
        ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
        h = [Line2D([], [], marker="o", ls="", color=c, label=s, markersize=6) for s, c in SOILCOL.items()]
        ax.legend(handles=h, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.215), ncol=2, fontsize=8.5)
        return _save(fig, out)


def fig3(DF, out, star=None):
    """Wheat biomass by neighbour and soil. star: treatments whose wheat biomass differs between soils (two-sided
    Mann–Whitney, Benjamini–Hochberg q < 0.05 across the 30 treatments); computed here when None."""
    md = pd.read_csv(os.path.join(DF, "metadata_bacteria_1314ASVs.csv")).dropna(subset=["Dry_wheat_mass"])
    md["soilL"] = md.soil.map(SOILCODE)
    fam = md[~md.species.isin(CTRL)].groupby("species").family.first()
    sp = list(fam.reset_index().sort_values(["family", "species"]).species)
    order = sp + ["wheat alone", "wheat x wheat"]
    if star is None:   # two-sided Mann-Whitney between soils, Benjamini-Hochberg across the 30 treatments, q < 0.05
        from statsmodels.stats.multitest import multipletests
        pv = [mannwhitneyu(*[md[(md.species == s) & (md.soilL == L)].Dry_wheat_mass for L in SOILCOL], alternative="two-sided").pvalue for s in order]
        q = multipletests(pv, method="fdr_bh")[1]
        star = [s for s, qq in zip(order, q) if qq < 0.05]
    ctrl = md[md.species.isin(CTRL)].groupby("species").Dry_wheat_mass.mean()
    with plt.rc_context(TNR):
        fig, ax = plt.subplots(figsize=(7.2, 3.6))
        for i, x in enumerate(order):
            for k, s in enumerate(SOILCOL):
                v = md[(md.species == x) & (md.soilL == s)].Dry_wheat_mass.values
                bp = ax.boxplot(v, positions=[i + (k - 0.5) * 0.38], widths=0.34, patch_artist=True, showfliers=True,
                                flierprops=dict(marker="o", markersize=1.8, markerfacecolor="black", markeredgecolor="none"),
                                medianprops=dict(color="black", lw=0.7), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
                bp["boxes"][0].set(facecolor=SOILCOL[s], edgecolor="black", lw=0.4)
            if x in star:
                top = md[md.species == x].Dry_wheat_mass.max(); ax.text(i, top + 0.04, "*", ha="center", va="bottom", fontsize=11, color="#c00000")
        ax.axhline(ctrl["wheat x wheat"], color="#1f3fbf", lw=0.8, ls="--"); ax.axhline(ctrl["wheat alone"], color="#c00000", lw=0.8, ls="--")
        famseq = [fam.get(x, "Control") for x in order]
        for j in range(1, len(order)):
            if famseq[j] != famseq[j - 1]: ax.axvline(j - 0.5, color="#aaaaaa", lw=0.5, ls="--", zorder=0)
        ymax = md.Dry_wheat_mass.max() * 1.12
        for f in ("Asteraceae", "Fabaceae", "Poaceae"):
            idx = [j for j, g in enumerate(famseq) if g == f]
            ax.text(np.mean(idx), ymax * 0.97, f, ha="center", va="top", fontsize=8.5)
        ax.set_ylim(0, ymax); ax.set_xlim(-0.6, len(order) - 0.4)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([("wheat × wheat" if o == "wheat x wheat" else o) for o in order], rotation=80, style="italic", fontsize=8)
        ax.set_ylabel("Wheat aboveground dry mass (g)")
        h = [Patch(fc=c, ec="black", lw=0.4, label=s) for s, c in SOILCOL.items()] + \
            [Line2D([], [], color="#c00000", ls="--", lw=0.8, label="mean, wheat alone"), Line2D([], [], color="#1f3fbf", ls="--", lw=0.8, label="mean, wheat × wheat")]
        ax.legend(handles=h, frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.50), fontsize=8)
        return _save(fig, out), star


def fig4(scores_csv, tests_csv, out, label_ends=3):
    """Fig. 4: species-centroid db-RDA of endosphere composition on species-mean wheat biomass (28 weeds per panel).
    scores_csv / tests_csv: fig4_scores.csv and fig4_direct_test.csv written by revision_analyses.R.
    Colours follow the groups tested as covariates (monocots: Poaceae, Juncaceae; legumes; other eudicots) and avoid
    the blue/orange used for the two soils elsewhere. The `label_ends` species at each end of the biomass axis are named."""
    sc = pd.read_csv(scores_csv); tt = pd.read_csv(tests_csv).set_index(["kingdom", "soil"])
    GRP = {"Poaceae": ("Poaceae (monocot)", "#1B7837"), "Juncaceae": ("Juncaceae (monocot)", "#A6DBA0"),
           "Fabaceae": ("Fabaceae (legume)", "#9C3D8F")}
    OTHER = ("other eudicots", "#E6E6E6")
    abbr = lambda s: s.split()[0][0] + ". " + " ".join(s.split()[1:])
    with plt.rc_context(TNR):
        fig, axs = plt.subplots(2, 2, figsize=(7.0, 6.4))
        for k, (ax, (kg, so)) in enumerate(zip(axs.flat, [(kg, so) for kg in ("bacteria", "fungi") for so in SOILCOL])):
            d = sc[(sc.kingdom == kg) & (sc.soil == so)].copy(); r = tt.loc[(kg, so)]
            if np.corrcoef(d.CAP1, d.biomass)[0, 1] < 0: d["CAP1"] = -d.CAP1          # biomass increases to the right
            b = d.biomass; size = 15 + 160 * (b - b.min()) / (b.max() - b.min())
            c = [GRP.get(f, OTHER)[1] for f in d.family]
            order = np.argsort(-size.values)                                            # large points first, small on top
            ax.scatter(d.CAP1.values[order], d.MDS1.values[order], s=size.values[order], c=[c[i] for i in order],
                       edgecolor="black", lw=0.4, zorder=3)
            ax.axhline(0, color="#dddddd", lw=0.5, zorder=0); ax.axvline(0, color="#dddddd", lw=0.5, zorder=0)
            if label_ends:
                ds = d.sort_values("CAP1"); yr = d.MDS1.max() - d.MDS1.min()
                for end, sub, ha, dx in (("low", ds.head(label_ends), "right", -5), ("high", ds.tail(label_ends), "left", 5)):
                    sub = sub.sort_values("MDS1"); last = None; dy = 0
                    for _, row in sub.iterrows():
                        dy = (dy + 7 if (last is not None and row.MDS1 - last < 0.07 * yr) else 0); last = row.MDS1
                        rad = np.sqrt(15 + 160 * (row.biomass - b.min()) / (b.max() - b.min())) / 2     # marker radius, points
                        ax.annotate(abbr(row.species), (row.CAP1, row.MDS1), xytext=(np.sign(dx) * (rad + 2), dy), textcoords="offset points",
                                    ha=ha, va="center", fontsize=6, style="italic", color="#333333", zorder=4)
            xr = d.CAP1.max() - d.CAP1.min(); ax.set_xlim(d.CAP1.min() - 0.42 * xr, d.CAP1.max() + 0.40 * xr); ax.margins(y=0.10)
            ax.set_xlabel(f"db-RDA1, wheat biomass ({d.eig_CAP1_pct.iloc[0]:.1f}%)")
            ax.set_ylabel(f"MDS1, unconstrained ({d.eig_MDS1_pct.iloc[0]:.1f}%)")
            ttl = f"{'Bacteria' if kg == 'bacteria' else 'Fungi'}, {'cereal-dominated' if so.startswith('Cereal') else 'mixed-system'} soil"
            pstr = f"P = {r.P:.3f}" if r.P >= 0.001 else "P < 0.001"
            ax.set_title(f"({'abcd'[k]}) {ttl}\nR² = {r.R2_pct:.1f}%, {pstr}, q = {r.q_BH:.2f}", loc="left")
        h = [Line2D([], [], ls="", marker="o", mfc=v[1], mec="black", mew=0.4, markersize=6, label=v[0]) for v in list(GRP.values()) + [OTHER]]
        fig.legend(handles=h, loc="lower center", ncol=4, frameon=False,
                   title="Neighbour group; point size, species-mean wheat biomass", bbox_to_anchor=(0.5, -0.01))
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        return _save(fig, out)
