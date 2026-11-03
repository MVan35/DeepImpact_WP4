"""
Fig. S15 -- weed effects on microbial alpha-diversity.

Observed richness and Shannon diversity per pot, summarised per neighbour
species and soil, for both kingdoms. Produces Fig. S15 and the accompanying
two-way ANOVA table.

Inputs (paper_1/data/matrice_phase1/df/):
  otu_bacteria_table_R_1314ASVs.csv   TAB-separated despite the .csv extension
  otu_fungi_table_R_647ASVs.csv       same
  metadata_bacteria_1314ASVs.csv
  metadata_fungi_647ASVs.csv

NOMENCLATURE. The bacterial metadata uses three superseded names. They are
harmonised on load to the accepted names used in Table S1 (which was built from
the fungal metadata -- verified: all 30 per-species replicate counts match
exactly), so both panels label the same plant the same way:
    Centaurea cyanus    -> Cyanus segetum
    Festuca arundinacea -> Lolium arundinaceum
    Bromus sterilis     -> Anisantha sterilis
The old names remain correct as genus-level synonyms and Table S1 keeps them in
its Genus column. Tables S3 and S4 still carry the superseded binomials and
need the same substitution.
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

DATA_DIR = os.path.join("paper_1", "data", "matrice_phase1", "df")
OUT_DIR = "results_meqo_convergence"

C_CER = "#B5651D"      # AF193 cereal-dominated  -- same pair as Fig. 4b
C_MIX = "#2E6FA7"      # AF190 mixed-system
META_GREY = "#8B949D"
CONTROLS = {"wheat alone", "wheat x wheat", "Wheat Alone", "Wheat x Wheat"}
SOIL_LABEL = {"resistant": "cereal", "Permissive": "mixed"}

# resistant = cereal-dominated (AF193); Permissive = mixed-system (AF190).
# Verified against Table S1: pairing Permissive with mixed-system matches all 30
# per-species replicate counts exactly, the alternative matches 5 of 30. Note
# results_v4/trait_microbiome_effect_sizes.csv has this SWAPPED in its
# soil_label column -- read `soil`, never `soil_label`.
AF = {"resistant": "AF193", "Permissive": "AF190"}

# Superseded binomials in the bacterial metadata -> the accepted names in Table S1
SYNONYMS = {"Centaurea cyanus": "Cyanus segetum",
            "Festuca arundinacea": "Lolium arundinaceum",
            "Bromus sterilis": "Anisantha sterilis"}


def _frame(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.8)
        ax.spines[sp].set_color("#3A3B38")


def load_alpha(data_dir=DATA_DIR):
    """Per-pot richness and Shannon for both kingdoms.

    sep=None + engine='python' sniffs the delimiter, so this works whether the
    OTU tables are tab- or comma-separated (they are tab-separated, despite the
    extension -- reading them with the default comma gives a single column and
    a silently degenerate result).
    """
    out = {}
    for kingdom, otu_f, meta_f in [
            ("bacteria", "otu_bacteria_table_R_1314ASVs.csv",
             "metadata_bacteria_1314ASVs.csv"),
            ("fungi", "otu_fungi_table_R_647ASVs.csv",
             "metadata_fungi_647ASVs.csv")]:
        otu = pd.read_csv(os.path.join(data_dir, otu_f), sep=None,
                          engine="python", index_col=0)
        meta = pd.read_csv(os.path.join(data_dir, meta_f), sep=None,
                           engine="python")
        samples = [s for s in meta["sample"] if s in otu.columns]
        M = otu[samples].T
        d = meta.set_index("sample").loc[samples].copy()
        d["richness"] = (M > 0).sum(1).values
        p = M.div(M.sum(1), axis=0).values
        with np.errstate(divide="ignore", invalid="ignore"):
            d["shannon"] = -np.nansum(np.where(p > 0, p * np.log(p), 0), axis=1)
        d["species"] = d["species"].replace(SYNONYMS)   # harmonise to Table S1
        d["soil2"] = d["soil"].map(SOIL_LABEL)
        d["is_control"] = d["species"].isin(CONTROLS)
        out[kingdom] = d.reset_index()
    return out


def alpha_anova(d, resp="richness"):
    """Two-way ANOVA with omega-squared effect sizes.

    Sequential (Type I) sums of squares, which is what the manuscript reports
    for the PERMANOVA and is kept consistent here. omega2 rather than eta2
    because it is the less biased estimator for a fixed-effects design and is
    what the text quotes.
    """
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    m = ols(f"{resp} ~ C(soil) * C(species)", data=d).fit()
    a = sm.stats.anova_lm(m, typ=1)
    ms_resid = a.loc["Residual", "sum_sq"] / a.loc["Residual", "df"]
    sst = a["sum_sq"].sum()
    a["omega2"] = [(r.sum_sq - r.df * ms_resid) / (sst + ms_resid)
                   if i != "Residual" else np.nan for i, r in a.iterrows()]
    return a


def species_means(d, resp="richness"):
    g = (d.groupby(["species", "soil2"])[resp]
           .agg(["mean", "std", "size"]).reset_index())
    g["se"] = g["std"] / np.sqrt(g["size"])
    return g


def fig_alpha_diversity(alpha, resp="richness", out_dir=OUT_DIR,
                        fname="FigS15_alpha_diversity.png",
                        width_in=7.2, height_in=7.4, scale=1.0, dpi=400):
    """Richness per neighbour species and soil, both kingdoms.

    Paired dot-and-interval rather than boxes: the quantity of interest is the
    species mean (n = 6-10 pots each), and 30 species x 2 soils x 2 kingdoms of
    boxes is unreadable at column width. The interval is +/- 1 SE of that mean;
    per-species n goes in the caption rather than onto 60 marks.

    Species are ordered by their mixed-system mean, so the soil crossover reads
    as points swapping order rather than as a colour wash. The two controls are
    drawn as reference lines, not as categories, because they are what the weed
    treatments are being compared against -- and a control sitting in the
    category axis would invite reading it as one more neighbour.

    hspace is large because panel (a)'s rotated species labels are ~90 pt tall
    and would otherwise collide with panel (b)'s title.
    """
    s = scale
    FS = {"base": 9 * s, "ann": 8 * s, "tick": 7 * s, "letter": 11 * s}
    fig, axes = plt.subplots(2, 1, figsize=(width_in, height_in),
                             gridspec_kw=dict(hspace=0.78, left=0.115,
                                              right=0.995, top=0.955,
                                              bottom=0.155))
    ylab = ("Observed richness (ASVs)" if resp == "richness"
            else "Shannon diversity")
    for ax, (kingdom, nm, marker) in zip(
            axes, [("bacteria", "Bacteria", "rpoB"), ("fungi", "Fungi", "ITS1")]):
        g = species_means(alpha[kingdom], resp)
        weeds = g[~g["species"].isin(CONTROLS)]
        order = (weeds[weeds.soil2 == "mixed"].set_index("species")["mean"]
                 .reindex(weeds["species"].unique()).sort_values().index.tolist())
        xpos = {sp: i for i, sp in enumerate(order)}
        for ctrl, ls in [("wheat alone", "--"), ("wheat x wheat", ":")]:
            for so, c in [("cereal", C_CER), ("mixed", C_MIX)]:
                row = g[(g["species"] == ctrl) & (g.soil2 == so)]
                if len(row):
                    ax.axhline(row["mean"].iloc[0], color=c, lw=0.8 * s,
                               ls=ls, alpha=0.55, zorder=1)
        for so, c, lab in [("cereal", C_CER, "Cereal-dominated (AF193)"),
                           ("mixed", C_MIX, "Mixed-system (AF190)")]:
            sub = (g[g.soil2 == so].set_index("species")
                   .reindex(order).reset_index())
            ax.errorbar([xpos[sp] for sp in sub["species"]], sub["mean"],
                        yerr=sub["se"], fmt="o", ms=3.4 * s, lw=0,
                        elinewidth=0.9 * s, capsize=0, color=c, zorder=3,
                        label=lab)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=90, fontsize=FS["tick"],
                           style="italic")
        ax.set_ylabel(ylab, fontsize=FS["base"])
        ax.set_title(f"{nm} ({marker}): neighbour identity shifts richness, "
                     f"and the shift differs by soil",
                     fontsize=FS["base"], loc="left", pad=6 * s)
        ax.margins(x=0.02, y=0.06)
        _frame(ax)
        ax.tick_params(labelsize=FS["tick"])
    axes[0].legend(frameon=False, fontsize=FS["ann"], loc="upper left",
                   handlelength=0.9, handletextpad=0.45, borderpad=0)
    axes[0].annotate("dashed: wheat alone   dotted: wheat \u00d7 wheat",
                     xy=(0.985, 0.03), xycoords="axes fraction", ha="right",
                     va="bottom", fontsize=FS["ann"] - 0.6, color=META_GREY)
    for ax, L in zip(axes, "ab"):
        bb = ax.get_position()
        fig.text(max(0.004, bb.x0 - 0.085), bb.y1 + 0.028, L,
                 fontsize=FS["letter"], fontweight="bold", ha="left", va="top")
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, fname), dpi=dpi, bbox_inches="tight",
                facecolor="white")
    return fig


def alpha_table(alpha, out_dir=OUT_DIR, fname="alpha_diversity_anova.csv"):
    """ANOVA terms for both kingdoms and both responses, as one tidy frame."""
    rows = []
    for kingdom in ("bacteria", "fungi"):
        for resp in ("richness", "shannon"):
            a = alpha_anova(alpha[kingdom], resp)
            for term, r in a.iterrows():
                rows.append(dict(kingdom=kingdom, response=resp,
                                 term=term.replace("C(", "").replace(")", ""),
                                 df=int(r["df"]), F=r.get("F"),
                                 p=r.get("PR(>F)"), omega2=r.get("omega2")))
    t = pd.DataFrame(rows)
    os.makedirs(out_dir, exist_ok=True)
    t.to_csv(os.path.join(out_dir, fname), index=False)
    return t
