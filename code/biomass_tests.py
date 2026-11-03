# -*- coding: utf-8 -*-
"""Wheat biomass statistics for Results 3.2.1 (Fig. 3, Table S3, Table S5 text), called from Paper_figures.ipynb.

  anova_tables(meta)        two-way ANOVA (soil x neighbour) on log biomass, all treatments and weeds only,
                            residual checks on raw and log scales
  between_soils(meta)       Mann-Whitney between soils for each of the 30 treatments, BH across the 30 (Fig. 3 asterisks)
  legume_test(meta)         species-mean wheat biomass under the four legumes vs the 24 other weeds
  table_s3(meta, out_xlsx)  each neighbour x soil against wheat x wheat and wheat alone: n, mean +/- SD,
                            difference with bootstrap 95% CI, two-sided Mann-Whitney P, BH q across the 56 contrasts
Soil codes: 'resistant' = cereal-dominated (AF193), 'Permissive' = mixed-system (AF190)."""
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

SOIL = {"resistant": "Cereal-dominated (AF193)", "Permissive": "Mixed-system (AF190)"}
CTRL = ["wheat x wheat", "wheat alone"]
LEGUMES = ["Lotus corniculatus", "Medicago lupulina", "Trigonella foenum-graecum", "Vicia cracca"]
SYN = {"Centaurea cyanus": "Cyanus segetum", "Festuca arundinacea": "Lolium arundinaceum", "Bromus sterilis": "Anisantha sterilis"}


def load(meta_csv):
    m = pd.read_csv(meta_csv).dropna(subset=["Dry_wheat_mass"]).copy()
    m["species"] = m.species.replace(SYN)
    m["logB"] = np.log(m.Dry_wheat_mass)
    return m


def _aov(d, y):
    f = smf.ols(f"{y} ~ C(soil) * C(species)", data=d).fit()
    a = sm.stats.anova_lm(f, typ=2)
    res = f.resid
    return a, dict(shapiro_W=stats.shapiro(res).statistic,
                   brown_forsythe_P=stats.levene(*[g[y].values for _, g in d.groupby(["soil", "species"])], center="median").pvalue)


def anova_tables(m):
    rows = []
    for set_, d in (("all treatments", m), ("weeds only", m[~m.species.isin(CTRL)])):
        for y in ("Dry_wheat_mass", "logB"):
            a, chk = _aov(d, y)
            for term in ("C(soil)", "C(species)", "C(soil):C(species)"):
                rows.append(dict(set=set_, scale="raw" if y == "Dry_wheat_mass" else "log", term=term.replace("C(", "").replace(")", ""),
                                 df=int(a.loc[term, "df"]), df_resid=int(a.loc["Residual", "df"]), F=a.loc[term, "F"], P=a.loc[term, "PR(>F)"], **chk))
    return pd.DataFrame(rows)


def between_soils(m):
    rows = []
    for sp, g in m.groupby("species"):
        a, b = g[g.soil == "resistant"].Dry_wheat_mass, g[g.soil == "Permissive"].Dry_wheat_mass
        rows.append(dict(species=sp, n_cereal=len(a), n_mixed=len(b), mean_cereal=a.mean(), mean_mixed=b.mean(),
                         P=stats.mannwhitneyu(a, b, alternative="two-sided").pvalue))
    t = pd.DataFrame(rows); t["q_BH"] = multipletests(t.P, method="fdr_bh")[1]
    return t.sort_values("P")


def legume_test(m):
    s = m[~m.species.isin(CTRL)].groupby(["soil", "species"]).Dry_wheat_mass.mean().reset_index()
    s["legume"] = s.species.isin(LEGUMES)
    rows = []
    for sc, g in s.groupby("soil"):
        rows.append(dict(soil=SOIL[sc], n_legumes=int(g.legume.sum()), n_other=int((~g.legume).sum()),
                         median_legumes=g[g.legume].Dry_wheat_mass.median(), median_other=g[~g.legume].Dry_wheat_mass.median(),
                         P=stats.mannwhitneyu(g[g.legume].Dry_wheat_mass, g[~g.legume].Dry_wheat_mass, alternative="two-sided").pvalue))
    f = smf.ols("Dry_wheat_mass ~ C(soil) + legume", data=s).fit()
    rows.append(dict(soil="both (linear model with soil)", n_legumes=int(s.legume.sum()), n_other=int((~s.legume).sum()),
                     median_legumes=s[s.legume].Dry_wheat_mass.median(), median_other=s[~s.legume].Dry_wheat_mass.median(),
                     P=f.pvalues["legume[T.True]"]))
    return pd.DataFrame(rows)


def _boot_ci(x, y, B=4999, seed=1):
    r = np.random.default_rng(seed)
    d = [r.choice(x, len(x)).mean() - r.choice(y, len(y)).mean() for _ in range(B)]
    return np.percentile(d, [2.5, 97.5])


def table_s3(m, out_xlsx=None):
    rows = []
    for sc in ("resistant", "Permissive"):
        b = m[m.soil == sc]; ww = b[b.species == "wheat x wheat"].Dry_wheat_mass.values; al = b[b.species == "wheat alone"].Dry_wheat_mass.values
        for sp, g in b.groupby("species"):
            x = g.Dry_wheat_mass.values; r = dict(soil=SOIL[sc], species=sp, n=len(x), mean_g=x.mean(), sd_g=x.std(ddof=1))
            if sp not in CTRL:
                lo, hi = _boot_ci(x, ww); lo2, hi2 = _boot_ci(x, al)
                r.update(delta_vs_wxw_g=x.mean() - ww.mean(), ci_lo_wxw=lo, ci_hi_wxw=hi,
                         P_vs_wxw=stats.mannwhitneyu(x, ww, alternative="two-sided").pvalue,
                         delta_vs_alone_g=x.mean() - al.mean(), ci_lo_alone=lo2, ci_hi_alone=hi2,
                         P_vs_alone=stats.mannwhitneyu(x, al, alternative="two-sided").pvalue)
            rows.append(r)
    t = pd.DataFrame(rows); k = t.P_vs_wxw.notna()
    t.loc[k, "q_vs_wxw"] = multipletests(t.loc[k, "P_vs_wxw"], method="fdr_bh")[1]
    t.loc[k, "q_vs_alone"] = multipletests(t.loc[k, "P_vs_alone"], method="fdr_bh")[1]
    t["is_control"] = ~k; t = t.sort_values(["soil", "is_control", "mean_g"], ascending=[True, False, False]).drop(columns="is_control")
    if out_xlsx:
        readme = pd.DataFrame({"Table S3": [
            "Wheat aboveground biomass under each neighbour, per soil (sequenced pots).",
            "Δ, difference in mean biomass (g) from the wheat × wheat control (W×W) or from wheat alone, with bootstrap 95% CI (4,999 resamples).",
            "P, two-sided Mann–Whitney U; q, Benjamini–Hochberg across the 56 neighbour × soil contrasts, separately for each control."]})
        with pd.ExcelWriter(out_xlsx) as w:
            readme.to_excel(w, sheet_name="README", index=False); t.round(4).to_excel(w, sheet_name="Table S3", index=False)
    return t


def scale_check(m):
    """Why log scale (Methods 2.7, Methods S3): Spearman correlation between group mean and SD of wheat biomass
    across the treatment x soil groups, on the raw and log scales."""
    rows = []
    for lab, v in (("raw", m.Dry_wheat_mass), ("log", np.log(m.Dry_wheat_mass))):
        g = v.groupby([m.soil, m.species]).agg(["mean", "std"])
        r = stats.spearmanr(g["mean"], g["std"])
        rows.append(dict(scale=lab, groups=len(g), spearman_rho_mean_sd=r[0], P=r[1]))
    return pd.DataFrame(rows)
