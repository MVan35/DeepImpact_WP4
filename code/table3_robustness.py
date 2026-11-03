# -*- coding: utf-8 -*-
"""Robustness of Table 2 (Notes S6): the Distance term under other dissimilarities, without the Lingoes correction,
without Juncus bufonius and as a monocot/eudicot contrast (traits permuted among species), and all traits on
species-centroid communities. db-RDA F computed as in trait_models_v2.R.
Usage (from this folder):
  python table3_robustness.py ../data/matrice_phase1/df species_traits.csv \
         ../results/supporting/trait_models/check_phylogenetic_distance.csv ../results/supporting"""
import sys, os
import numpy as np, pandas as pd
from scipy.spatial.distance import pdist, squareform

DF, TRAITS, DIST, OUT = sys.argv[1:5]
B = 999
rng = np.random.default_rng(7)

def gower(D, lingoes=True):
    n = len(D); J = np.eye(n) - 1 / n; G = J @ (-0.5 * D ** 2) @ J
    if lingoes:
        ev = np.linalg.eigvalsh(G)
        if ev.min() < -1e-10 * abs(ev).max():
            D2 = np.sqrt(D ** 2 + 2 * (-ev.min())); np.fill_diagonal(D2, 0); G = J @ (-0.5 * D2 ** 2) @ J
    return G

def Fstat(x, G):
    X = np.column_stack([np.ones(len(G)), x]); Q, _ = np.linalg.qr(X)
    s = np.sum(Q * (G @ Q)); tot = np.trace(G); p = X.shape[1] - 1
    return (s / p) / ((tot - s) / (len(G) - p - 1))

def species_perm(sp, val, G):
    u = pd.Series(val, index=sp).groupby(level=0).first()
    x = pd.Series(sp).map(u).values.astype(float); Fo = Fstat(x, G)
    Fp = np.array([Fstat(pd.Series(sp).map(dict(zip(u.index, rng.permutation(u.values)))).values.astype(float), G) for _ in range(B)])
    return f"({Fo:.3f}, {(np.sum(Fp >= Fo) + 1) / (B + 1):.3f})"

ST = pd.read_csv(TRAITS).set_index("species"); cd = pd.read_csv(DIST).set_index("species")
mono = {s: int(cd.loc[s, "recomputed"] < 200) for s in cd.index}
TR = {"Distance": cd.recomputed, "Segetal": ST.Segetal, "Nitrogen": ST.EllenN, "Light": ST.EllenL,
      "SRL": ST.SRL, "RDMC": ST.RDMC, "RootShoot": ST.RootShoot}
rows, cen_rows = [], []
for kg, tag in (("bacteria", "1314"), ("fungi", "647")):
    otu = pd.read_csv(os.path.join(DF, f"otu_{kg}_table_R_{tag}ASVs.csv"), sep="\t", index_col=0).T.astype(float)
    md = pd.read_csv(os.path.join(DF, f"metadata_{kg}_{tag}ASVs.csv")).set_index("sample")
    for code, sl in (("resistant", "cereal"), ("Permissive", "mixed")):
        m = md.loc[md.index.intersection(otu.index)]
        m = m[(m.soil == code) & ~m.species.isin(["wheat alone", "wheat x wheat"])]
        X = otu.loc[m.index].values; sp = m.species.values
        Dbc = squareform(pdist(X, "braycurtis"))
        H = np.sqrt(X / X.sum(1, keepdims=True))
        dist = pd.Series(sp).map(cd.recomputed).values
        nj = ~pd.Series(sp).isin(["Juncus bufonius"]).values
        r = dict(kingdom=kg, soil=sl, n=len(sp))
        r["Distance BC+Lingoes"] = species_perm(sp, dist, gower(Dbc))
        r["Distance BC no correction"] = species_perm(sp, dist, gower(Dbc, False))
        r["Distance Hellinger"] = species_perm(sp, dist, gower(squareform(pdist(H, "euclidean")), False))
        r["Distance Jaccard (presence)"] = species_perm(sp, dist, gower(squareform(pdist(X > 0, "jaccard"))))
        r["Monocot vs eudicot (binary)"] = species_perm(sp, pd.Series(sp).map(mono).values, gower(Dbc))
        r["Distance, Juncus removed"] = species_perm(sp[nj], dist[nj], gower(Dbc[np.ix_(nj, nj)]))
        rows.append(r)
        Xr = otu.loc[m.index]; Xr = Xr.div(Xr.sum(1), axis=0); C = Xr.groupby(m.species.values).mean()
        c = dict(kingdom=kg, soil=sl)
        for t, s in TR.items():
            v = s.reindex(C.index).astype(float); ok = v.notna().values
            G = gower(squareform(pdist(C.values[ok], "braycurtis"))); x = v.values[ok]
            Fo = Fstat(x, G); Fp = np.array([Fstat(rng.permutation(x), G) for _ in range(B)])
            c[t] = f"{Fo:.2f} ({(np.sum(Fp >= Fo) + 1) / (B + 1):.3f})"
        cen_rows.append(c)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "table3_robustness.csv"), index=False)
pd.DataFrame(cen_rows).to_csv(os.path.join(OUT, "table3_species_centroid_check.csv"), index=False)
print(pd.DataFrame(rows).to_string(index=False)); print(pd.DataFrame(cen_rows).to_string(index=False))
