# -*- coding: utf-8 -*-
"""Overlap of the FDR differential-abundance core (Fig. 2) with the TITAN2 indicators of wheat biomass (z+ and z-).
Universe = the eligible set of each kingdom x soil (present in >= 10% of pots with biomass, mean relative abundance >= 1e-4), as in
the ranking test. Expected overlap = n_core * n_set / N; P = hypergeometric upper tail; pooled P = exact distribution of the summed
count over the four kingdom x soil cells (convolution), by exact convolution.
Usage: python core_vs_indicators.py <matrix_dir> <Table_S5_biomass_associated_taxa.csv> <Fig2_FDR_ASV_core_values.csv> <out.csv>"""
import sys
import numpy as np, pandas as pd
from scipy.stats import hypergeom
DF, S5, CORE, OUT = sys.argv[1:5]
s5, core = pd.read_csv(S5), pd.read_csv(CORE)
code = {"Cereal-dominated (AF193)": "resistant", "Mixed-system (AF190)": "Permissive"}
rows, conv, obs = [], {k: np.array([1.0]) for k in ("z+", "z-")}, {k: 0 for k in ("z+", "z-")}
for kg, n in (("bacteria", 1314), ("fungi", 647)):
    X0 = pd.read_csv(f"{DF}/otu_{kg}_table_R_{n}ASVs.csv", sep="\t", index_col=0).T.astype(float)
    md = pd.read_csv(f"{DF}/metadata_{kg}_{n}ASVs.csv").set_index("sample")
    for so, sc in code.items():
        m = md.loc[md.index.intersection(X0.index)]; m = m[(m.soil == sc) & m.Dry_wheat_mass.notna()]
        X = X0.loc[m.index]; X = X.div(X.sum(1), axis=0)
        elig = set(X.columns[((X > 0).mean() >= 0.10) & (X.mean() >= 1e-4)]); N = len(elig)
        c = set(core[(core.kingdom == kg) & (core.soil == so)].ASV) & elig
        g = s5[(s5.kingdom == kg) & (s5.soil == so)]
        sets = {"z+": set(g[g.TITAN2_direction == "z+"].ASV) & elig,
                "z-": set(g[g.TITAN2_direction == "z-"].ASV) & elig}
        row = dict(kingdom=kg, soil=so, eligible=N, core_in_eligible=len(c))
        for k, A in sets.items():
            o = len(c & A); obs[k] += o
            conv[k] = np.convolve(conv[k], hypergeom(N, len(A), len(c)).pmf(np.arange(len(c) + 1)))
            row.update({f"{k}_n": len(A), f"{k}_shared": o, f"{k}_expected": round(len(c) * len(A) / N, 2),
                        f"{k}_P": round(hypergeom(N, len(A), len(c)).sf(o - 1), 3) if o else 1.0,
                        f"{k}_ASVs": ", ".join(sorted(x.replace("Cluster_", "ASV_") for x in c & A))})
        rows.append(row)
T = pd.DataFrame(rows)
T.loc[len(T)] = dict(kingdom="all", soil="pooled", **{f"{k}_shared": obs[k] for k in obs},
                     **{f"{k}_expected": round(T[f"{k}_expected"].sum(), 2) for k in obs}, **{f"{k}_P": round(conv[k][obs[k]:].sum(), 3) for k in obs})
T.to_csv(OUT, index=False)
print(T.to_string(index=False))
