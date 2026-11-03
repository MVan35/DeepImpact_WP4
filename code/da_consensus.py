# -*- coding: utf-8 -*-
"""
High-confidence neighbour-enriched taxa (Fig. 2, Table S3, Notes S3, Results 3.1.2).

A taxon is retained for a (soil, neighbour) pair when at least two of DESeq2, edgeR and
metagenomeSeq call it ENRICHED under the neighbour relative to wheat x wheat:
    DESeq2         p < 0.05      and log2FC(neighbour/control) >  2
    edgeR          PValue < 0.01 and logFC(neighbour/control)  >  2
    metagenomeSeq  adj. p < 0.05 and coefficient(neighbour)    >  0
These are the thresholds of the original analysis (Analyses_bacteria_v1.2.ipynb, cell 77).

What this module changes relative to that analysis, and why
-----------------------------------------------------------
1. Direction is taken from the data, not assumed. The three packages export fold
   changes with different signs (DESeq2: log2(control/neighbour) when the control is
   named first; edgeR: log2(neighbour/control); metagenomeSeq: neighbour coefficient).
   The original code applied the DESeq2 rule to edgeR, so edgeR contributed DEPLETED
   taxa. Here every file is correlated with log2(neighbour/control) computed from the
   abundance table and flipped if the median correlation is negative; a file whose
   |median r| < 0.5 is refused (wrong ASV numbering or wrong pairing).
2. Neighbour labels are checked. In the DESeq2/metagenomeSeq exports two Poaceae carry
   wrong names ("Festuca myuros" is Lolium arundinaceum and "Viola arvensis" is Vulpia
   myuros; both confirmed by row counts and raw-abundance correlation). Every comparison
   must correlate best with the neighbour it is mapped to, otherwise the run stops.
3. metagenomeSeq rows are keyed "species.Cluster_n"; they are split before matching.
4. The edgeR contrasts are the re-export from edger_export.R (identical to the old files
   where those were complete; the mixed-system bacterial file previously lacked the
   wheat x wheat contrast).
"""
import glob, os, re
import numpy as np
import pandas as pd

CONTROL = "wheat x wheat"
SYN = {"Centaurea cyanus": "Cyanus segetum", "Bromus sterilis": "Anisantha sterilis",
       "Festuca arundinacea": "Lolium arundinaceum",
       "Festuca myuros": "Lolium arundinaceum",   # mislabelled in the DESeq2 / metagenomeSeq exports
       "Viola arvensis": "Vulpia myuros"}        # idem
SOILCODE = {"Resistant": "resistant", "Permissive": "Permissive"}
SOILLABEL = {"Resistant": "cereal-dominated", "Permissive": "mixed-system"}
THR = {"DESeq2": ("p", 0.05, 2.0), "edgeR": ("p", 0.01, 2.0), "metagenomeSeq": ("p", 0.05, 0.0)}
FILES = {"DESeq2": ("DESeq2_species", "DESeq_{kg}_{soil}_species_vs_wheat_x_wheat_pval.csv"),
         "edgeR": ("edgeR_species", "edgeR_{kg}_{soil}_vs_wheat_x_wheat.csv"),
         "metagenomeSeq": ("metagenomeSeq_species", "MetaSeq_{kg}_{soil}_CSS_wheat_x_wheat.csv")}


def norm_name(s):
    s = str(s).replace(".", " ").strip()
    for a, b in (("bursa pastoris", "bursa-pastoris"), ("speculum veneris", "speculum-veneris"),
                 ("foenum graecum", "foenum-graecum")):
        s = s.replace(a, b)
    return SYN.get(s, s)


def load_tables(matrix_dir):
    T = {}
    for kg, tag in (("bacteria", "1314"), ("fungi", "647")):
        otu = pd.read_csv(os.path.join(matrix_dir, f"otu_{kg}_table_R_{tag}ASVs.csv"), sep="\t", index_col=0).T.astype(float)
        md = pd.read_csv(os.path.join(matrix_dir, f"metadata_{kg}_{tag}ASVs.csv")).set_index("sample")
        md["species"] = md["species"].map(norm_name)
        otu = otu.loc[otu.index.intersection(md.index)]
        tax = pd.read_csv(os.path.join(matrix_dir, f"tax_{kg}_R_{tag}ASVs.csv"), sep="\t", index_col=0)
        T[kg] = dict(X=otu.div(otu.sum(1), axis=0), md=md, tax=tax)
    return T


def raw_lfc(T, kg, soil, weed):
    X, md = T[kg]["X"], T[kg]["md"]
    s = md[md.soil == SOILCODE[soil]]
    ia = s.index[s.species == weed].intersection(X.index)
    ic = s.index[s.species == CONTROL].intersection(X.index)
    if len(ia) == 0 or len(ic) == 0:
        return None
    return np.log2((X.loc[ia].mean() * 7000 + 0.5) / (X.loc[ic].mean() * 7000 + 0.5))


def read_export(path, method):
    """Long table: weed, Cluster, lfc (as exported), p."""
    d = pd.read_csv(path)
    if method == "metagenomeSeq":
        key = d.iloc[:, 0].astype(str)
        d, key = d[key.str.contains(".Cluster_", regex=False)], key[key.str.contains(".Cluster_", regex=False)]  # drop "<sp>.NA" rows
        return pd.DataFrame({"weed": key.str.split(".Cluster_").str[0].map(norm_name),
                             "Cluster": "Cluster_" + key.str.split(".Cluster_").str[1],
                             "lfc": d.iloc[:, 1].astype(float), "p": d["adjPvalues"].astype(float)})
    a = d["comparison"].str.split(" vs ", expand=True)
    a0, a1 = a[0].map(norm_name), a[1].map(norm_name)
    if not ((a0 == CONTROL) | (a1 == CONTROL)).all():
        raise ValueError(f"{os.path.basename(path)}: rows without the {CONTROL} contrast")
    lfc = d["log2FoldChange" if method == "DESeq2" else "logFC"].astype(float)
    # express everything as if the control were named first, so one sign applies per file
    lfc = np.where(a0 == CONTROL, lfc, -lfc)
    return pd.DataFrame({"weed": np.where(a0 == CONTROL, a1, a0), "Cluster": d["taxa"].astype(str), "lfc": lfc,
                         "p": d["pvalue" if method == "DESeq2" else "PValue"].astype(float)})


def orient_and_check(T, kg, soil, df, fname, min_r=0.5, label_margin=0.1):
    """Sign that turns df.lfc into log2(neighbour/control); verifies every neighbour label."""
    weeds = sorted(set(df.weed.unique()) - {"wheat alone"}); raws = {w: raw_lfc(T, kg, soil, w) for w in T[kg]["md"].species.unique()
                                             if w != CONTROL and w != "wheat alone"}
    rs, audit = [], []
    for w in weeds:
        g = df[df.weed == w]
        if w not in raws or raws[w] is None:
            raise ValueError(f"{fname}: neighbour '{w}' not in the metadata")
        cor = {}
        for cand, r in raws.items():
            if r is None: continue
            t = r.reindex(g.Cluster.values).values; ok = ~np.isnan(t)
            if ok.sum() >= 10: cor[cand] = np.corrcoef(g.lfc.values[ok], t[ok])[0, 1]
        if w not in cor: continue
        best = max(cor, key=lambda c: abs(cor[c]))
        # A label is wrong when another neighbour fits clearly better. Small exports (a few
        # dozen rows, mostly taxa absent from the control) correlate with every neighbour,
        # so near-ties are not evidence of a mix-up; the known mix-up differs by ~0.3.
        ok_ = abs(cor[w]) >= abs(cor[best]) - label_margin
        audit.append(dict(file=fname, neighbour=w, n_rows=len(g), r=round(cor[w], 3), best_match=best,
                          r_best=round(cor[best], 3), label_ok=ok_))
        rs.append(cor[w])
    med = float(np.median(rs))
    if abs(med) < min_r:
        raise ValueError(f"{fname}: median r with raw log2FC = {med:.2f}; file does not match the abundance table")
    bad = [a for a in audit if not a["label_ok"]]
    if bad:
        raise ValueError(f"{fname}: neighbour labels do not match the data: {bad[:3]}")
    return (1 if med > 0 else -1), med, audit


def enriched(T, roots, kg, soil):
    """All (method, neighbour, Cluster) calls enriched under the neighbour, plus a file audit."""
    calls, audit = [], []
    for m, (sub, pat) in FILES.items():
        f = next((os.path.join(r, sub, pat.format(kg=kg, soil=soil)) for r in roots
                  if os.path.exists(os.path.join(r, sub, pat.format(kg=kg, soil=soil)))), None)
        if f is None:
            raise FileNotFoundError(f"{m} export for {kg} {soil} not found under {roots}")
        df = read_export(f, m)
        df = df[df.weed != "wheat alone"]          # control-vs-control contrast, not a neighbour
        sign, med, a = orient_and_check(T, kg, soil, df, os.path.basename(f))
        audit += [dict(a, method=m, kingdom=kg, soil=soil, sign_applied=sign, median_r=round(med, 3), path=f) for a in a]
        df["lfc_n_over_c"] = sign * df.lfc
        _, pthr, lthr = THR[m]
        k = (df.p < pthr) & (df.lfc_n_over_c > lthr)
        calls.append(df.loc[k, ["weed", "Cluster"]].assign(Method=m))
    return pd.concat(calls, ignore_index=True).drop_duplicates(), pd.DataFrame(audit)


def consensus(matrix_dir, roots, min_methods=2):
    """Returns (long, audit): long has one row per method call for every retained
    (Cluster, Plant_Species, Soil) — the `filtered_df` layout the Fig. 2 builders expect."""
    T = load_tables(matrix_dir); out, aud = [], []
    for kg in ("bacteria", "fungi"):
        for soil in ("Resistant", "Permissive"):
            c, a = enriched(T, roots, kg, soil); aud.append(a)
            n = c.groupby(["Cluster", "weed"]).Method.nunique()
            keep = n[n >= min_methods].reset_index()[["Cluster", "weed"]]
            c = c.merge(keep, on=["Cluster", "weed"])
            out.append(c.rename(columns={"weed": "Plant_Species"}).assign(Kingdom=kg, Soil=soil))
    return pd.concat(out, ignore_index=True), pd.concat(aud, ignore_index=True), T


def table_s3(long, T):
    """One row per ASV, as in the published Table S3, with soil-specific neighbour lists."""
    rows = []
    for (kg, cl), g in long.groupby(["Kingdom", "Cluster"]):
        tax = T[kg]["tax"].reindex([cl]).iloc[0]
        soils = sorted({SOILLABEL[s] for s in g.Soil}, key=lambda s: s != "cereal-dominated")
        treat = "; ".join(f"{SOILLABEL[s]}: " + ", ".join(sorted(set(gg.Plant_Species)))
                          for s, gg in sorted(g.groupby("Soil"), key=lambda x: x[0] != "Resistant"))
        pair_methods = g.groupby(["Soil", "Plant_Species"]).Method.apply(lambda x: ", ".join(sorted(set(x))))
        rows.append(dict(Kingdom=kg, Cluster=cl, Phylum=tax.get("Phylum"), Genus=tax.get("Genus"),
                         Species=tax.get("Species"), Soil=", ".join(soils),
                         Methods=", ".join(sorted(set(g.Method))),
                         Method_combinations="; ".join(sorted(set(pair_methods))),
                         n_neighbours=g.groupby("Soil").Plant_Species.nunique().sum(),
                         Plant_Treatments=treat))
    t = pd.DataFrame(rows)
    t["_n"] = t.Cluster.str.extract(r"(\d+)").astype(int)
    return t.sort_values(["Kingdom", "_n"]).drop(columns="_n").reset_index(drop=True)


def summary(long):
    """Counts quoted in Results 3.1.2 / Notes S3."""
    pairs = long.drop_duplicates(["Kingdom", "Soil", "Plant_Species", "Cluster"])
    s = pairs.groupby(["Kingdom", "Soil"]).agg(ASVs=("Cluster", "nunique"), pairs=("Cluster", "size"),
                                               neighbours=("Plant_Species", "nunique")).reset_index()
    s["Soil"] = s.Soil.map(SOILLABEL)
    per = (long.groupby(["Kingdom", "Soil", "Plant_Species", "Cluster"]).Method
           .apply(lambda x: ", ".join(sorted(set(x)))).reset_index())
    combo = per.groupby(["Kingdom", "Method"]).size().unstack(fill_value=0)
    return s, combo


# --------------------------------------------------------------------- Fig. 2
# Genus -> putative function and colours, copied verbatim from
# Analyses_bacteria_v1.2.ipynb (cell 92) so the figure keeps its legend.
FUNCTION_BY_GENUS = {'Rhizobium': 'Nitrogen Fixer', 'Bradyrhizobium': 'Nitrogen Fixer', 'Pseudomonas': 'PGPR + Biocontrol', 'Bacillus': 'PGPR + Biocontrol', 'Azospirillum': 'PGPR', 'Streptomyces': 'Biocontrol agent', 'Trichoderma': 'Biocontrol + Nutrient cycling', 'Burkholderia': 'Stress tolerance / PGPR', 'Paenibacillus': 'PGPR + Biocontrol', 'Enterobacter': 'PGPR', 'Azotobacter': 'Nitrogen Fixer + PGPR', 'Arthrobacter': 'Stress-tolerant', 'Serratia': 'Biocontrol', 'Frankia': 'Nitrogen Fixer', 'Mesorhizobium': 'Nitrogen Fixer', 'Sinorhizobium': 'Nitrogen Fixer', 'Azoarcus': 'Nitrogen Fixer / PGPR', 'Pantoea': 'PGPR', 'Klebsiella': 'Nitrogen Fixer', 'Herbaspirillum': 'Nitrogen Fixer / PGPR', 'Devosia': 'Rhizosphere beneficial', 'Sphingomonas': 'Stress-tolerant', 'Flavobacterium': 'PGPR', 'Acinetobacter': 'PGPR', 'Microbacterium': 'PGPR', 'Variovorax': 'Plant-associated PGPR', 'Rhodococcus': 'Biocontrol', 'Lysobacter': 'Biocontrol', 'Niastella': 'Potential PGPR', 'Polaromonas': 'Mixed (some PGPR)', 'Hyphomicrobium': 'Nitrogen cycle (unclear PGPR)', 'Agrobacterium': 'Plant-associated PGPR', 'Afipia': 'Potential PGPR', 'Nitrobacter': 'Nitrogen cycle (unclear PGPR)', 'Massilia': 'PGPR', 'Mycobacterium': 'Decomposer'}
FUNCTION_PALETTE = {'Nitrogen Fixer': '#1f77b4', 'PGPR + Biocontrol': '#ff7f0e', 'PGPR': '#2ca02c', 'Stress tolerance / PGPR': '#d62728', 'Biocontrol agent': '#9467bd', 'Biocontrol + Nutrient cycling': '#8c564b', 'Decomposer / Biocontrol': '#e377c2', 'Decomposer': '#7f7f7f', 'Potential PGPR': '#bcbd22', 'Mixed (some PGPR)': '#17becf', 'Nitrogen Fixer + PGPR': '#aec7e8', 'Rhizosphere beneficial': '#c49c94', 'Stress-tolerant': '#f7b6d2', 'Nitrogen cycle (unclear PGPR)': '#dbdb8d', 'Plant-associated PGPR': '#9edae5', 'Biocontrol': '#ff9896', 'Nitrogen Fixer / PGPR': '#98df8a', 'PGPR + Nitrogen Fixer': '#c5b0d5'}
BACT_PHYLUM_COLOURS = {'Proteobacteria': '#332288', 'Bacteroidetes': '#999933', 'Firmicutes': '#CC6677', 'Actinobacteria': '#88CCEE', 'unidentified': '#AAAAAA', 'Chloroflexi': '#882255', 'Multi-affiliation': '#000000'}


def write_table_s3(long, T, path):
    """Table S3 workbook: README, one sheet per kingdom, summary, method combinations."""
    t = table_s3(long, T); s, combo = summary(long)
    with pd.ExcelWriter(path) as xw:
        pd.DataFrame({"Table S3": [
            "High-confidence neighbour-enriched ASVs: called enriched under the neighbour relative to wheat x wheat "
            "by at least two of DESeq2 (p < 0.05, log2FC > 2), edgeR (p < 0.01, logFC > 2) and metagenomeSeq "
            "(adj. p < 0.05, coefficient > 0) for the same neighbour and soil.",
            "Fold-change direction of every export was verified against log2(neighbour/control) from the abundance "
            "table (check_da_consensus_files.csv).",
            "Soil: soil(s) in which the ASV passed. Plant_Treatments: neighbours per soil. Method_combinations: the "
            "method sets supporting the (neighbour, soil) calls.",
            "Code: da_consensus.py (Paper_figures_v2.ipynb section 5); edgeR contrasts from edger_export.R "
            "(Paper_analyses_v2.Rmd section 13.1)."]}).to_excel(xw, sheet_name="README", index=False)
        for kg, sh in (("bacteria", "bacterial_clusters"), ("fungi", "fungi_clusters")):
            t[t.Kingdom == kg].drop(columns="Kingdom").to_excel(xw, sheet_name=sh, index=False)
        s.to_excel(xw, sheet_name="summary", index=False)
        combo.reset_index().to_excel(xw, sheet_name="method_combinations", index=False)
    return t


def fig2_inputs(long, T, min_combos_bacteria=5):
    """One row per (ASV, neighbour, soil) — the builders count rows, so this makes point
    size the number of enriched ASVs (the earlier frame had one row per method call).
    Bacterial genera enriched in fewer than `min_combos_bacteria` neighbour x soil
    combinations are left out of panel (a) for legibility; all are in Table S3."""
    u = long.drop_duplicates(["Kingdom", "Soil", "Plant_Species", "Cluster"])
    fam = {}
    for kg in ("bacteria", "fungi"):
        fam.update(dict(zip(T[kg]["md"].species, T[kg]["md"].family)))
    fb = u[u.Kingdom == "bacteria"][["Cluster", "Plant_Species", "Soil"]].reset_index(drop=True)
    ff = u[u.Kingdom == "fungi"][["Cluster", "Plant_Species", "Soil"]].reset_index(drop=True)
    gb = T["bacteria"]["tax"].Genus.to_dict()
    combos = fb.assign(G=fb.Cluster.map(gb)).drop_duplicates(["G", "Plant_Species", "Soil"]).groupby("G").size()
    keep = set(combos[combos >= min_combos_bacteria].index)
    return dict(bact=fb[fb.Cluster.map(gb).isin(keep)], bact_all=fb, fungi=ff, family=fam,
                genus_b=gb, phylum_b=T["bacteria"]["tax"].Phylum.to_dict(),
                genus_f=T["fungi"]["tax"].Genus.to_dict(), phylum_f=T["fungi"]["tax"].Phylum.to_dict(),
                n_genera_shown=len(keep), n_genera_all=int(combos.size))


def match_fig2_panels(fig_a, fig_b):
    """Give panel (b) the same row spacing (inches per y unit) and the same horizontal frame as (a).
    Both builders draw on an 18-in figure with bottom = 0.22 and top = 0.96, so panel (b), which has fewer
    genera, got more room per row. Here (b)'s figure height is set so that one row is the same physical
    height as in (a), the absolute bottom and top margins (x-axis labels, title room) are kept, and the axes
    get (a)'s left and right edges so the neighbour columns line up when the panels are stacked."""
    ax_a, ax_b = fig_a.axes[0], fig_b.axes[0]
    Wa, Ha = fig_a.get_size_inches(); pa = ax_a.get_position()
    per_unit = pa.height * Ha / np.ptp(ax_a.get_ylim())
    bottom_in, top_in = pa.y0 * Ha, (1 - pa.y1) * Ha
    Hb = bottom_in + top_in + per_unit * np.ptp(ax_b.get_ylim())
    fig_b.set_size_inches(Wa, Hb)
    fig_b.subplots_adjust(left=pa.x0, right=pa.x1, bottom=bottom_in / Hb, top=1 - top_in / Hb)
    fs_a = next((t.get_fontsize() for t in ax_a.texts if t.get_bbox_patch() is not None), None)
    separate_group_labels(fig_b, fontsize=fs_a)
    return Hb


def save_fig2(fig_a, fig_b, out_prefix, dpi=600, composite_dpi=300):
    """Save both panels with a common horizontal crop (so columns align) and a stacked composite PNG."""
    from matplotlib.transforms import Bbox
    from PIL import Image, ImageDraw, ImageFont
    Image.MAX_IMAGE_PIXELS = None
    boxes = []
    for f in (fig_a, fig_b):
        f.canvas.draw()
        boxes.append(f.get_tightbbox(f.canvas.get_renderer()))
    x0 = min(b.x0 for b in boxes) - 0.1; x1 = max(b.x1 for b in boxes) + 0.1
    paths = []
    for tag, f, b in (("a", fig_a, boxes[0]), ("b", fig_b, boxes[1])):
        bb = Bbox([[x0, b.y0 - 0.1], [x1, b.y1 + 0.1]])
        for ext in ("png", "pdf", "svg"):
            f.savefig(f"{out_prefix}{tag}.{ext}", dpi=dpi, bbox_inches=bb)
        paths.append(f"{out_prefix}{tag}.png")
    ims = [Image.open(p).convert("RGB") for p in paths]
    k = composite_dpi / dpi
    ims = [im.resize((round(im.width * k), round(im.height * k)), Image.LANCZOS) for im in ims]
    gap = round(0.25 * composite_dpi)
    C = Image.new("RGB", (max(i.width for i in ims), sum(i.height for i in ims) + gap), "white")
    C.paste(ims[0], (0, 0)); C.paste(ims[1], (0, ims[0].height + gap))
    d = ImageDraw.Draw(C)
    try:
        fnt = ImageFont.truetype("Arial Bold.ttf", round(0.35 * composite_dpi))
    except Exception:
        try:
            fnt = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", round(0.35 * composite_dpi))
        except Exception:
            fnt = ImageFont.load_default()
    d.text((10, 5), "a", fill="black", font=fnt); d.text((10, ims[0].height + gap + 5), "b", fill="black", font=fnt)
    C.save(f"{out_prefix}.png", dpi=(composite_dpi, composite_dpi))
    return f"{out_prefix}.png", C.size


def separate_group_labels(fig, fontsize=None, min_gap_pt=None, max_iter=50):
    """Phylum labels (boxed texts at group mid-points) of adjacent one-row groups collide once rows are
    spaced as in panel (a). Push neighbouring labels apart symmetrically until their centres are at least
    `min_gap_pt` apart (default 2.2 x font size), keeping them inside the y-limits. Optionally harmonise size."""
    ax = fig.axes[0]
    T = sorted([t for t in ax.texts if t.get_bbox_patch() is not None], key=lambda t: t.get_position()[1])
    if fontsize:
        for t in T: t.set_fontsize(fontsize)
    fs = T[0].get_fontsize() if T else 16
    gap = (min_gap_pt or 2.2 * fs) / 72.0                                   # inches
    per_unit = ax.get_position().height * fig.get_size_inches()[1] / np.ptp(ax.get_ylim())
    g = gap / per_unit                                                      # data units
    y = np.array([t.get_position()[1] for t in T], float)
    lo, hi = ax.get_ylim()[0] + 0.5 * g, ax.get_ylim()[1] - 0.5 * g
    for _ in range(max_iter):
        moved = False
        for i in range(len(y) - 1):
            d = y[i + 1] - y[i]
            if d < g:
                s = (g - d) / 2; y[i] -= s; y[i + 1] += s; moved = True
        y = np.clip(y, lo, hi)
        if not moved: break
    for t, yy in zip(T, y):
        t.set_position((t.get_position()[0], yy))
    return dict(zip([t.get_text() for t in T], np.round(y, 2)))


# ------------------------------------------------ overlap with the biomass indicators (Notes S4, Table S14)
SOIL_CODE = {"Resistant": "resistant", "Permissive": "Permissive"}          # metadata soil codes
S5_SOIL = {"Resistant": "Cereal-dominated (AF193)", "Permissive": "Mixed-system (AF190)"}


def eligible_set(T, kg, soil, prev_min=0.10, abund_min=1e-4):
    """Taxa eligible for TITAN2 / mEQO in one soil (same filter as candidate_pool.R and Table S8):
    present in >= 10% of pots with biomass and mean relative abundance >= 1e-4."""
    md = T[kg]["md"]; X = T[kg]["X"]
    s = md[(md.soil == SOIL_CODE[soil]) & md.Dry_wheat_mass.notna()].index.intersection(X.index)
    Xs = X.loc[s]
    return set(Xs.columns[((Xs > 0).mean() >= prev_min) & (Xs.mean() >= abund_min)])


def indicator_overlap(long, T, table_s5, biomass_up=None):
    """How many consensus-enriched ASVs are also TITAN2 indicators (z+, z-) or mEQO selections.
    Test: one-sided hypergeometric within each soil's eligible set (the population both indicator methods
    draw from), BH across all cells. `table_s5` is Table_S5_biomass_associated_taxa.csv.
    `biomass_up`, optional: set of (soil label, neighbour) whose wheat biomass exceeded wheat x wheat
    (Table S4, p < 0.05), used to report how many consensus calls come from such neighbours."""
    from scipy.stats import hypergeom
    from statsmodels.stats.multitest import multipletests
    LB = long.drop_duplicates(["Kingdom", "Soil", "Plant_Species", "Cluster"])
    rows, det = [], []
    for kg in ("bacteria", "fungi"):
        for soil in ("Resistant", "Permissive"):
            E = eligible_set(T, kg, soil)
            calls = LB[(LB.Kingdom == kg) & (LB.Soil == soil)]
            C = set(calls.Cluster); CE = C & E
            cell = table_s5[(table_s5.kingdom == kg) & (table_s5.soil == S5_SOIL[soil])]
            up = (np.mean([(SOILLABEL[soil].capitalize(), w) in biomass_up for w in calls.Plant_Species])
                  if biomass_up is not None else np.nan)
            for lab, S in (("TITAN2 z+", set(cell[cell.TITAN2_direction == "z+"].ASV)),
                           ("TITAN2 z-", set(cell[cell.TITAN2_direction == "z-"].ASV)),
                           ("mEQO selected", set(cell[cell.mEQO == "selected"].ASV))):
                SE = S & E; ov = CE & SE
                exp = len(CE) * len(SE) / len(E)
                p = hypergeom(len(E), len(SE), len(CE)).sf(len(ov) - 1) if SE else np.nan
                rows.append(dict(Kingdom=kg, Soil=SOILLABEL[soil], eligible=len(E), consensus_ASVs=len(C),
                                 consensus_in_eligible=len(CE), indicator_set=lab, set_size=len(S),
                                 set_in_eligible=len(SE), shared=len(C & S), shared_in_eligible=len(ov),
                                 expected=round(exp, 2), fold=round(len(ov) / exp, 2) if exp else np.nan,
                                 p_hypergeometric=p, share_calls_biomass_raising=round(up, 2)))
                for a in sorted(C & S):
                    g = calls[calls.Cluster == a]
                    det.append(dict(Kingdom=kg, Soil=SOILLABEL[soil], indicator_set=lab, ASV=a,
                                    Genus=T[kg]["tax"].Genus.get(a), in_eligible=a in E,
                                    also_mEQO=a in set(cell[cell.mEQO == "selected"].ASV),
                                    n_neighbours=g.Plant_Species.nunique(),
                                    neighbours=", ".join(sorted(set(g.Plant_Species)))))
    s = pd.DataFrame(rows)
    s["q_BH"] = multipletests(s.p_hypergeometric.fillna(1), method="fdr_bh")[1]
    return s, pd.DataFrame(det)


def biomass_raising_neighbours(T, alpha=0.05, control=CONTROL):
    """(soil label, neighbour) pairs whose wheat biomass exceeds wheat x wheat: two-sided Mann-Whitney U,
    uncorrected, as in Table S4 (which this reproduces), and a higher mean."""
    from scipy.stats import mannwhitneyu
    md = T["bacteria"]["md"].dropna(subset=["Dry_wheat_mass"])
    out = set()
    for soil, code in SOIL_CODE.items():
        m = md[md.soil == code]; c = m[m.species == control].Dry_wheat_mass.values
        for w, g in m.groupby("species"):
            if w in (control, "wheat alone"): continue
            x = g.Dry_wheat_mass.values
            if mannwhitneyu(x, c, alternative="two-sided").pvalue < alpha and x.mean() > c.mean():
                out.add((SOILLABEL[soil].capitalize(), w))
    return out
