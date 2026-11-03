# ============================================================================
#  Fungal genera dot-plot — revised, parameterised (single notebook cell)
#  Matches the bacteria figure: same y-axis spacing, same fixes.
#  Fixes: (1) figure-size & (2) y-axis-size knobs, (3) larger markers,
#         (4) no empty leading column, (5) phylum labels outside the box,
#         (6) corrected + de-capitalised species names on the x-axis.
#  Coloured by PHYLUM (fungal convention), not function.
#  Assumes these already exist in your notebook:
#     filtered_df, cluster_2_genus_d, cluster_2_phylum_d, plant_2_family
# ============================================================================
import pandas as pd
import matplotlib.pyplot as plt

# fungal phylum palette (edit colours here if you like)
FUNGI_PHYL_COL = {
    "Ascomycota":"#87CEEB","Basidiomycota":"#32CD32","Mortierellomycota":"#9370DB",
    "Glomeromycota":"#8B4513","Mucoromycota":"#228B22","Chytridiomycota":"#FF69B4",
    "Rozellomycota":"#BA55D3","Olpidiomycota":"#FFD700","unidentified":"#A9A9A9",
    "Multi-affiliation":"#000000",
}

# ---- species-name corrections (key = raw name in YOUR data, lower-cased) ----
SPECIES_RENAME = {
    "centaurea cyanus":    "Cyanus segetum",
    "bromus sterilis":     "Anisantha sterilis",
    "festuca arundinacea": "Lolium arundinaceum",   # only if present in your data
    "wheat alone":         "wheat alone",
    "wheat x wheat":       "wheat x wheat",
}
def _skey(raw):
    return str(raw).strip().lower().replace("_", " ")
def format_species_name(raw):
    """Correct known misnamings; else capitalise first letter only
       ('daucus carota' -> 'Daucus carota'; hyphen epithets stay lower)."""
    k = _skey(raw)
    if k in SPECIES_RENAME:                    return SPECIES_RENAME[k]
    if k in ("wheat alone", "wheat x wheat"):  return k
    return k.capitalize()

def create_combined_fungi_plot(
        multi_method_data, cluster_2_genus_d, cluster_2_phylum_d, plant_2_family,
        phyl_col_d=None,
        # -------------------- TUNABLE KNOBS ---------------------------------
        fig_width=32, fig_height=18,     # (1) figure size, inches
        y_fontsize=17,                   # (2) y-axis (genus) label font size
        x_fontsize=20,                   # x-axis (species) label font size
        point_scale=100, point_base=150, # (3) marker size = count*scale + base (orig 100 / 90)
        phylum_gap=1.2,                  # (5) phylum-label distance right of box
        species_spacing=3,
        legend_x=1.13,                   # legend distance right of box (raise if labels collide)
        out_prefix=None):                # if set: saves png(600dpi)/pdf/svg/csv
    if phyl_col_d is None:
        phyl_col_d = FUNGI_PHYL_COL

    # family resolver keyed on RAW species name (independent of display rename)
    fam_by_key = {_skey(k): v for k, v in plant_2_family.items()}

    gd = []
    for _, r in multi_method_data.iterrows():
        c = r['Cluster']
        gd.append({"raw":    r['Plant_Species'],
                   "Genus":  cluster_2_genus_d.get(c, "Unknown"),
                   "Phylum": cluster_2_phylum_d.get(c, "unidentified"),
                   "Soil":   r['Soil']})
    df = pd.DataFrame(gd)

    # (6) correct + de-capitalise species names; resolve family from raw key
    df["Species"] = df["raw"].apply(format_species_name)
    df["Family"]  = df["raw"].apply(lambda s: fam_by_key.get(_skey(s), "Unknown"))
    miss = df.loc[df["Family"] == "Unknown", "raw"].unique()
    if len(miss):
        print("WARNING: no family mapping for:", list(miss))

    df_plot = (df.groupby(["Family","Species","Genus","Phylum","Soil"])
                 .size().reset_index(name="Cluster_Count"))
    df_plot["Plant_Group"] = df_plot["Family"] + " - " + df_plot["Species"]
    ordered = df_plot.sort_values(["Family","Species"])["Plant_Group"].unique()

    # (4) first plant at position 0 -> no empty leading column
    plant_positions, pos = {}, 0
    for pl in ordered:
        plant_positions[pl] = pos; pos += species_spacing
    df_plot["x_pos"] = df_plot.apply(
        lambda r: plant_positions[r["Plant_Group"]] - 0.5 if r["Soil"] == "Resistant"
                  else plant_positions[r["Plant_Group"]] + 0.5, axis=1)

    # genera ordered within phylum (same y-spacing logic as bacteria fig)
    uniq_phyla = sorted(df_plot['Phylum'].unique())
    phyla_genera = {p: [] for p in uniq_phyla}
    for _, r in df_plot.drop_duplicates(['Genus','Phylum']).iterrows():
        if r['Genus'] not in phyla_genera[r['Phylum']]:
            phyla_genera[r['Phylum']].append(r['Genus'])
    for p in phyla_genera: phyla_genera[p].sort()
    gpp, yp, phylum_positions = {}, 0, {}
    for p in uniq_phyla:
        phylum_positions[p] = yp
        for g in phyla_genera[p]:
            gpp[(g, p)] = yp; yp += 1
    df_plot['y_pos'] = df_plot.apply(lambda r: gpp.get((r['Genus'], r['Phylum']), 0), axis=1)

    # (1) figure size
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    soil_edge = {"Resistant": "#404040", "Permissive": "#A0A0A0"}
    for soil, edge in soil_edge.items():
        sd = df_plot[df_plot["Soil"] == soil]
        ax.scatter(sd["x_pos"], sd["y_pos"],
                   s=sd["Cluster_Count"] * point_scale + point_base,   # (3) marker size
                   c=[phyl_col_d.get(p, "#7f7f7f") for p in sd["Phylum"]],
                   marker='o', edgecolor=edge, linewidth=2, alpha=0.8, label=soil)

    # y ticks
    yt, ytl = [], []
    for p in uniq_phyla:
        for g in phyla_genera[p]:
            yt.append(gpp[(g, p)]); ytl.append(g)
    ax.set_yticks(yt); ax.set_yticklabels(ytl, fontsize=y_fontsize, fontstyle='italic')
    ax.set_ylabel("Fungal Genus", fontsize=y_fontsize + 7)

    # x ticks — corrected names, italic
    xt  = [plant_positions[pl] for pl in ordered]
    xtl = [pl.split(" - ")[1] for pl in ordered]
    ax.set_xticks(xt); ax.set_xticklabels(xtl, rotation=90, fontsize=x_fontsize, fontstyle='italic')

    # phylum band separators
    for p in uniq_phyla[1:]:
        ax.axhline(phylum_positions[p] - 0.5, color='gray', lw=2, alpha=0.7)

    # (5) limits first, then phylum labels just OUTSIDE the box
    xmax = max(plant_positions.values())
    ax.set_xlim(-species_spacing * 0.6, xmax + species_spacing * 0.6)
    ax.set_ylim(-1, yp)
    x_right = ax.get_xlim()[1]
    for p in uniq_phyla:
        n = len(phyla_genera[p])
        if n == 0: continue
        centre = phylum_positions[p] + n / 2 - 0.5
        ax.text(x_right + phylum_gap, centre, p, fontsize=y_fontsize, fontweight='bold',
                va='center', ha='left', clip_on=False,
                bbox=dict(facecolor=phyl_col_d.get(p, "#7f7f7f"), alpha=0.2,
                          edgecolor=phyl_col_d.get(p, "#7f7f7f"), boxstyle="round,pad=0.5"))

    # family separators
    cur = None
    for pl in ordered:
        fam = pl.split(' - ')[0]
        if fam != cur:
            ax.axvline(plant_positions[pl] - 1.5, color='gray', lw=2, alpha=0.6); cur = fam

    # legend (push further right than bacteria fig so it clears phylum labels)
    soil_leg = [
        plt.Line2D([0],[0], marker='o', color='w', markerfacecolor='gray',
                   markeredgecolor='#404040', markeredgewidth=2, markersize=13, label='Cereal-dominated Soil'),
        plt.Line2D([0],[0], marker='o', color='w', markerfacecolor='gray',
                   markeredgecolor='#A0A0A0', markeredgewidth=2, markersize=13, label='Mixed-system Soil')]
    phyl_leg = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=13, label=p)
                for p, c in phyl_col_d.items() if p in uniq_phyla]
    ax.legend(handles=soil_leg + phyl_leg, bbox_to_anchor=(legend_x, 1), loc="upper left",
              title="Soil Type & Phylum", title_fontsize=16, fontsize=13, frameon=True)

    ax.grid(axis='y', linestyle='--', alpha=0.3)
    fig.subplots_adjust(right=0.80, bottom=0.22, top=0.96)

    if out_prefix:
        fig.savefig(out_prefix + ".png", dpi=600, bbox_inches='tight')
        fig.savefig(out_prefix + ".pdf", bbox_inches='tight')
        fig.savefig(out_prefix + ".svg", bbox_inches='tight')
        df_plot.to_csv(out_prefix + ".csv", index=False)
    return fig, df_plot

# ---------------------------------------------------------------------------
# NOTE: this file was originally a notebook CELL, and ended with a top-level
# call plus plt.show(). That call executed on `import`, so
#   from fungi_dotplot_cell import create_combined_fungi_plot
# failed with `NameError: name 'filtered_df' is not defined` before the caller
# had a chance to supply anything. The call has been removed; the function is
# invoked from Paper_figures.ipynb instead.
#
# Required arguments, in order:
#   multi_method_data, cluster_2_genus_d, cluster_2_phylum_d, plant_2_family
# Note the FOURTH is plant_2_family (a species -> family lookup), not a plain
# ordering vector.
# ---------------------------------------------------------------------------
