# =============================================================
#  Bacteria genera dot-plot  —  revised, parameterised version
#  Fixes: (1) easy figure-size & y-axis-size knobs, (2) larger markers,
#         (3) no leading empty column, (4) phylum labels outside the box,
#         (5) corrected & de-capitalised species names on the x-axis.
# =============================================================
import pandas as pd
import matplotlib.pyplot as plt

# ---- species name corrections -------------------------------------------
# key = the raw name as it appears in YOUR data (lower-cased, spaces);
# value = the correct display name. Add/adjust entries as needed.
SPECIES_RENAME = {
    "centaurea cyanus":       "Cyanus segetum",
    "bromus sterilis":        "Anisantha sterilis",
    "festuca arundinacea":    "Lolium arundinaceum",   # if present in your data
    "wheat alone":            "wheat alone",
    "wheat x wheat":          "wheat x wheat",
}
def _skey(raw):
    return str(raw).strip().lower().replace("_", " ")
def format_species_name(raw):
    """Correct known misnamings; otherwise capitalise first letter only
    ('daucus carota' -> 'Daucus carota'; hyphen epithets stay lower)."""
    k = _skey(raw)
    if k in SPECIES_RENAME:
        return SPECIES_RENAME[k]
    if k in ("wheat alone", "wheat x wheat"):
        return k
    return k.capitalize()

def create_combined_bacteria_plot(
        filtered_df, cluster_2_genus_d, cluster_2_phylum_d, plant_2_family,
        beneficial_genera_d, function_palette, phyl_col_d,
        default_genus_color="#A9A9A9",
        # -------- TUNABLE KNOBS ---------------------------------------
        fig_width=32, fig_height=18,   # <<< figure size (inches)
        y_fontsize=17,                 # <<< y-axis (genus) label font size
        x_fontsize=22,                 # x-axis (species) label font size
        point_scale=60, point_base=200,# marker size = count*scale + base  (was 30 / 100)
        phylum_gap=1.2,                # distance of phylum labels to the RIGHT of the box
        species_spacing=3,
        out_prefix=None):              # if set, saves png/pdf/svg with this path prefix
    # ---- build genus-level long table ----
    gd = []
    for _, r in filtered_df.iterrows():
        c = r['Cluster']
        g  = cluster_2_genus_d.get(c, "Unknown")
        ph = cluster_2_phylum_d.get(c, "unidentified")
        gd.append({"Plant Species": r['Plant_Species'], "Genus": g, "Phylum": ph,
                   "Function": beneficial_genera_d.get(g, "Unknown"), "Soil": r['Soil']})
    df = pd.DataFrame(gd)

    # ---- FIX (5): correct + de-capitalise species names ----
    df["Species"] = df["Plant Species"].apply(format_species_name)
    df["Family"]  = df["Species"].map(lambda s: plant_2_family.get(s, "Unknown"))

    df_plot = (df.groupby(["Family", "Species", "Genus", "Phylum", "Function", "Soil"])
                 .size().reset_index(name="Cluster_Count"))
    df_plot["Plant_Group"] = df_plot["Family"] + " - " + df_plot["Species"]
    ordered = df_plot.sort_values(["Family", "Species"])["Plant_Group"].unique()

    # ---- FIX (3): first plant starts at position 0 (no empty leading column) ----
    plant_positions = {}
    pos = 0
    for pl in ordered:
        plant_positions[pl] = pos
        pos += species_spacing
    df_plot["x_pos"] = df_plot.apply(
        lambda r: plant_positions[r["Plant_Group"]] - 0.5 if r["Soil"] == "Resistant"
                  else plant_positions[r["Plant_Group"]] + 0.5, axis=1)

    # ---- order genera within phylum ----
    uniq_phyla = sorted(df_plot['Phylum'].unique())
    phyla_genera = {p: [] for p in uniq_phyla}
    for _, r in df_plot.drop_duplicates(['Genus', 'Phylum']).iterrows():
        if r['Genus'] not in phyla_genera[r['Phylum']]:
            phyla_genera[r['Phylum']].append(r['Genus'])
    for p in phyla_genera:
        phyla_genera[p].sort()
    gpp = {}; yp = 0; phylum_positions = {}
    for p in uniq_phyla:
        phylum_positions[p] = yp
        for g in phyla_genera[p]:
            gpp[(g, p)] = yp; yp += 1
    df_plot['y_pos'] = df_plot.apply(lambda r: gpp.get((r['Genus'], r['Phylum']), 0), axis=1)

    # ---- FIX (1): figure size via fig_width/fig_height; use fig,ax ----
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    soil_edge = {"Resistant": "#404040", "Permissive": "#A0A0A0"}
    for soil, edge in soil_edge.items():
        sd = df_plot[df_plot["Soil"] == soil]
        colors = [function_palette.get(f, default_genus_color) for f in sd["Function"]]
        ax.scatter(sd["x_pos"], sd["y_pos"],
                   s=sd["Cluster_Count"] * point_scale + point_base,   # <<< FIX (2): marker size
                   c=colors, marker='o', edgecolor=edge, linewidth=2, alpha=0.85, label=soil)

    # y ticks
    yt, ytl = [], []
    for p in uniq_phyla:
        for g in phyla_genera[p]:
            yt.append(gpp[(g, p)]); ytl.append(g)
    ax.set_yticks(yt); ax.set_yticklabels(ytl, fontsize=y_fontsize, fontstyle='italic')
    ax.set_ylabel("Bacterial Genus", fontsize=y_fontsize + 7)

    # x ticks — corrected names, one per plant, italic
    xt  = [plant_positions[pl] for pl in ordered]
    xtl = [pl.split(" - ")[1] for pl in ordered]
    ax.set_xticks(xt); ax.set_xticklabels(xtl, rotation=90, fontsize=x_fontsize, fontstyle='italic')

    # phylum band separators
    for p in uniq_phyla[1:]:
        ax.axhline(phylum_positions[p] - 0.5, color='gray', lw=2, alpha=0.7)

    # ---- FIX (4): set limits first, then place phylum labels just OUTSIDE the box ----
    xmax = max(plant_positions.values())
    ax.set_xlim(-species_spacing * 0.6, xmax + species_spacing * 0.6)  # tight left edge
    ax.set_ylim(-1, yp)
    x_right = ax.get_xlim()[1]
    for p in uniq_phyla:
        n = len(phyla_genera[p])
        if n == 0:
            continue
        centre = phylum_positions[p] + n / 2 - 0.5
        ax.text(x_right + phylum_gap, centre, p, fontsize=y_fontsize - 1, fontweight='bold',
                va='center', ha='left', clip_on=False,
                bbox=dict(facecolor=phyl_col_d.get(p, "#7f7f7f"), alpha=0.2,
                          edgecolor=phyl_col_d.get(p, "#7f7f7f"), boxstyle="round,pad=0.5"))

    # family separators
    cur = None
    for pl in ordered:
        fam = pl.split(' - ')[0]
        if fam != cur:
            ax.axvline(plant_positions[pl] - 1.5, color='gray', lw=2, alpha=0.6); cur = fam

    # legends (marker size in legend bumped to match)
    soil_leg = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                   markeredgecolor='#404040', markeredgewidth=2, markersize=13, label='Cereal-dominated Soil'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                   markeredgecolor='#A0A0A0', markeredgewidth=2, markersize=13, label='Mixed-system Soil')]
    func_leg = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=13, label=f)
                for f, c in function_palette.items() if f in df_plot["Function"].values]
    if "Unknown" in df_plot["Function"].values:
        func_leg.append(plt.Line2D([0], [0], marker='o', color='w',
                        markerfacecolor=default_genus_color, markersize=13, label="Unknown Function"))
    ax.legend(handles=soil_leg + func_leg, bbox_to_anchor=(1.10, 1), loc="upper left",
              title="Soil Type & Function", title_fontsize=16, fontsize=13, frameon=True)

    ax.grid(axis='y', linestyle='--', alpha=0.3)
    fig.subplots_adjust(right=0.82, bottom=0.22, top=0.96)

    if out_prefix:
        fig.savefig(out_prefix + ".png", dpi=600, bbox_inches='tight')
        fig.savefig(out_prefix + ".pdf", bbox_inches='tight')
        fig.savefig(out_prefix + ".svg", bbox_inches='tight')
        df_plot.to_csv(out_prefix + ".csv", index=False)
    return fig, df_plot
