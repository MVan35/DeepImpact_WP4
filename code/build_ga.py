"""Graphical abstract v2 (rev B) — Vanhove et al., neighbouring weeds x wheat root microbiota.
Redesign answering S. Cordeau's review comments:
  (1) soil descriptors (texture, pH, organic matter, N) added to capsule 1
  (2) capsule 2 shows the experimental design explicitly, with "x 10 replicates"
  (3) capsule 3 lists the measured growth variables and spells out what
      "collective, not single-taxon" means
  (4) variance bar is direct-labelled: no unexplained colour squares
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Ellipse, Circle, FancyArrow

# ---------------------------------------------------------------- style
BROWN   = "#7A5324"   # soil legacy
GREEN   = "#2E5B23"   # weed identity
AMBER   = "#B8860B"   # wheat performance
BLUE    = "#3B6FB5"   # soil x weed interaction
GREY    = "#CBCBC6"   # unexplained
INK     = "#1A1A1A"
CARD    = "#F4F1EA"
CARD2   = "#FFFFFF"
TEXTURE = ["#3E2A11", "#9A7B4F", "#DCCBA8"]   # clay / silt / sand (single-hue ramp)

FS_BASE, FS_ANN, FS_TICK = 8.0, 7.0, 6.0
KX = 4.25 / 7.09          # x-units per y-unit for physically square shapes

SOIL = {
    "AF193": dict(name="Cereal-dominated soil", code="AF193",
                  hist="continuous cereal rotation",
                  clay=16.7, silt=58.9, sand=24.5, ph=5.9, om=27, n=1.5, cn=10.5),
    "AF190": dict(name="Mixed-system soil", code="AF190",
                  hist="diversified rotation, legume cover crops, livestock",
                  clay=18.9, silt=59.9, sand=21.3, ph=5.9, om=33, n=2.0, cn=9.7),
}
VAR = {  # sequential PERMANOVA, Table 1
    "Bacteria ($\\it{rpoB}$)": [("Soil legacy", 2.6, BROWN), ("Weed identity", 10.6, GREEN),
                                ("Soil x weed", 8.5, BLUE), ("Unexplained", 78.3, GREY)],
    "Fungi (ITS1)":            [("Soil legacy", 1.7, BROWN), ("Weed identity", 11.0, GREEN),
                                ("Soil x weed", 9.5, BLUE), ("Unexplained", 77.9, GREY)],
}


# ---------------------------------------------------------------- glyphs
def wheat(ax, x, y, h=6.0, col=AMBER, lw=0.9):
    """Small upright cereal tuft."""
    for dx, hh in ((-0.9, 0.82), (0.0, 1.0), (0.9, 0.80)):
        ax.plot([x, x + dx * KX * h * 0.30], [y, y + h * hh], color=col, lw=lw,
                solid_capstyle="round", zorder=3)
        ax.add_patch(Ellipse((x + dx * KX * h * 0.30, y + h * hh + h * 0.07),
                             width=h * 0.16 * KX, height=h * 0.20,
                             facecolor=col, edgecolor="none", zorder=3))


def weed(ax, x, y, h=4.6, col=GREEN, lw=0.9):
    """Small broadleaf neighbour."""
    ax.plot([x, x], [y, y + h], color=col, lw=lw, solid_capstyle="round", zorder=3)
    for s in (-1, 1):
        ax.add_patch(Ellipse((x + s * h * 0.26 * KX, y + h * 0.72),
                             width=h * 0.42 * KX, height=h * 0.30,
                             angle=0, facecolor=col, edgecolor="none", zorder=3))


def pot_top(ax, cx, cy, r=8.5, nb=4, lw=0.9, label_neighbours=True):
    """Top view of one experimental unit: central wheat + nb neighbours."""
    ax.add_patch(Ellipse((cx, cy), width=2 * r * KX, height=2 * r,
                         facecolor=CARD2, edgecolor="#8A8A82", lw=lw, zorder=2))
    ax.add_patch(Ellipse((cx, cy), width=1.62 * r * KX, height=1.62 * r,
                         facecolor="#EDE7DB", edgecolor="#B9B4A8", lw=0.6, zorder=2))
    wheat(ax, cx, cy - 2.6, h=6.4)
    for dx, dy in ((-1, 0), (1, 0), (0, 1), (0, -1)):
        weed(ax, cx + dx * r * 0.60 * KX, cy + dy * r * 0.52 - 2.0, h=4.4)


def bullet(ax, x, y, text, size=FS_TICK, color=INK, weight="normal", dy=0):
    ax.text(x, y + dy, text, ha="left", va="center", fontsize=size,
            color=color, weight=weight, zorder=4)


def capsule(ax, x0, x1, y0, y1, number, title, colour):
    """Header bar + body card. Returns (bx0, bx1, by0, by1) of the body area."""
    hh = 6.6
    ax.add_patch(FancyBboxPatch((x0, y1 - hh), x1 - x0, hh,
                                boxstyle="round,pad=0,rounding_size=1.0",
                                facecolor=colour, edgecolor="none", zorder=3))
    ax.text(x0 + 1.6, y1 - hh / 2, f"{number}   {title}", ha="left", va="center",
            fontsize=FS_BASE, color="white", weight="bold", zorder=4)
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - hh - y0 - 0.6,
                                boxstyle="round,pad=0,rounding_size=1.0",
                                facecolor=CARD, edgecolor="#DAD4C6", lw=0.7, zorder=1))
    return x0 + 1.6, x1 - 1.6, y0 + 1.0, y1 - hh - 1.6


def texture_bar(ax, x0, w, y, h, s, names=False):
    """Clay/silt/sand stacked strip for one soil."""
    xx = x0
    for frac, col, nm in zip((s["clay"], s["silt"], s["sand"]), TEXTURE, ("clay", "silt", "sand")):
        wseg = w * frac / 100
        ax.add_patch(Rectangle((xx, y), wseg, h, facecolor=col, edgecolor="white",
                               lw=0.5, zorder=3))
        lbl = f"{frac:.0f}%"
        ax.text(xx + wseg / 2, y + h / 2, lbl,
                ha="center", va="center", fontsize=FS_TICK,
                color="white" if col == TEXTURE[0] else INK, zorder=4)
        xx += wseg


def chevron(ax, x, y, w=2.6, h=3.4, col="#B4AE9F"):
    ax.add_patch(FancyArrow(x, y, w, 0, width=h * 0.42, head_width=h,
                            head_length=w * 0.55, length_includes_head=True,
                            facecolor=col, edgecolor="none", zorder=3))


# ---------------------------------------------------------------- figure
def key_column(ax, x0, y_mid, entries, dy=4.0, sw=2.1):
    """Frameless swatch-first key: binds each colour to its factor name."""
    y = y_mid + dy * (len(entries) - 1) / 2
    for nm, col in entries:
        ax.add_patch(Rectangle((x0, y - sw / 2 * (7.28 / 5.6) * 0.62), sw * (5.6 / 7.28) * 0.62 * 2.0,
                               sw, facecolor=col, edgecolor="white", lw=0.5, zorder=4))
        ax.text(x0 + 2.9, y + sw / 2 - 0.05, nm, ha="left", va="center",
                fontsize=FS_TICK, color=INK, zorder=4)
        y -= dy


def fits(fig, ax, cards):
    """Texts whose bbox escapes the capsule rect they sit in."""
    r = fig.canvas.get_renderer()
    bad = []
    for t in fig.findobj(mpl.text.Text):
        if not (t.get_text().strip() and t.get_visible()):
            continue
        bb = t.get_window_extent(r)
        x, y = t.get_position()
        for (cx0, cx1, cy0, cy1) in cards:
            if cx0 - 1 <= x <= cx1 + 1 and cy0 - 1 <= y <= cy1 + 1:
                d0 = ax.transData.transform((cx0, cy0)); d1 = ax.transData.transform((cx1, cy1))
                if bb.x0 < d0[0] - 1 or bb.x1 > d1[0] + 1 or bb.y0 < d0[1] - 1 or bb.y1 > d1[1] + 1:
                    bad.append((t.get_text()[:38], round(x, 1), round(y, 1)))
                break
    return bad


def build(path="graphical_abstract_v2.png"):
    mpl.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42,
                         "svg.fonttype": "none"})
    fig = plt.figure(figsize=(7.28, 5.6), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 100, 100, facecolor="white", edgecolor="none", zorder=0))

    ax.text(2.2, 98.2, "Neighbouring weeds filter the wheat root microbiota more strongly than the soil does",
            ha="left", va="top", fontsize=FS_BASE, color=INK, weight="bold", zorder=4)

    TOP, BOT = 94.4, 43.0
    xs = [(2.2, 33.6), (35.9, 64.3), (66.6, 97.8)]
    cards = []

    # ============================================================ capsule 1
    bx0, bx1, by0, by1 = capsule(ax, *xs[0], BOT, TOP, "1", "Soil legacy sets the pool", BROWN)
    cards.append((xs[0][0], xs[0][1], BOT, TOP))
    ax.text(bx0, by1 - 0.2, "Two field soils, contrasting management history",
            ha="left", va="top", fontsize=FS_TICK, color="#5A554A", style="italic", zorder=4)
    y = by1 - 3.4
    for si, (key, hist_lines) in enumerate((("AF193", ["continuous cereal rotation"]),
                                            ("AF190", ["diversified rotation, legume cover crops,",
                                                       "livestock integration"]))):
        s = SOIL[key]
        ax.text(bx0, y, s["name"], ha="left", va="center", fontsize=FS_ANN,
                color=INK, weight="bold", zorder=4)
        hy = y - 2.7
        for hl in hist_lines:
            ax.text(bx0, hy, hl, ha="left", va="center", fontsize=FS_TICK,
                    color="#5A554A", zorder=4)
            hy -= 2.3
        texture_bar(ax, bx0, bx1 - bx0 - 0.4, hy - 2.9, 2.9, s)
        hy -= 2.9
        if si == 0:   # bind texture colours to names once, directly under the first bar
            xw = 0.0
            for nm, col in zip(("clay", "silt", "sand"), TEXTURE):
                ax.add_patch(Rectangle((bx0 + xw, hy - 3.0), 1.4, 1.5, facecolor=col,
                                       edgecolor="#B9B4A8", lw=0.4, zorder=4))
                ax.text(bx0 + xw + 1.9, hy - 2.25, nm, ha="left", va="center",
                        fontsize=FS_TICK, color="#5A554A", zorder=4)
                xw += 6.6
            ax.text(bx0 + xw - 0.4, hy - 2.25, "(% mineral)", ha="left", va="center",
                    fontsize=FS_TICK, color="#7C776C", zorder=4)
            hy -= 3.2
        ax.text(bx0, hy - 3.0, f"silt loam  ·  pH {s['ph']:.1f}  ·  C:N {s['cn']:.1f}",
                ha="left", va="center", fontsize=FS_TICK, color=INK, zorder=4)
        ax.text(bx0, hy - 5.3,
                f"organic matter {s['om']} g kg$^{{-1}}$  ·  total N {s['n']:.1f} g kg$^{{-1}}$",
                ha="left", va="center", fontsize=FS_TICK, color=INK, zorder=4)
        y = hy - 7.6
    ax.text(bx0, by0 + 2.4, "Same texture class and pH — the legacy",
            ha="left", va="center", fontsize=FS_TICK, color=INK, weight="bold", zorder=4)
    ax.text(bx0, by0 + 0.2, "differs in organic matter and nitrogen",
            ha="left", va="center", fontsize=FS_TICK, color=INK, weight="bold", zorder=4)

    chevron(ax, xs[0][1] + 0.3, (TOP + BOT) / 2)

    # ============================================================ capsule 2
    bx0, bx1, by0, by1 = capsule(ax, *xs[1], BOT, TOP, "2", "Weed neighbour filters it", GREEN)
    cards.append((xs[1][0], xs[1][1], BOT, TOP))
    ax.text(bx0, by1 - 0.2, "Controlled plant-matrix experiment",
            ha="left", va="top", fontsize=FS_TICK, color="#5A554A", style="italic", zorder=4)

    cx, cy = bx0 + 7.8, by1 - 12.2
    for k in (2, 1):
        ax.add_patch(Ellipse((cx + 1.8 * k * KX, cy - 1.6 * k), width=2 * 7.4 * KX,
                             height=2 * 7.4, facecolor="#FBFAF7", edgecolor="#B9B4A8",
                             lw=0.6, ls=(0, (2.2, 1.8)), zorder=1))
    pot_top(ax, cx, cy, r=7.4)
    ax.annotate("x 10 replicates", xy=(cx + 4.2 * KX, cy - 8.4),
                xytext=(bx0 + 13.2, by1 - 22.2), fontsize=FS_ANN, color=INK,
                weight="bold", ha="left", va="center", zorder=5,
                arrowprops=dict(arrowstyle="-", lw=0.7, color="#8A8A82",
                                shrinkA=0, shrinkB=2))
    ax.text(bx0 + 15.4, by1 - 4.4, "1 wheat plant,\ncentre", ha="left", va="top",
            fontsize=FS_TICK, color=AMBER, weight="bold", zorder=4)
    ax.text(bx0 + 15.4, by1 - 11.0, "4 neighbours,\nequidistant", ha="left", va="top",
            fontsize=FS_TICK, color=GREEN, weight="bold", zorder=4)

    yy = by1 - 25.2
    for txt, w, c in (("28 weed species  +  2 controls", "bold", INK),
                      ("controls: wheat alone, wheat x wheat", "normal", "#5A554A"),
                      ("x 2 soils  x  10 replicates  =  600 pots", "bold", INK)):
        ax.text(bx0, yy, txt, ha="left", va="center", fontsize=FS_TICK,
                color=c, weight=w, zorder=4)
        yy -= 2.5
    yy -= 1.5
    ax.text(bx0, yy, "Unrelated weeds enrich the same taxa", ha="left", va="center",
            fontsize=FS_TICK + 0.6, color=GREEN, weight="bold", zorder=4)
    yy -= 2.6
    for txt in ("38 ASVs enriched vs wheat x wheat;",
                "29 of them under weeds of ≥ 3 plant families",
                "7 also vs wheat alone: not plant size"):
        ax.text(bx0, yy, txt, ha="left", va="center", fontsize=FS_TICK,
                color="#5A554A", zorder=4)
        yy -= 2.5

    chevron(ax, xs[1][1] + 0.3, (TOP + BOT) / 2)

    # ============================================================ capsule 3
    bx0, bx1, by0, by1 = capsule(ax, *xs[2], BOT, TOP, "3", "Community and growth", AMBER)
    cards.append((xs[2][0], xs[2][1], BOT, TOP))
    ax.text(bx0, by1 - 0.2, "Measured on every wheat plant", ha="left", va="top",
            fontsize=FS_TICK, color="#5A554A", style="italic", zorder=4)
    yy = by1 - 4.2
    for txt in ("aboveground dry biomass (g)",
                "species means (28 weeds)",
                "neighbour biomass $B_{tot}$ (g)"):
        ax.text(bx0 + 0.6, yy, "•", ha="left", va="center", fontsize=FS_ANN, color=AMBER, zorder=4)
        ax.text(bx0 + 2.2, yy, txt, ha="left", va="center", fontsize=FS_TICK, color=INK, zorder=4)
        yy -= 2.8
    wheat(ax, bx1 - 3.2, by1 - 11.4, h=8.6, lw=1.0)

    yy -= 2.0
    ax.text(bx0, yy, "No single taxon tracks growth", ha="left", va="center",
            fontsize=FS_ANN, color=INK, weight="bold", zorder=4)
    yy -= 2.8
    ax.text(bx0, yy, "none after correction, at species level", ha="left", va="center",
            fontsize=FS_TICK, color="#5A554A", zorder=4)
    yy -= 4.2
    ax.text(bx0, yy, "Bacterial composition covaries with it", ha="left", va="center",
            fontsize=FS_ANN, color=INK, weight="bold", zorder=4)
    yy -= 2.8
    ax.text(bx0, yy, "R² = 5.3 % (adj. 1.7 %) across neighbours in", ha="left", va="center",
            fontsize=FS_TICK, color="#5A554A", zorder=4)
    yy -= 2.5
    ax.text(bx0, yy, "cereal-dominated soil, weak in mixed-system;", ha="left", va="center",
            fontsize=FS_TICK, color="#5A554A", zorder=4)
    yy -= 2.5
    ax.text(bx0, yy, "unrelated for fungi", ha="left", va="center", fontsize=FS_TICK,
            color="#5A554A", zorder=4)
    yy -= 4.2
    ax.text(bx0, yy, "Weed traits do not predict it", ha="left", va="center",
            fontsize=FS_ANN, color=INK, weight="bold", zorder=4)
    yy -= 2.8
    ax.text(bx0, yy, "7 traits = random profiles; only", ha="left", va="center",
            fontsize=FS_TICK, color="#5A554A", zorder=4)
    yy -= 2.5
    ax.text(bx0, yy, "monocot vs eudicot structures composition", ha="left", va="center",
            fontsize=FS_TICK, color="#5A554A", zorder=4)

    # ============================================== variance panel
    ax.text(2.2, 37.0, "What explains wheat root microbiota composition",
            ha="left", va="center", fontsize=FS_BASE, color=INK, weight="bold", zorder=4)
    ax.text(2.2, 33.6, "sequential PERMANOVA on Bray–Curtis dissimilarity; each bar = 100 % of the variance",
            ha="left", va="center", fontsize=FS_TICK, color="#5A554A", zorder=4)

    BX0, BXW, hbar = 27.0, 70.0, 5.0
    rows = list(VAR.items())
    ybar = [24.0, 13.6]
    key_column(ax, 3.0, (ybar[0] + hbar / 2 + ybar[1] + hbar / 2) / 2,
               [(nm, col) for nm, _, col in VAR[rows[0][0]]], dy=4.3)

    for (row_name, segs), yb in zip(rows, ybar):
        ax.text(BX0, yb + hbar + 1.5, row_name, ha="left", va="center",
                fontsize=FS_TICK, color=INK, weight="bold", zorder=4)
        xx = BX0
        for nm, pct, col in segs:
            w = BXW * pct / 100
            ax.add_patch(Rectangle((xx, yb), w, hbar, facecolor=col,
                                   edgecolor="white", lw=0.7, zorder=3))
            if w > 5.0:
                ax.text(xx + w / 2, yb + hbar / 2, f"{pct:.1f} %", ha="center", va="center",
                        fontsize=FS_TICK, color="white" if col != GREY else "#4A473F",
                        weight="bold", zorder=4)
            else:   # soil-legacy segment is too narrow to hold its own value
                ax.annotate(f"{pct:.1f} %", xy=(xx + w * 0.5, yb + hbar / 2),
                            xytext=(BX0 - 2.0, yb + hbar / 2), fontsize=FS_TICK, color=col,
                            weight="bold", ha="right", va="center", zorder=5,
                            arrowprops=dict(arrowstyle="-", lw=0.6, color=col,
                                            shrinkA=1, shrinkB=0))
            xx += w

    ax.text(2.6, 5.8,
            "Weed identity outweighs soil legacy 4–7x, and the soil x weed interaction is nearly as large:",
            ha="left", va="center", fontsize=FS_ANN, color=INK, zorder=4)
    ax.text(2.6, 3.3, "each weed's effect on the wheat root microbiota depends on the soil it grows in.",
            ha="left", va="center", fontsize=FS_ANN, color=INK, zorder=4)
    ax.text(2.6, 1.0, "n = 455 (bacteria, rpoB) / 449 (fungi, ITS1) wheat root endospheres; 30 neighbour treatments x 2 soils x 10 replicates.",
            ha="left", va="center", fontsize=FS_TICK, color="#7C776C", zorder=4)

    fig.savefig(path, dpi=300)
    fig.savefig(path.replace(".png", ".pdf"))
    return fig, ax, cards


if __name__ == "__main__":
    fig, ax, cards = build()
