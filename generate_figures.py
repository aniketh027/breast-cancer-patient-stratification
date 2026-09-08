"""
generate_figures.py
Generates Figures 1, 2, and 3 for the research paper:
  Fig 1 - Dual-pipeline workflow diagram
  Fig 2 - Side-by-side module size bar chart (colour-coded by pathway)
  Fig 3 - 4x4 Jaccard heatmap (STRING vs Disease PPI modules)

Data verified against BC_Pipeline_Clean_FINAL notebook outputs.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import os

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR  = os.path.join(PROJECT_DIR, "results", "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ==============================================================================
# FIGURE 1 - Dual-pipeline workflow diagram
# ==============================================================================

def draw_fig1():
    fig, ax = plt.subplots(figsize=(10, 12))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis("off")
    fig.patch.set_facecolor("#0d1117")
    ax.set_facecolor("#0d1117")

    # Colour palette
    COL_SHARED  = "#58a6ff"
    COL_STRING  = "#3fb950"
    COL_DISEASE = "#f0883e"
    COL_COMPARE = "#bc8cff"
    TEXT_COL    = "#f0f6fc"

    def add_box(ax, x, y, w, h, text, colour, fontsize=10, bold=False):
        box = FancyBboxPatch(
            (x - w/2, y - h/2), w, h,
            boxstyle="round,pad=0.15",
            facecolor=colour + "22",
            edgecolor=colour,
            linewidth=2.0,
            zorder=3
        )
        ax.add_patch(box)
        weight = "bold" if bold else "normal"
        ax.text(x, y, text, ha="center", va="center",
                fontsize=fontsize, color=TEXT_COL, fontweight=weight, zorder=4)

    def add_arrow(ax, x1, y1, x2, y2, colour="#8b949e"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                     arrowprops=dict(arrowstyle="-|>", color=colour,
                                     lw=2.0, mutation_scale=16),
                     zorder=2)

    # Title
    ax.text(5, 11.5, "End-to-End Dual-Pipeline Workflow\nfor Breast Cancer Gene Module Discovery",
            ha="center", va="center", fontsize=14, color=TEXT_COL, fontweight="bold")

    # Shared: Seed Ingestion (25 genes from notebook)
    add_box(ax, 5, 10.3, 3.8, 0.65,
            "Step 1: Seed Gene Ingestion\n(25 BRCA seed genes from GWAS)",
            COL_SHARED, fontsize=11, bold=True)

    # Fork arrows
    add_arrow(ax, 5, 9.95, 2.7, 9.25)
    add_arrow(ax, 5, 9.95, 7.3, 9.25)

    # Pipeline labels
    ax.text(2.7, 9.6, "Pipeline A\n(STRING)", ha="center", va="center",
            fontsize=12, color=COL_STRING, fontweight="bold")
    ax.text(7.3, 9.6, "Pipeline B\n(Disease PPI)", ha="center", va="center",
            fontsize=12, color=COL_DISEASE, fontweight="bold")

    # Pipeline A steps (left column, x=2.7)
    steps_a = [
        (8.5, "Step 2: STRING v12.0\nNetwork (16,155 nodes)"),
        (7.4, "Step 3: RWR Propagation (R)\n(r = 0.3, ~16K genes)"),
        (6.3, "Step 4: Proximity Testing\n(NetColoc, FDR < 0.05)"),
        (5.2, "Step 5: Subgraph\n(252 genes, 1,014 edges)"),
        (4.1, "Step 6: Leiden Clustering\n(4 modules, Q = 0.34)"),
        (3.0, "Step 7: Enrichment\n(KEGG + Reactome)"),
    ]
    for y, label in steps_a:
        add_box(ax, 2.7, y, 3.3, 0.7, label, COL_STRING, fontsize=9)
    for i in range(len(steps_a) - 1):
        add_arrow(ax, 2.7, steps_a[i][0] - 0.40, 2.7, steps_a[i+1][0] + 0.40, COL_STRING)

    # Pipeline B steps (right column, x=7.3)
    steps_b = [
        (8.5, "Step 8: Disease PPI\n(PSICQUIC+BioGRID, 9,504)"),
        (7.4, "Step 9: RWR Propagation (R)\n(r = 0.3, ~9.5K genes)"),
        (6.3, "Step 10: Proximity Testing\n(NetColoc, FDR < 0.05)"),
        (5.2, "Step 11: Subgraph\n(203 genes, 356 edges)"),
        (4.1, "Step 12: Leiden Clustering\n(4 modules, Q = 0.58)"),
        (3.0, "Step 14: Enrichment\n(KEGG + Reactome)"),
    ]
    for y, label in steps_b:
        add_box(ax, 7.3, y, 3.3, 0.7, label, COL_DISEASE, fontsize=9)
    for i in range(len(steps_b) - 1):
        add_arrow(ax, 7.3, steps_b[i][0] - 0.40, 7.3, steps_b[i+1][0] + 0.40, COL_DISEASE)

    # Comparative Analysis (bottom centre)
    add_arrow(ax, 2.7, 2.55, 5, 1.85, COL_COMPARE)
    add_arrow(ax, 7.3, 2.55, 5, 1.85, COL_COMPARE)
    add_box(ax, 5, 1.4, 4.5, 0.7,
            "Step 13: Comparative Analysis\n(ARI=0.37, NMI=0.55, Jaccard=0.38)",
            COL_COMPARE, fontsize=10, bold=True)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COL_SHARED+"44", edgecolor=COL_SHARED, label="Shared"),
        mpatches.Patch(facecolor=COL_STRING+"44", edgecolor=COL_STRING, label="Pipeline A"),
        mpatches.Patch(facecolor=COL_DISEASE+"44", edgecolor=COL_DISEASE, label="Pipeline B"),
        mpatches.Patch(facecolor=COL_COMPARE+"44", edgecolor=COL_COMPARE, label="Comparative"),
    ]
    ax.legend(handles=legend_elements, loc="lower center", fontsize=10,
              facecolor="#161b22", edgecolor="#30363d", labelcolor=TEXT_COL,
              framealpha=0.9, ncol=4, bbox_to_anchor=(0.5, -0.02))

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, "fig1_workflow.png")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[OK] Saved: {path}")


# ==============================================================================
# FIGURE 2 - Side-by-side module size bar chart
# ==============================================================================

def draw_fig2():
    # Module sizes from notebook (Cell 16 & 28 outputs)
    # STRING: Cluster 0=94, 1=84, 2=34, 3=25  (total 237 genes)
    # Disease: Cluster 0=55, 1=50, 2=31, 3=25  (total 161 genes)
    string_sizes  = [94, 84, 34, 25]
    disease_sizes = [55, 50, 31, 25]

    # Dominant pathway annotations from enrichment (Cell 18 & 32 outputs)
    # STRING Module 0: mTOR, MAPK, C-type lectin → PI3K/AKT/mTOR Signaling
    # STRING Module 1: Homologous recombination, Fanconi anemia → DNA Repair
    # STRING Module 2: TP53 Transcriptional Regulation (from Reactome)
    # STRING Module 3: Meiotic Recombination, Base Excision Repair
    # Disease Module 0: Homologous recomb, DNA repair, G2/M → DNA Damage Response
    # Disease Module 1: Homologous recombination → Homologous Recombination
    # Disease Module 2: TRAIL, mTOR, AMPK, Apoptosis → mTOR/Apoptosis
    # Disease Module 3: Endometrial cancer, PI3K-Akt, RTK → RTK/PI3K Signaling

    string_pathways = [
        "PI3K/AKT/mTOR\nSignaling",
        "DNA Repair /\nHom. Recombination",
        "TP53 Transcriptional\nRegulation",
        "Meiotic Recomb. /\nBase Excision Repair",
    ]
    disease_pathways = [
        "DNA Damage Response /\nCheckpoints",
        "Homologous\nRecombination",
        "mTOR / Apoptosis\nSignaling",
        "RTK / PI3K\nSignaling",
    ]

    pathway_colours = {
        "PI3K/AKT/mTOR\nSignaling":               "#58a6ff",
        "DNA Repair /\nHom. Recombination":         "#f97583",
        "TP53 Transcriptional\nRegulation":         "#d2a8ff",
        "Meiotic Recomb. /\nBase Excision Repair":  "#79c0ff",
        "DNA Damage Response /\nCheckpoints":       "#f97583",
        "Homologous\nRecombination":                "#ff7b72",
        "mTOR / Apoptosis\nSignaling":              "#58a6ff",
        "RTK / PI3K\nSignaling":                    "#7ee787",
    }

    fig, axes = plt.subplots(2, 1, figsize=(9, 10), sharex=True)
    fig.patch.set_facecolor("#0d1117")

    fig.suptitle("Fig. 2: Module Size Comparison\nSTRING vs Disease PPI Pipelines",
                 color="#f0f6fc", fontsize=14, fontweight="bold", y=0.98)

    for ax_idx, (ax, sizes, pathways, title, title_col) in enumerate(zip(
            axes,
            [string_sizes, disease_sizes],
            [string_pathways, disease_pathways],
            ["Pipeline A (STRING) -- 237 genes, 4 modules, Q = 0.34",
             "Pipeline B (Disease PPI) -- 161 genes, 4 modules, Q = 0.58"],
            ["#3fb950", "#f0883e"])):

        ax.set_facecolor("#161b22")
        colours = [pathway_colours[p] for p in pathways]
        x = np.arange(len(sizes))
        bars = ax.bar(x, sizes, width=0.5, color=colours, edgecolor="#30363d",
                      linewidth=1.5, zorder=3)

        # Value labels on bars
        for bar, val in zip(bars, sizes):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
                    str(val), ha="center", va="bottom",
                    color="#f0f6fc", fontsize=13, fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels([f"Module {i}" for i in range(len(sizes))],
                           color="#f0f6fc", fontsize=12)
        ax.set_title(title, color=title_col, fontsize=12, fontweight="bold", pad=10)
        ax.set_ylabel("Number of Genes", color="#f0f6fc", fontsize=11)
        ax.tick_params(axis="y", colors="#8b949e", labelsize=11)
        ax.tick_params(axis="x", colors="#8b949e")
        ax.spines["bottom"].set_color("#30363d")
        ax.spines["left"].set_color("#30363d")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="#21262d", linestyle="--", alpha=0.5, zorder=0)
        ax.set_ylim(0, 110)

        # Pathway annotation below each bar
        for i, pathway in enumerate(pathways):
            ax.text(i, -12, pathway, ha="center", va="top",
                    fontsize=9, color=colours[i], fontstyle="italic")

    fig.tight_layout(rect=[0, 0.12, 1, 0.95])

    # Unified legend
    unique_pathways = list(dict.fromkeys(list(string_pathways) + list(disease_pathways)))
    legend_patches = [mpatches.Patch(color=pathway_colours[p], label=p.replace("\n", " "))
                      for p in unique_pathways]
    fig.legend(handles=legend_patches, loc="lower center", ncol=2,
               fontsize=10, facecolor="#161b22", edgecolor="#30363d",
               labelcolor="#f0f6fc", framealpha=0.9,
               bbox_to_anchor=(0.5, 0.0))
    path = os.path.join(OUTPUT_DIR, "fig2_module_sizes.png")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[OK] Saved: {path}")


# ==============================================================================
# FIGURE 3 - 4x4 Jaccard heatmap
# ==============================================================================

def draw_fig3():
    # Jaccard matrix from notebook Step 13 output (module_pair_similarity.tsv)
    # Rows = STRING modules 0-3, Cols = Disease PPI modules 0-3
    # STRING 0 vs Disease: 0=0.0, 1=0.0, 2=0.0684, 3=0.0259
    # STRING 1 vs Disease: 0=0.0451, 1=0.0635, 2=0.0, 3=0.0
    # STRING 2 vs Disease: 0=0.0114, 1=0.0, 2=0.0, 3=0.0
    # STRING 3 vs Disease: 0=0.0526, 1=0.0274, 2=0.0, 3=0.0
    jaccard = np.array([
        [0.0000, 0.0000, 0.0684, 0.0259],
        [0.0451, 0.0635, 0.0000, 0.0000],
        [0.0114, 0.0000, 0.0000, 0.0000],
        [0.0526, 0.0274, 0.0000, 0.0000],
    ])

    fig, ax = plt.subplots(figsize=(9, 8))
    fig.patch.set_facecolor("#0d1117")
    ax.set_facecolor("#0d1117")

    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "dark_blue", ["#161b22", "#1f3a5f", "#2d6abf", "#58a6ff"], N=256
    )

    im = ax.imshow(jaccard, cmap=cmap, aspect="auto", vmin=0, vmax=0.08)

    # Annotate cells with Jaccard values
    for i in range(4):
        for j in range(4):
            val = jaccard[i, j]
            text_col = "#f0f6fc" if val > 0.03 else "#8b949e"
            ax.text(j, i, f"{val:.4f}", ha="center", va="center",
                    fontsize=14, color=text_col, fontweight="bold")

    ax.set_xticks(range(4))
    ax.set_yticks(range(4))
    ax.set_xticklabels([f"Disease M{i}" for i in range(4)], color="#f0f6fc", fontsize=12)
    ax.set_yticklabels([f"STRING M{i}" for i in range(4)], color="#f0f6fc", fontsize=12)
    ax.set_xlabel("Disease PPI Modules", color="#f0f6fc", fontsize=13, labelpad=12)
    ax.set_ylabel("STRING Modules", color="#f0f6fc", fontsize=13, labelpad=12)
    ax.set_title("Fig. 3: Pairwise Jaccard Similarity Between STRING and Disease PPI Modules",
                 color="#f0f6fc", fontsize=14, fontweight="bold", pad=18)

    ax.tick_params(axis="both", colors="#8b949e", labelsize=12)
    for spine in ax.spines.values():
        spine.set_color("#30363d")

    # Colourbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Jaccard Index", color="#f0f6fc", fontsize=12)
    cbar.ax.tick_params(colors="#8b949e", labelsize=11)
    cbar.outline.set_edgecolor("#30363d")

    # Highlight best matches (from notebook: Disease 0->STRING 3, Disease 1->STRING 1,
    # Disease 2->STRING 0, Disease 3->STRING 0)
    best_matches = [(0, 2), (1, 1), (2, 0), (3, 0)]
    for si, di in best_matches:
        if jaccard[si, di] > 0:
            rect = plt.Rectangle((di - 0.5, si - 0.5), 1, 1,
                                 fill=False, edgecolor="#f0883e",
                                 linewidth=2.5, linestyle="--", zorder=5)
            ax.add_patch(rect)

    # Overlap gene counts as secondary annotation
    overlap_counts = np.array([
        [0, 0, 8, 3],
        [6, 8, 0, 0],
        [1, 0, 0, 0],
        [4, 2, 0, 0],
    ])
    for i in range(4):
        for j in range(4):
            if overlap_counts[i, j] > 0:
                ax.text(j, i + 0.3, f"({overlap_counts[i,j]} genes)",
                        ha="center", va="center",
                        fontsize=9, color="#8b949e", fontstyle="italic")

    # Legend
    highlight = mpatches.Patch(facecolor="none", edgecolor="#f0883e",
                                linestyle="--", linewidth=2, label="Best match pair")
    ax.legend(handles=[highlight], loc="upper right", fontsize=11,
              facecolor="#161b22", edgecolor="#30363d", labelcolor="#f0f6fc")

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, "fig3_jaccard_heatmap.png")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[OK] Saved: {path}")


# ==============================================================================
# Run all
# ==============================================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Generating paper figures ...")
    print("=" * 60)
    draw_fig1()
    draw_fig2()
    draw_fig3()
    print("\n[OK] All figures generated in:", OUTPUT_DIR)
