# ============================================================================
# step7_enrichment.py
# Pathway Enrichment Analysis per Leiden Module
#
# For each module identified in step6_clustering.py, runs Enrichr enrichment
# against KEGG and Reactome databases and saves results and plots per module.
#
# Notes:
#   1. Skips plotting/saving when a cluster has no enrichment results
#   2. Uses only Adjusted P-value < 0.05 for significance (removes
#      non-significant terms before saving/plotting)
#   3. Saves to per-module folders (no overwriting)
#
# Input : results/modules/modules_BC_leiden.tsv
# Output: results/enrichment/module_*/  (one folder per module)
#         results/enrichment/all_modules_enrichment.csv  (combined)
# ============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

try:
    import gseapy as gp
except ImportError:
    raise ImportError("gseapy not installed. Run: pip install gseapy")

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_DIR    = os.path.dirname(os.path.abspath(__file__))

MODULES_FILE   = os.path.join(PROJECT_DIR, "results", "modules", "modules_BC_leiden.tsv")
ENRICHMENT_DIR = os.path.join(PROJECT_DIR, "results", "enrichment")

GENE_SETS      = ["KEGG_2021_Human", "Reactome_2022"]
FDR_CUTOFF     = 0.05
MIN_GENES      = 10
TOP_N          = 10

os.makedirs(ENRICHMENT_DIR, exist_ok=True)

# ── 1. Load clustering results ────────────────────────────────────────────────
print("Loading module assignments...")
modules_df = pd.read_csv(MODULES_FILE, sep='\t')
print(f"  Total genes loaded : {len(modules_df)}")
print(f"  Modules found      : {modules_df['cluster'].nunique()}")
print(f"  Columns            : {list(modules_df.columns)}")

# ── 2. Run enrichment per module ──────────────────────────────────────────────
all_results = []

for cluster_id in sorted(modules_df['cluster'].unique()):

    module_genes = (
        modules_df[modules_df['cluster'] == cluster_id]['gene']
        .dropna()
        .astype(str)
        .tolist()
    )

    print(f"\n{'='*60}")
    print(f"Module {cluster_id} — {len(module_genes)} genes")
    print(f"{'='*60}")

    if len(module_genes) < MIN_GENES:
        print(f"  Skipping — fewer than {MIN_GENES} genes")
        continue

    # Output folder for this module
    module_dir = os.path.join(ENRICHMENT_DIR, f"module_{cluster_id}")
    os.makedirs(module_dir, exist_ok=True)

    # ── Run Enrichr ───────────────────────────────────────────────────────────
    try:
        enr = gp.enrichr(
            gene_list   = module_genes,
            gene_sets   = GENE_SETS,
            organism    = "human",
            outdir      = os.path.join(module_dir, "enrichr_raw"),
            cutoff      = 1.0,        # get ALL terms, filter below
            no_plot     = True,
        )
        enrich_df = enr.results.copy()
    except Exception as e:
        print(f"  Enrichr failed for module {cluster_id}: {e}")
        continue

    # 1. Skip empty results
    if enrich_df is None or enrich_df.empty:
        print(f"  No enrichment results → skip")
        continue

    # 2. Use Adjusted P-value only — remove non-significant terms
    enrich_df = enrich_df[enrich_df["Adjusted P-value"] < FDR_CUTOFF].copy()

    if enrich_df.empty:
        print(f"  No significant pathways (Adjusted P-value < {FDR_CUTOFF}) → skip")
        continue

    # Add module info
    enrich_df['module'] = cluster_id
    enrich_df['n_genes_in_module'] = len(module_genes)

    # 3. Save per-module (no overwriting)
    enrich_df.to_csv(os.path.join(module_dir, "enrichment_full.csv"), index=False)

    # ── Extract top pathways ──────────────────────────────────────────────────
    top_kegg = (
        enrich_df[enrich_df['Gene_set'] == 'KEGG_2021_Human']
        .sort_values('Adjusted P-value')
        .head(TOP_N)
        .copy()
    )

    top_reactome = (
        enrich_df[enrich_df['Gene_set'] == 'Reactome_2022']
        .sort_values('Adjusted P-value')
        .head(TOP_N)
        .copy()
    )

    # Print summary
    print(f"\n  Top KEGG pathways:")
    if not top_kegg.empty:
        for _, row in top_kegg.head(5).iterrows():
            print(f"    {row['Term'][:60]:<60}  adj_p={row['Adjusted P-value']:.2e}")
    else:
        print("    None found")

    print(f"\n  Top Reactome pathways:")
    if not top_reactome.empty:
        for _, row in top_reactome.head(5).iterrows():
            print(f"    {row['Term'][:60]:<60}  adj_p={row['Adjusted P-value']:.2e}")
    else:
        print("    None found")

    # ── Save CSVs ─────────────────────────────────────────────────────────────
    top_kegg.to_csv(os.path.join(module_dir, "top_kegg.csv"), index=False)
    top_reactome.to_csv(os.path.join(module_dir, "top_reactome.csv"), index=False)

    # ── Plot KEGG (only if significant terms exist) ──────────────────────────
    if not top_kegg.empty:
        top_kegg['-log10(padj)'] = -np.log10(
            top_kegg['Adjusted P-value'].replace(0, 1e-300)
        )
        top_kegg['Term'] = top_kegg['Term'].apply(
            lambda x: x if len(x) <= 55 else x[:52] + '...'
        )
        plt.figure(figsize=(10, 6))
        sns.barplot(data=top_kegg, y='Term', x='-log10(padj)', color='skyblue')
        plt.title(f"Module {cluster_id} — Top {TOP_N} KEGG Pathways ({len(module_genes)} genes)")
        plt.xlabel('-log10(Adjusted P-value)')
        plt.ylabel('Pathway')
        plt.tight_layout()
        plt.savefig(os.path.join(module_dir, "top_kegg.png"), dpi=150)
        plt.close()

    # ── Plot Reactome (only if significant terms exist) ──────────────────────
    if not top_reactome.empty:
        top_reactome['-log10(padj)'] = -np.log10(
            top_reactome['Adjusted P-value'].replace(0, 1e-300)
        )
        top_reactome['Term'] = top_reactome['Term'].apply(
            lambda x: x if len(x) <= 55 else x[:52] + '...'
        )
        plt.figure(figsize=(10, 6))
        sns.barplot(data=top_reactome, y='Term', x='-log10(padj)', color='lightgreen')
        plt.title(f"Module {cluster_id} — Top {TOP_N} Reactome Pathways ({len(module_genes)} genes)")
        plt.xlabel('-log10(Adjusted P-value)')
        plt.ylabel('Pathway')
        plt.tight_layout()
        plt.savefig(os.path.join(module_dir, "top_reactome.png"), dpi=150)
        plt.close()

    # Collect for combined output
    all_results.append(enrich_df)
    print(f"\n  ✅ Module {cluster_id} saved to: {module_dir}")

# ── 3. Save combined results across all modules ──────────────────────────────
if all_results:
    combined_df = pd.concat(all_results, ignore_index=True)
    combined_df = combined_df.sort_values(['module', 'Adjusted P-value'])
    combined_path = os.path.join(ENRICHMENT_DIR, "all_modules_enrichment.csv")
    combined_df.to_csv(combined_path, index=False)
    print(f"\n{'='*60}")
    print(f"ENRICHMENT COMPLETE")
    print(f"{'='*60}")
    print(f"  Modules processed    : {modules_df['cluster'].nunique()}")
    print(f"  Total significant    : {len(combined_df)}")
    print(f"  Combined results     : {combined_path}")
    print(f"\n  Per-module folders:")
    for cluster_id in sorted(modules_df['cluster'].unique()):
        print(f"    results/enrichment/module_{cluster_id}/")
    print(f"\n✅ Step 7 complete")
else:
    print("\n⚠️  No enrichment results generated.")