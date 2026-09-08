# ============================================================================
# step14_enrichment.py
# Pathway Enrichment for Disease PPI Modules
#
# Runs KEGG and Reactome enrichment on each disease PPI Leiden module
# from step12. Produces per-module CSV files and bar plots, plus a
# combined enrichment summary across all modules.
#
# Notes:
#   1. Skips plotting/saving when a cluster has no enrichment results
#   2. Uses only Adjusted P-value < 0.05 for significance (removes
#      non-significant terms before saving/plotting)
#   3. Saves to per-module folders (no overwriting)
#
# References:
#   [1] Kanehisa, M. & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes
#       and Genomes. Nucleic Acids Research, 28(1), 27–30.
#       https://doi.org/10.1093/nar/28.1.27
#
#   [2] Jassal, B. et al. (2020). The Reactome pathway knowledgebase.
#       Nucleic Acids Research, 48(D1), D498–D503.
#       https://doi.org/10.1093/nar/gkz1031
# ============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

try:
    import gseapy as gp
except ImportError:
    raise ImportError("gseapy is not installed. Run: pip install gseapy")

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

MODULES_FILE   = os.path.join(PROJECT_DIR, "results", "modules", "modules_BC_psicquic.tsv")
ENRICHMENT_DIR = os.path.join(PROJECT_DIR, "results", "enrichment_disease")

GENE_SETS = ["KEGG_2021_Human", "Reactome_2022"]
P_CUTOFF  = 0.05

os.makedirs(ENRICHMENT_DIR, exist_ok=True)

print("=" * 60)
print("STEP 14: Pathway Enrichment for Disease PPI Modules")
print("=" * 60)

# ── 1. Load module assignments ───────────────────────────────────────────────
modules_df = pd.read_csv(MODULES_FILE, sep="\t")
if "module_id" not in modules_df.columns and "cluster" in modules_df.columns:
    modules_df["module_id"] = modules_df["cluster"]
module_ids = sorted(modules_df["module_id"].unique())

print(f"\nLoaded {len(modules_df)} genes across {len(module_ids)} modules")

# ── 2. Per-module enrichment ─────────────────────────────────────────────────
all_results = []

for mod_id in module_ids:
    mod_genes = modules_df[modules_df["module_id"] == mod_id]["gene"].tolist()

    print(f"\n{'─' * 50}")
    print(f"Module {mod_id}: {len(mod_genes)} genes")
    print(f"  Genes: {', '.join(sorted(mod_genes)[:15])}{'...' if len(mod_genes) > 15 else ''}")

    if len(mod_genes) < 3:
        print("  ⚠ Too few genes for enrichment — skipping")
        continue

    # ── Run Enrichr (no cutoff here — filter manually after) ─────────────────
    mod_dir = os.path.join(ENRICHMENT_DIR, f"module_{mod_id}")
    os.makedirs(mod_dir, exist_ok=True)

    try:
        enr = gp.enrichr(
            gene_list=mod_genes,
            gene_sets=GENE_SETS,
            organism="human",
            outdir=os.path.join(mod_dir, "enrichr_raw"),
            cutoff=1.0,        # get ALL terms, filter below
            no_plot=True,
        )
        enrich_df = enr.results
    except Exception as e:
        print(f"  ❌ Enrichr failed: {e}")
        continue

    # 1. Skip empty results
    if enrich_df is None or enrich_df.empty:
        print(f"  No enrichment results → skip")
        continue

    # 2. Use Adjusted P-value only — remove non-significant terms
    enrich_df = enrich_df[enrich_df["Adjusted P-value"] < P_CUTOFF].copy()

    if enrich_df.empty:
        print(f"  No significant pathways (Adjusted P-value < {P_CUTOFF}) → skip")
        continue

    # Add module metadata
    enrich_df["module_id"]   = mod_id
    enrich_df["module_size"] = len(mod_genes)
    all_results.append(enrich_df)

    # 3. Save per-module (no overwriting — each module has its own folder)
    enrich_df.to_csv(os.path.join(mod_dir, "enrichment_full.csv"), index=False)

    # ── Per-database CSV and plots ──────────────────────────────────────────
    for db_name in GENE_SETS:
        db_df = (enrich_df[enrich_df["Gene_set"] == db_name]
                 .sort_values("Adjusted P-value")
                 .head(10)
                 .copy())

        if db_df.empty:
            continue

        short_name = db_name.split("_")[0].lower()
        db_df.to_csv(os.path.join(mod_dir, f"top_{short_name}.csv"), index=False)

        # Plot
        db_df["-log10(padj)"] = -np.log10(
            db_df["Adjusted P-value"].replace(0, 1e-300))
        db_df["Term"] = db_df["Term"].apply(
            lambda x: x if len(str(x)) <= 55 else str(x)[:52] + "...")

        fig, ax = plt.subplots(figsize=(10, max(3, len(db_df) * 0.4)))
        color = "skyblue" if "KEGG" in db_name else "lightgreen"
        sns.barplot(data=db_df, y="Term", x="-log10(padj)", color=color, ax=ax)
        ax.set_title(f"Module {mod_id} — Top {db_name} Pathways")
        ax.set_xlabel("-log10(Adjusted P-value)")
        ax.set_ylabel("")
        plt.tight_layout()
        plt.savefig(os.path.join(mod_dir, f"top_{short_name}.png"), dpi=150)
        plt.close()

        # Print top 5
        print(f"\n  Top 5 {db_name}:")
        for _, row in db_df.head(5).iterrows():
            print(f"    {row['Term']:<50} padj={row['Adjusted P-value']:.2e}")

# ── 3. Combined enrichment results ───────────────────────────────────────────
if all_results:
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(
        os.path.join(ENRICHMENT_DIR, "all_modules_enrichment.csv"),
        index=False
    )
    print(f"\n{'=' * 60}")
    print(f"ENRICHMENT SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Modules with enrichment : {combined['module_id'].nunique()}")
    print(f"  Total significant terms : {len(combined)}")
    print(f"  Combined CSV saved      : {os.path.join(ENRICHMENT_DIR, 'all_modules_enrichment.csv')}")
else:
    print("\n⚠ No enrichment results found for any module")

print(f"\n✅ Step 14 complete — Disease PPI enrichment done")
print(f"   Results in: {ENRICHMENT_DIR}")
