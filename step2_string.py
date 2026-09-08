# ============================================================================
# step2_string.py
# Full STRING PPI Network Construction
#
# Downloads the complete STRING v12.0 human PPI network, filters by
# combined score >= 0.7, maps ENSP protein IDs to gene symbols, removes
# duplicates and self-loops, and saves the final edge list.
#
# This approach uses the FULL human interactome (not seed-gene-centric API)
# which is the standard approach in network medicine.
#
# Required files (auto-downloaded if not present):
#   resources/raw/string/9606.protein.links.v12.0.txt.gz
#   resources/raw/string/9606.protein.info.v12.0.txt.gz
# ============================================================================

import os
import requests
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_DIR  = os.path.dirname(os.path.abspath(__file__))

STRING_DIR   = os.path.join(PROJECT_DIR, "resources", "raw", "string")
LINKS_FILE   = os.path.join(STRING_DIR,  "9606.protein.links.v12.0.txt.gz")
INFO_FILE    = os.path.join(STRING_DIR,  "9606.protein.info.v12.0.txt.gz")

OUTPUT_DIR   = os.path.join(PROJECT_DIR, "results", "networks")
OUTPUT_FILE  = os.path.join(OUTPUT_DIR,  "string_bg.tsv")
PLOT_FILE    = os.path.join(OUTPUT_DIR,  "ppi_network_BC.png")

SCORE_CUTOFF = 0.7   # combined score threshold (STRING scores are 0–1 after /1000)

LINKS_URL = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
INFO_URL  = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"

os.makedirs(STRING_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Helper: download file if not already present ──────────────────────────────
def download_if_missing(url, dest_path):
    if os.path.exists(dest_path):
        print(f"  Already exists: {dest_path}")
        return
    print(f"  Downloading: {url}")
    print(f"  Saving to  : {dest_path}")
    response = requests.get(url, stream=True)
    response.raise_for_status()
    with open(dest_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)
    print(f"  Download complete.")

# ── 1. Download STRING files ───────────────────────────────────────────────────
print("Checking STRING raw files...")
download_if_missing(LINKS_URL, LINKS_FILE)
download_if_missing(INFO_URL,  INFO_FILE)

# ── 2. Load full STRING links and filter by score ─────────────────────────────
print("\nLoading STRING protein links...")
df = pd.read_csv(LINKS_FILE, sep=" ", compression="gzip")

# STRING scores are integers 0–1000 → divide by 1000 to get 0–1
df["combined_score"] = df["combined_score"] / 1000.0

print(f"  Total edges (before filter) : {len(df)}")
df = df[df["combined_score"] >= SCORE_CUTOFF]
print(f"  Edges after score >= {SCORE_CUTOFF}    : {len(df)}")

# ── 3. Strip species prefix from protein IDs (9606.ENSP... → ENSP...) ─────────
df["protein1"] = df["protein1"].str.split(".").str[1]
df["protein2"] = df["protein2"].str.split(".").str[1]
df = df.reset_index(drop=True)

# ── 4. Load protein info and build ENSP → gene symbol mapping ─────────────────
print("\nLoading STRING protein info for ENSP → gene symbol mapping...")
info = pd.read_csv(INFO_FILE, sep="\t", compression="gzip")
info["protein_id"] = info["#string_protein_id"].str.split(".").str[1]

# Build mapping dict
protein_to_gene = dict(zip(info["protein_id"], info["preferred_name"]))

# Remove entries where preferred_name is an ENSG/ENSP ID (not a real gene name)
protein_to_gene = {k: v for k, v in protein_to_gene.items() if "ENS" not in v}
print(f"  Proteins with valid gene symbols: {len(protein_to_gene)}")

# ── 5. Map protein IDs to gene symbols ────────────────────────────────────────
print("\nMapping protein IDs to gene symbols...")
edges = df.copy()
edges["geneA"] = edges["protein1"].map(protein_to_gene)
edges["geneB"] = edges["protein2"].map(protein_to_gene)

# Drop edges where either gene could not be mapped
edges = edges.dropna(subset=["geneA", "geneB"])

# Drop self-loops (same gene on both sides)
edges = edges[edges["geneA"] != edges["geneB"]]

print(f"  Edges after gene mapping & self-loop removal: {len(edges)}")

# ── 6. Build final edge list ───────────────────────────────────────────────────
gene_edges = edges[["geneA", "geneB", "combined_score"]].copy()
gene_edges = gene_edges.rename(columns={
    "combined_score": "weight",
    "geneA": "nodeA",
    "geneB": "nodeB"
})

# Remove duplicate edges (A-B and B-A are the same undirected edge)
gene_edges["sorted_pair"] = gene_edges.apply(
    lambda x: tuple(sorted([x["nodeA"], x["nodeB"]])), axis=1
)
gene_edges = gene_edges.drop_duplicates("sorted_pair")
gene_edges = gene_edges.drop(columns="sorted_pair")
gene_edges = gene_edges.reset_index(drop=True)

# Add metadata columns to match your pipeline format
gene_edges["source"]   = "STRING"
gene_edges["evidence"] = "combined_score"
gene_edges["version"]  = "v12.0"

print(f"  Final edges (after dedup)    : {len(gene_edges)}")
print(f"  Unique genes in network      : {len(set(gene_edges['nodeA']) | set(gene_edges['nodeB']))}")

# ── 7. Save edge list ─────────────────────────────────────────────────────────
gene_edges.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"\n✅ Saved: {OUTPUT_FILE}")

# ── 8. Quick plot of seed gene subnetwork ─────────────────────────────────────
# (plotting full network is too large — plot seed-gene neighbourhood instead)
print("\nGenerating seed gene network plot...")
seeds_df = pd.read_csv(
    os.path.join(PROJECT_DIR, "results", "genes", "seeds_BC.tsv"), sep="\t"
)
seed_genes = set(seeds_df["gene_symbol"].tolist())

# Filter edges to seed genes only for plotting
seed_edges = gene_edges[
    gene_edges["nodeA"].isin(seed_genes) & gene_edges["nodeB"].isin(seed_genes)
]

G_plot = nx.from_pandas_edgelist(seed_edges, source="nodeA", target="nodeB", edge_attr="weight")

plt.figure(figsize=(12, 12))
nx.draw_networkx(G_plot, node_size=200, with_labels=True, font_size=8, width=1, alpha=0.7)
plt.title("STRING PPI Network — Breast Cancer Seed Genes (score ≥ 0.7)")
plt.savefig(PLOT_FILE, dpi=150, bbox_inches="tight")
plt.close()
print(f"✅ Saved: {PLOT_FILE}")

# ── 9. Summary ────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"STRING NETWORK SUMMARY")
print(f"{'='*60}")
print(f"  Source         : STRING v12.0 (full human interactome)")
print(f"  Species        : Homo sapiens (9606)")
print(f"  Score cutoff   : >= {SCORE_CUTOFF}")
print(f"  Total edges    : {len(gene_edges)}")
print(f"  Unique genes   : {len(set(gene_edges['nodeA']) | set(gene_edges['nodeB']))}")
print(f"\n✅ Step 2 complete — saved to:\n   {OUTPUT_FILE}")