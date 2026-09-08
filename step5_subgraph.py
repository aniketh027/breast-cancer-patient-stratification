# ============================================================================
# step5_subgraph.py
# Subgraph Extraction for Proximity-Significant RWR Genes
#
# Builds a subgraph from the STRING PPI network using only genes that:
#   1. Passed the Benjamini-Hochberg FDR proximity test (adj_pvalue_bh < 0.05)
#   2. Have at least one edge in STRING where BOTH endpoints are in the
#      proximity-significant gene set.
#
# Rule: An edge (A, B) is included only if BOTH gene A AND gene B are present
# in the significant gene set. If either endpoint is absent, the edge is
# dropped. Genes with no surviving edges are excluded from the subgraph.
# ============================================================================

import pandas as pd
import networkx as nx
import os

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_DIR    = os.path.dirname(os.path.abspath(__file__))

NETWORK_FILE   = os.path.join(PROJECT_DIR, "results", "networks", "string_bg.tsv")
PROXIMITY_FILE = os.path.join(PROJECT_DIR, "results", "genes",   "proximity_BC_rwr.tsv")
OUTPUT_EDGES   = os.path.join(PROJECT_DIR, "results", "networks", "string_subgraph_BC_rwr.tsv")

# ── 1. Load proximity results and filter significant genes ────────────────────
# Uses the 'significant' column from step4 (adj_pvalue_bh < 0.05).
print("Loading proximity results...")
prox_df = pd.read_csv(PROXIMITY_FILE, sep='\t')

sig_df    = prox_df[prox_df['significant'] == True].copy()
sig_genes = set(sig_df['gene_symbol'].tolist())

print(f"  Total genes tested                    : {len(prox_df)}")
print(f"  Significant (adj_pvalue_bh < 0.05)    : {len(sig_genes)}")

# ── 2. Load STRING PPI network ────────────────────────────────────────────────
print("\nLoading STRING PPI network...")
net_df = pd.read_csv(NETWORK_FILE, sep='\t')
print(f"  Total edges in STRING   : {len(net_df)}")

# ── 3. Filter edges — BOTH endpoints must be in the significant gene set ──────
# Rule: edge (A, B) is kept only if A ∈ sig_genes AND B ∈ sig_genes.
# This ensures the subgraph contains only well-connected significant genes.
# Genes whose ALL neighbours are outside the significant set are automatically
# excluded because none of their edges survive this filter.
print("\nFiltering edges (both endpoints must be in significant gene set)...")
mask      = net_df['nodeA'].isin(sig_genes) & net_df['nodeB'].isin(sig_genes)
sub_edges = net_df[mask].copy()

print(f"  Edges after filtering   : {len(sub_edges)}")

# ── 4. Identify which genes actually appear in the subgraph ───────────────────
nodes_in_subgraph = set(sub_edges['nodeA']).union(set(sub_edges['nodeB']))
excluded = sig_genes - nodes_in_subgraph

print(f"  Genes in subgraph       : {len(nodes_in_subgraph)}")
print(f"  Significant genes with no internal edges (excluded): {len(excluded)}")
if excluded:
    print(f"    Excluded: {', '.join(sorted(excluded))}")

# ── 5. Save subgraph edge list ────────────────────────────────────────────────
os.makedirs(os.path.dirname(OUTPUT_EDGES), exist_ok=True)
sub_edges.to_csv(OUTPUT_EDGES, sep='\t', index=False)

# ── 6. Summary ───────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"SUBGRAPH EXTRACTION RESULTS")
print(f"{'='*60}")
print(f"  Input  : {len(sig_genes)} proximity-significant genes (BH FDR < 0.05)")
print(f"  Output : {len(nodes_in_subgraph)} genes, {len(sub_edges)} edges")
print(f"\n✅ Step 5 subgraph complete — saved to:\n   {OUTPUT_EDGES}")