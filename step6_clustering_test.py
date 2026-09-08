# ============================================================================
# step6_clustering.py
# Leiden Clustering on the Proximity-Filtered STRING Subgraph
#
# Fixes vs original:
#   - Cluster IDs remapped to 0,1,2... after dropping small clusters
#     (original kept gapped Leiden IDs e.g. 0,3,7 which confuses step7)
#   - Eval summary CSV saved to disk (modularity, cluster counts)
#     (original only printed modularity to console — lost on re-run)
#
# Design choices kept (deliberately better than friend's script):
#   - RWR-derived edge weights: mean(rwr_a, rwr_b) per edge
#     (friend uses raw PPI weights; ours bake in disease relevance)
#   - Pre-filter to proximity-significant genes before clustering
#     (friend clusters all LCC nodes; ours removes low-signal genes first)
#   - resolution_parameter = 1.0  (correct for STRING-scale network)
#
# References:
#   [1] Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain
#       to Leiden: guaranteeing well-connected communities.
#       Scientific Reports, 9(1), 5233.
#       https://doi.org/10.1038/s41598-019-41695-z
# ============================================================================

import pandas as pd
import networkx as nx
import igraph as ig
import leidenalg
import os

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

SUBGRAPH_FILE  = os.path.join(PROJECT_DIR, "results", "networks", "string_subgraph_BC_rwr.tsv")
PROXIMITY_FILE = os.path.join(PROJECT_DIR, "results", "genes",   "proximity_BC_rwr.tsv")
OUTPUT_FILE    = os.path.join(PROJECT_DIR, "results", "modules", "modules_BC_leiden.tsv")
EVAL_FILE      = os.path.join(PROJECT_DIR, "results", "modules", "clustering_eval_BC.tsv")

MIN_CLUSTER_SIZE     = 20
RESOLUTION_PARAMETER = 1   # higher = more, smaller clusters (correct for large STRING network)
SEED                 = 42
N_ITERATIONS         = 50

os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

print("=" * 60)
print("STEP 6: Leiden Clustering on STRING Subgraph")
print("=" * 60)

# ── 1. Load data ──────────────────────────────────────────────────────────────
edges   = pd.read_csv(SUBGRAPH_FILE,  sep='\t')
rwr_top = pd.read_csv(PROXIMITY_FILE, sep='\t')

print(f"\n  Subgraph edges loaded : {len(edges)}")
print(f"  Genes in proximity file: {len(rwr_top)}")

# RWR score lookup {gene -> rwr_score}
rwr_lookup = dict(zip(rwr_top["gene_symbol"], rwr_top["rwr_score"]))

# Only cluster genes that passed the proximity significance test
sig_mask       = rwr_top["significant"] == True
selected_nodes = set(rwr_top.loc[sig_mask, "gene_symbol"])
print(f"  Proximity-significant genes: {len(selected_nodes)}")

# ── 2. Build subgraph with RWR-based edge weights ─────────────────────────────
# Edge weight = mean RWR score of the two endpoint genes.
# This makes community detection sensitive to disease relevance, not just
# raw PPI confidence — edges between high-scoring genes are weighted heavier.
G = nx.Graph()

for _, row in edges.iterrows():
    a, b = row["nodeA"], row["nodeB"]
    if a in selected_nodes and b in selected_nodes:
        rwr_edge_weight = (rwr_lookup.get(a, 0.0) + rwr_lookup.get(b, 0.0)) / 2.0
        G.add_edge(a, b, weight=rwr_edge_weight)

print(f"\n  Graph nodes : {G.number_of_nodes()}")
print(f"  Graph edges : {G.number_of_edges()}")

if G.number_of_nodes() == 0:
    raise RuntimeError("Graph is empty — check that SUBGRAPH_FILE and PROXIMITY_FILE overlap.")

# ── 3. Keep largest connected component ───────────────────────────────────────
lcc_nodes = max(nx.connected_components(G), key=len)
G         = G.subgraph(lcc_nodes).copy()

print(f"  LCC nodes   : {G.number_of_nodes()}")
print(f"  LCC edges   : {G.number_of_edges()}")

# ── 4. Convert to igraph ──────────────────────────────────────────────────────
node_list       = list(G.nodes())
mapping         = {node: i for i, node in enumerate(node_list)}
reverse_mapping = {i: node for node, i in mapping.items()}

edges_with_weight = [
    (mapping[u], mapping[v], data["weight"])
    for u, v, data in G.edges(data=True)
]

ig_graph              = ig.Graph(
    n        = len(node_list),
    edges    = [(u, v) for u, v, _ in edges_with_weight],
    directed = False
)
ig_graph.es["weight"] = [w for _, _, w in edges_with_weight]
ig_graph.vs["name"]   = node_list

# ── 5. Leiden clustering ──────────────────────────────────────────────────────
leiden_partition = leidenalg.find_partition(
    ig_graph,
    leidenalg.RBConfigurationVertexPartition,
    weights              = "weight",
    resolution_parameter = RESOLUTION_PARAMETER,
    seed                 = SEED,
    n_iterations         = N_ITERATIONS
)

n_raw      = len(leiden_partition)
modularity = leiden_partition.modularity
print(f"\n  Clusters found (raw) : {n_raw}")
print(f"  Modularity           : {modularity:.4f}")

# ── 6. Extract results ────────────────────────────────────────────────────────
leiden_clusters = {
    reverse_mapping[node_id]: comm
    for node_id, comm in enumerate(leiden_partition.membership)
}

cluster_df = pd.DataFrame([
    {
        "gene"      : gene,
        "cluster"   : cluster,
        "rwr_score" : rwr_lookup.get(gene, 0.0),
    }
    for gene, cluster in leiden_clusters.items()
])

# ── 7. Drop clusters smaller than MIN_CLUSTER_SIZE ───────────────────────────
sizes         = cluster_df["cluster"].value_counts()
valid_raw_ids = sizes[sizes >= MIN_CLUSTER_SIZE].index
cluster_df    = cluster_df[cluster_df["cluster"].isin(valid_raw_ids)].copy()

# ── 8. Remap cluster IDs to 0, 1, 2... (FIX: removes gaps after filtering) ───
# Sort by descending size so module_0 is always the largest module
size_order = (
    cluster_df["cluster"]
    .value_counts()
    .sort_values(ascending=False)
    .index
)
remap      = {old: new for new, old in enumerate(size_order)}
cluster_df["cluster"] = cluster_df["cluster"].map(remap)

# Sort: by cluster ID, then by descending RWR score within each cluster
cluster_df = cluster_df.sort_values(
    ["cluster", "rwr_score"], ascending=[True, False]
).reset_index(drop=True)

n_valid       = cluster_df["cluster"].nunique()
n_genes_valid = len(cluster_df)
n_dropped     = G.number_of_nodes() - n_genes_valid

print(f"\n  Clusters kept (>= {MIN_CLUSTER_SIZE} genes): {n_valid}")
print(f"  Genes in valid clusters            : {n_genes_valid}")
print(f"  Genes dropped (small clusters)     : {n_dropped}")

# ── 9. Cluster summary ────────────────────────────────────────────────────────
print(f"\n  Cluster summary:")
summary = cluster_df.groupby("cluster").agg(
    n_genes  = ("gene",      "count"),
    mean_rwr = ("rwr_score", "mean"),
    max_rwr  = ("rwr_score", "max")
)
print(summary.to_string())

# ── 10. Save outputs ──────────────────────────────────────────────────────────
cluster_df.to_csv(OUTPUT_FILE, sep='\t', index=False)

# Eval summary (FIX: saved to disk, not just printed)
eval_df = pd.DataFrame([{
    "step"              : "step6",
    "network"           : "STRING_subgraph",
    "disease"           : "BC",
    "resolution"        : RESOLUTION_PARAMETER,
    "modularity"        : round(modularity, 6),
    "n_clusters_raw"    : n_raw,
    "n_clusters_valid"  : n_valid,
    "n_genes_valid"     : n_genes_valid,
    "n_genes_dropped"   : n_dropped,
    "min_cluster_size"  : MIN_CLUSTER_SIZE,
}])
eval_df.to_csv(EVAL_FILE, sep='\t', index=False)

print(f"\n✅ Step 6 complete")
print(f"   Modules  → {OUTPUT_FILE}")
print(f"   Eval     → {EVAL_FILE}")