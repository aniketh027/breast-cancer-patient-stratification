# ============================================================================
# step6_clustering.py
# Leiden Clustering on the Proximity-Filtered STRING Subgraph
#
# Loads the subgraph built in step5_subgraph.py and runs Leiden clustering.
# Clusters with fewer than 20 genes are dropped as too small to be meaningful.
# ============================================================================

import pandas as pd
import networkx as nx
import igraph as ig
import leidenalg
import os

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_DIR  = os.path.dirname(os.path.abspath(__file__))

SUBGRAPH_FILE  = os.path.join(PROJECT_DIR, "results", "networks", "string_subgraph_BC_rwr.tsv")
PROXIMITY_FILE = os.path.join(PROJECT_DIR, "results", "genes",   "proximity_BC_rwr.tsv")
OUTPUT_FILE    = os.path.join(PROJECT_DIR, "results", "modules", "modules_BC_leiden.tsv")

MIN_CLUSTER_SIZE = 20   # drop clusters smaller than this

os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
edges   = pd.read_csv(SUBGRAPH_FILE,  sep='\t')
rwr_top = pd.read_csv(PROXIMITY_FILE, sep='\t')

# RWR score lookup {gene -> rwr_score}
rwr_lookup = dict(zip(rwr_top["gene_symbol"], rwr_top["rwr_score"]))

# Only use proximity-significant genes
sig_mask = rwr_top["significant"] == True
selected_nodes = set(rwr_top.loc[sig_mask, "gene_symbol"])

# ── Build subgraph with RWR-based edge weights ────────────────────────────────
G = nx.Graph()

for _, row in edges.iterrows():
    a, b = row["nodeA"], row["nodeB"]
    if a in selected_nodes and b in selected_nodes:
        rwr_a = rwr_lookup.get(a, 0.0)
        rwr_b = rwr_lookup.get(b, 0.0)
        rwr_edge_weight = (rwr_a + rwr_b) / 2.0   # mean RWR score as edge weight
        G.add_edge(a, b, weight=rwr_edge_weight)

print(f"Subgraph nodes : {G.number_of_nodes()}")
print(f"Subgraph edges : {G.number_of_edges()}")

# ── Keep largest connected component ─────────────────────────────────────────
lcc_nodes = max(nx.connected_components(G), key=len)
G = G.subgraph(lcc_nodes).copy()

print(f"LCC nodes : {G.number_of_nodes()}")
print(f"LCC edges : {G.number_of_edges()}")

# ── Convert to igraph ─────────────────────────────────────────────────────────
node_list       = list(G.nodes())
mapping         = {node: i for i, node in enumerate(node_list)}
reverse_mapping = {i: node for node, i in mapping.items()}

edges_with_weight = [
    (mapping[u], mapping[v], data["weight"])
    for u, v, data in G.edges(data=True)
]

ig_graph = ig.Graph(
    n        = len(node_list),
    edges    = [(u, v) for u, v, _ in edges_with_weight],
    directed = False
)
ig_graph.es["weight"] = [w for _, _, w in edges_with_weight]
ig_graph.vs["name"]   = node_list

# ── Leiden clustering ─────────────────────────────────────────────────────────
leiden_partition = leidenalg.find_partition(
    ig_graph,
    leidenalg.RBConfigurationVertexPartition,
    weights          = "weight",
    resolution_parameter = 1.0,
    seed             = 42,
    n_iterations     = 50
)

print(f"\nClusters found : {len(leiden_partition)}")
print(f"Modularity     : {leiden_partition.modularity:.4f}")

# ── Extract results ───────────────────────────────────────────────────────────
leiden_clusters = {
    reverse_mapping[node_id]: comm
    for node_id, comm in enumerate(leiden_partition.membership)
}

cluster_df = pd.DataFrame([
    {
        "gene"     : gene,
        "cluster"  : cluster,
        "rwr_score": rwr_lookup.get(gene, 0.0),
    }
    for gene, cluster in leiden_clusters.items()
])

# ── Drop clusters with fewer than 20 genes ────────────────────────────────────
# FIX 1: threshold raised from 5 → 20 as instructed
valid_clusters = cluster_df["cluster"].value_counts()
valid_clusters = valid_clusters[valid_clusters >= MIN_CLUSTER_SIZE].index
cluster_df     = cluster_df[cluster_df["cluster"].isin(valid_clusters)]

cluster_df = cluster_df.sort_values(["cluster", "rwr_score"], ascending=[True, False])

print(f"\nCluster summary (clusters with >= {MIN_CLUSTER_SIZE} genes):")
print(cluster_df.groupby("cluster").agg(
    n_genes  = ("gene",      "count"),
    mean_rwr = ("rwr_score", "mean"),
    max_rwr  = ("rwr_score", "max")
).to_string())

# ── Save ──────────────────────────────────────────────────────────────────────
# FIX 2: save cluster_df (filtered) instead of df (all nodes in G)
cluster_df.to_csv(OUTPUT_FILE, sep='\t', index=False)

print(f"\nClustering complete.")
print(f"Leiden clusters (>= {MIN_CLUSTER_SIZE} genes): {cluster_df['cluster'].nunique()}")
print(f"Total genes saved : {len(cluster_df)}")
print(f"\n✅ Step 6 complete — saved to:\n   {OUTPUT_FILE}")