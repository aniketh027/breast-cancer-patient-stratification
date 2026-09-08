# ============================================================================
# step12_disease_clustering.py
# Leiden Clustering directly on Disease PPI Network
#
# Since the disease PPI (IntAct + BioGRID via PSICQUIC) is already
# disease-specific, we skip RWR/proximity/subgraph and cluster directly.
# Uses:
#   - Full disease PPI from step8
#   - LCC only
#   - Minimum cluster size of 20 genes
#
# References:
#   [1] Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain
#       to Leiden: guaranteeing well-connected communities.
#       Scientific Reports, 9(1), 5233.
#       https://doi.org/10.1038/s41598-019-41695-z
#
#   [2] Traag, V. A. (2024). leidenalg: Leiden algorithm for community
#       detection. Python package.
#       https://github.com/vtraag/leidenalg
#
#   [3] Csárdi, G., Nepusz, T., et al. (2024). igraph: Network analysis
#       and visualization. R/Python package.
#       https://doi.org/10.5281/zenodo.7682609
#       https://github.com/igraph/python-igraph
# ============================================================================

import pandas as pd
import networkx as nx
import igraph as ig
import leidenalg
import os

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

SUBGRAPH_FILE  = os.path.join(PROJECT_DIR, "results", "networks", "disease_subgraph_BC_rwr.tsv")
PROXIMITY_FILE = os.path.join(PROJECT_DIR, "results", "genes",    "proximity_BC_disease_rwr.tsv")
OUTPUT_FILE    = os.path.join(PROJECT_DIR, "results", "modules",  "modules_BC_psicquic.tsv")

MIN_CLUSTER_SIZE = 20

os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

# ── 1. Load disease subgraph + proximity data ───────────────────────────────
print("=" * 60)
print("STEP 12: Leiden Clustering on Disease PPI (filtered)")
print("=" * 60)

edges = pd.read_csv(SUBGRAPH_FILE, sep='\t')
print(f"\n  Disease subgraph edges: {len(edges)}")

# Load RWR scores for edge weights
prox_df = pd.read_csv(PROXIMITY_FILE, sep='\t')
rwr_lookup = dict(zip(prox_df["gene_symbol"], prox_df["rwr_score"]))

# Only use significant genes
sig_mask = prox_df["significant"] == True
selected_nodes = set(prox_df.loc[sig_mask, "gene_symbol"])

# ── 2. Build graph ───────────────────────────────────────────────────────────
G = nx.Graph()
for _, row in edges.iterrows():
    a, b = row["nodeA"], row["nodeB"]
    if a in selected_nodes and b in selected_nodes:
        rwr_a = rwr_lookup.get(a, 0.0)
        rwr_b = rwr_lookup.get(b, 0.0)
        rwr_edge_weight = (rwr_a + rwr_b) / 2.0
        G.add_edge(a, b, weight=rwr_edge_weight)

print(f"  Graph nodes: {G.number_of_nodes()}")
print(f"  Graph edges: {G.number_of_edges()}")

# ── 3. Keep Largest Connected Component ──────────────────────────────────────
lcc_nodes = max(nx.connected_components(G), key=len)
G = G.subgraph(lcc_nodes).copy()

print(f"  LCC nodes: {G.number_of_nodes()}")
print(f"  LCC edges: {G.number_of_edges()}")

# ── 4. Convert to igraph ─────────────────────────────────────────────────────
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

# ── 5. Leiden clustering [1] ─────────────────────────────────────────────────
leiden_partition = leidenalg.find_partition(
    ig_graph,
    leidenalg.RBConfigurationVertexPartition,
    weights              = "weight",
    resolution_parameter = 0.8,
    seed                 = 42,
    n_iterations         = 50
)

print(f"\n  Clusters found: {len(leiden_partition)}")
print(f"  Modularity    : {leiden_partition.modularity:.4f}")

# ── 6. Extract results ───────────────────────────────────────────────────────
leiden_clusters = {
    reverse_mapping[node_id]: comm
    for node_id, comm in enumerate(leiden_partition.membership)
}

cluster_df = pd.DataFrame([
    {"gene": gene, "cluster": cluster}
    for gene, cluster in leiden_clusters.items()
])

# ── 7. Drop clusters < 20 genes ─────────────────────────────────────────────
valid_clusters = cluster_df["cluster"].value_counts()
valid_clusters = valid_clusters[valid_clusters >= MIN_CLUSTER_SIZE].index
cluster_df     = cluster_df[cluster_df["cluster"].isin(valid_clusters)]

cluster_df = cluster_df.sort_values(["cluster", "gene"]).reset_index(drop=True)

# ── 8. Summary ───────────────────────────────────────────────────────────────
print(f"\nCluster summary (clusters with >= {MIN_CLUSTER_SIZE} genes):")
if len(cluster_df) > 0:
    summary = cluster_df["cluster"].value_counts().sort_index()
    for c, n in summary.items():
        print(f"  Cluster {c}: {n} genes")
    print(f"\n  Total clusters: {len(summary)}")
    print(f"  Total genes:    {len(cluster_df)}")
else:
    print("  ⚠ No clusters with >= 20 genes")

# ── 9. Save ──────────────────────────────────────────────────────────────────
cluster_df.to_csv(OUTPUT_FILE, sep='\t', index=False)

print(f"\n✅ Step 12 complete — saved to:\n   {OUTPUT_FILE}")
