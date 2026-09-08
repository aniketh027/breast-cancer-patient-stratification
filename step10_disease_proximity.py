# ============================================================================
# step10_disease_proximity.py
# Network Proximity Testing on Disease PPI (NetColoc)
# Same methodology as step4 but on the disease PPI network.
#
# References:
#   [1] Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
#       rate. JRSS-B, 57(1), 289–300.
#   [2] Seabold, S., & Perktold, J. (2010). Statsmodels. SciPy 2010.
#   [3] Wright, S. et al. (2021). NetColoc. GitHub.
# ============================================================================

import pandas as pd
import numpy as np
import networkx as nx
from scipy import stats
from statsmodels.stats.multitest import multipletests
from netcoloc import netprop, netprop_zscore
import os

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

NETWORK_FILE   = os.path.join(PROJECT_DIR, "results", "networks", "ppi_disease_BC.tsv")
SEEDS_FILE     = os.path.join(PROJECT_DIR, "results", "genes", "seeds_BC.tsv")
RWR_FILE       = os.path.join(PROJECT_DIR, "results", "genes", "expanded_BC_disease_rwr.tsv")
OUTPUT_FILE    = os.path.join(PROJECT_DIR, "results", "genes", "proximity_BC_disease_rwr.tsv")

FDR_THRESHOLD  = 0.05

os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

# ── 1. Load Disease PPI and build graph ──────────────────────────────────────
print("=" * 60)
print("STEP 10: Proximity Testing on Disease PPI (NetColoc)")
print("=" * 60)

edges = pd.read_csv(NETWORK_FILE, sep="\t")
G = nx.from_pandas_edgelist(edges, source="nodeA", target="nodeB")
int_nodes = list(G.nodes())
print(f"\nNetwork: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# ── 2. Precompute individual heats matrix ────────────────────────────────────
print("\nCalculating w_prime (normalized adjacency)...")
w_prime = netprop.get_normalized_adjacency_matrix(G, conserve_heat=True)
print("Computing individual heats matrix (one-time, ~2-5 mins)...")
w_double_prime = netprop.get_individual_heats_matrix(w_prime, alpha=0.5)
print("Done.")

# ── 3. Load seed genes ──────────────────────────────────────────────────────
seeds_df   = pd.read_csv(SEEDS_FILE, sep="\t")
seed_genes = list(set(seeds_df["gene_symbol"]) & set(int_nodes))
print(f"\nSeed genes in network: {len(seed_genes)} / {len(seeds_df)}")

# ── 4. Calculate per-gene proximity z-scores (degree-binned, 1000 perms) ────
print("Calculating proximity z-scores (1000 permutations)...")
z_scores, Fnew, Fnew_rand = netprop_zscore.calculate_heat_zscores(
    w_double_prime,
    int_nodes,
    dict(G.degree()),
    seed_genes,
    num_reps         = 1000,
    minimum_bin_size = 100
)
print("Done.")

nan_z = z_scores.isna().sum()
print(f"NaN z-scores in full network : {nan_z} / {len(z_scores)}")

# ── 5. Load RWR genes and compute raw p-values ──────────────────────────────
rwr_df    = pd.read_csv(RWR_FILE, sep="\t")
rwr_genes = set(rwr_df["gene_symbol"]) & set(int_nodes)

records = []
for gene in rwr_genes:
    if gene in z_scores.index:
        z = z_scores[gene]
        if np.isnan(z):
            continue
        rwr_row = rwr_df.loc[rwr_df["gene_symbol"] == gene]
        records.append({
            "gene_symbol" : gene,
            "z_score"     : z,
            "raw_pvalue"  : stats.norm.sf(z),
            "rwr_score"   : rwr_row["score"].values[0],
            "rwr_rank"    : rwr_row["rank"].values[0]
        })

df = pd.DataFrame(records).sort_values("z_score", ascending=False).reset_index(drop=True)

print(f"\nRWR genes with valid z-scores : {len(df)}")
print(f"NaN p-values                  : {df['raw_pvalue'].isna().sum()}")

# ── 6. Benjamini-Hochberg FDR correction ─────────────────────────────────────
df_valid = df.dropna(subset=["raw_pvalue"]).copy()
df_nan   = df[df["raw_pvalue"].isna()].copy()

print(f"Genes entering BH correction  : {len(df_valid)}")

pvals = np.clip(df_valid["raw_pvalue"].to_numpy(dtype=float), 0, 1)

reject, p_adj, _, _ = multipletests(pvals, alpha=FDR_THRESHOLD, method="fdr_bh")

df_valid["adj_pvalue_bh"] = p_adj
df_valid["significant"]   = reject

if len(df_nan) > 0:
    df_nan["adj_pvalue_bh"] = np.nan
    df_nan["significant"]   = False

df = pd.concat([df_valid, df_nan], ignore_index=True) \
       .sort_values("z_score", ascending=False) \
       .reset_index(drop=True)

# ── 7. Save ──────────────────────────────────────────────────────────────────
df.to_csv(OUTPUT_FILE, sep="\t", index=False)

# ── 8. Summary ───────────────────────────────────────────────────────────────
n_sig   = df["significant"].sum()
n_total = len(df)

print(f"\n{'='*60}")
print(f"PROXIMITY TESTING RESULTS (Disease PPI)")
print(f"{'='*60}")
print(f"  Total candidates tested              : {n_total}")
print(f"  Significant (adj_pvalue_bh < {FDR_THRESHOLD})   : {n_sig}")
print(f"  Not significant                      : {n_total - n_sig}")

print(f"\nTop 20 significant genes (sorted by adj_pvalue_bh):")
top20 = df[df["significant"]].sort_values("adj_pvalue_bh").head(20)
print(top20[[
    "gene_symbol", "z_score", "raw_pvalue", "adj_pvalue_bh", "rwr_score"
]].to_string(index=False))

print(f"\n✅ Step 10 complete — results saved to:\n   {OUTPUT_FILE}")
