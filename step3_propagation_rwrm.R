# ============================================================================
# step3_propagation_multiplex.R
# Random Walk with Restart — Multiplex Edition (RWR-M)
# Merges your Step 3 pipeline with the Valdeolivas et al. (2017) RWR-M method.
#
# What changed vs the original step3_propagation.R:
#   - Supra-adjacency matrix is now built by hand (friend's approach) instead
#     of relying on RandomWalkRestartMH internals. This lets you add more
#     network layers later (pathway, co-expression, etc.) just by extending
#     Network_List.
#   - Inter-layer jump probability (delta) is now a tunable parameter.
#   - Restart probability raised to r = 0.7 (friend's setting, closer to seeds).
#   - Seed scoring distributes weight uniformly across all layers (tau/L per layer).
#   - Final gene scores are the geometric mean across layers, matching the paper.
#   - Seeds are still pinned to score = 1 and added back to the output, keeping
#     your downstream steps (4-9) fully compatible.
#
# References:
#   Valdeolivas, A. et al. (2019). Random Walk with Restart on Multiplex and
#   Heterogeneous Biological Networks. Bioinformatics, 35(3), 497-505.
#   https://doi.org/10.1093/bioinformatics/btx612
# ============================================================================

library(igraph)
library(Matrix)

# ── 0. Project root & parameters ─────────────────────────────────────────────
project_dir <- "C:/Users/anike/OneDrive/Documents/Prob stats project/project phase 1 final/Project"

# RWR-M parameters (from friend's Parameters_AD.txt, now applied to BC)
r     <- 0.7   # Global restart probability (friend uses 0.7; original step3 used 0.3)
delta <- 0.5   # Inter-layer jump probability (only matters when L > 1)
# tau: per-layer restart weight vector — must sum to L (one entry per layer).
# With a single layer, tau = c(1) means all restart weight goes to that layer.
# If you add a second layer later, change to e.g. c(1, 1).
tau   <- c(1)
K     <- 20    # Top-K genes to report in the console summary

cat('── RWR-M Parameters ─────────────────────────────\n')
cat(sprintf('  r     = %.2f  (restart probability)\n', r))
cat(sprintf('  delta = %.2f  (inter-layer jump)\n', delta))
cat(sprintf('  tau   = %s   (layer weights)\n', paste(tau, collapse=', ')))
cat(sprintf('  K     = %d    (top-K summary)\n', K))
cat('─────────────────────────────────────────────────\n\n')

# ── 1. Load network layers ────────────────────────────────────────────────────
# Network_List drives which layers are loaded.
# Currently monoplex (STRING only) — identical scope to original step3.
# To add layers later, append names here and add a matching read block below.
Network_List <- c("STRING_PPI")
L            <- length(Network_List)

cat(sprintf('Loading %d network layer(s)...\n', L))

load_layer <- function(layer_name, project_dir) {
  if (layer_name == "STRING_PPI") {
    path <- file.path(project_dir, "results/networks/string_bg.tsv")
    df   <- read.table(path, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    g    <- graph_from_data_frame(d = df[, c('nodeA', 'nodeB', 'weight')],
                                  directed = FALSE)
    g    <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                     edge.attr.comb = list(weight = 'max'))
    cat(sprintf('  [STRING_PPI] %d nodes, %d edges\n', vcount(g), ecount(g)))
    return(g)
    
    # ── Placeholder for future layers ─────────────────────────────────────────
    # } else if (layer_name == "PATHWAY") {
    #   path <- file.path(project_dir, "results/networks/pathways_BC.tsv")
    #   df   <- read.table(path, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    #   g    <- graph_from_data_frame(d = df[, c('nodeA', 'nodeB')], directed = FALSE)
    #   g    <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    #   cat(sprintf('  [PATHWAY] %d nodes, %d edges\n', vcount(g), ecount(g)))
    #   return(g)
    #
    # } else if (layer_name == "COEXPRESSION") {
    #   path <- file.path(project_dir, "results/networks/coexpression_BC.tsv")
    #   df   <- read.table(path, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    #   g    <- graph_from_data_frame(d = df[, c('nodeA', 'nodeB', 'weight')],
    #                                 directed = FALSE)
    #   g    <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
    #                    edge.attr.comb = list(weight = 'max'))
    #   cat(sprintf('  [COEXPRESSION] %d nodes, %d edges\n', vcount(g), ecount(g)))
    #   return(g)
    # ──────────────────────────────────────────────────────────────────────────
    
  } else {
    stop(sprintf('Unknown layer name: "%s". Add a read block in load_layer().', layer_name))
  }
}

List_Layers <- lapply(Network_List, load_layer, project_dir = project_dir)
names(List_Layers) <- Network_List

# ── 2. Build unified node pool across all layers ──────────────────────────────
# Every gene that appears in ANY layer goes into the pool.
pool_nodes        <- unique(unlist(lapply(List_Layers, function(g) V(g)$name)))
pool_nodes_sorted <- sort(pool_nodes)
N                 <- length(pool_nodes_sorted)
cat(sprintf('\nUnion node pool: %d genes across %d layer(s)\n', N, L))

# ── 3. Add missing nodes to each layer as isolated vertices ───────────────────
# The supra-adjacency matrix requires every layer to have identical node sets.
List_Layers_Full <- lapply(List_Layers, function(g) {
  missing <- pool_nodes_sorted[!(pool_nodes_sorted %in% V(g)$name)]
  if (length(missing) > 0) g <- add_vertices(g, length(missing), name = missing)
  return(g)
})

# Sanity check: all layers must now have exactly N nodes
layer_sizes <- sapply(List_Layers_Full, vcount)
if (!all(layer_sizes == N)) {
  stop('Node count mismatch across layers after padding. Check load_layer().')
}
cat('All layers padded to the same node set. ✅\n')

# ── 4. Build the supra-adjacency matrix (friend's approach) ───────────────────
# The supra-adjacency matrix is an (N*L) x (N*L) block matrix:
#   - Diagonal blocks: (1 - delta) * adjacency of each layer
#   - Off-diagonal blocks: (delta / (L-1)) * identity  [inter-layer jumps]
# When L = 1 there are no off-diagonal blocks; delta has no effect.
cat('\nBuilding supra-adjacency matrix...\n')

I_N <- Diagonal(N, x = 1)                          # N×N identity
SupraAdj <- Matrix(0, nrow = N * L, ncol = N * L, sparse = TRUE)

row_names <- character()
col_names <- character()

for (i in seq_len(L)) {
  A_i <- as_adjacency_matrix(List_Layers_Full[[i]], attr = "weight", sparse = TRUE)
  # Sort rows/cols alphabetically so all layers align on the same node order
  A_i <- A_i[order(rownames(A_i)), order(colnames(A_i))]
  
  layer_nodes <- paste(colnames(A_i), i, sep = '_')
  col_names   <- c(col_names, layer_nodes)
  row_names   <- c(row_names, layer_nodes)
  
  ri <- (1 + (i - 1) * N) : (N + (i - 1) * N)
  
  # Intra-layer block
  SupraAdj[ri, ri] <- (1 - delta) * A_i
  
  # Inter-layer coupling blocks
  for (j in seq_len(L)) {
    if (j != i) {
      rj <- (1 + (j - 1) * N) : (N + (j - 1) * N)
      SupraAdj[ri, rj] <- (delta / (L - 1)) * I_N
    }
  }
}

rownames(SupraAdj) <- row_names
colnames(SupraAdj) <- col_names

# Column-normalise (transition matrix)
col_sums           <- Matrix::colSums(SupraAdj)
col_sums[col_sums == 0] <- 1    # avoid division-by-zero for isolated nodes
SupraAdj_Norm      <- t(t(SupraAdj) / col_sums)
cat('Supra-adjacency matrix built and normalised. ✅\n')

# ── 5. Load and validate seed genes ──────────────────────────────────────────
seeds_df      <- read.table(
  file.path(project_dir, 'results/genes/seeds_BC.tsv'),
  sep = '\t', header = TRUE, stringsAsFactors = FALSE
)
all_seeds     <- seeds_df$gene_symbol
valid_seeds   <- intersect(all_seeds, pool_nodes_sorted)
missing_seeds <- setdiff(all_seeds, pool_nodes_sorted)

cat(sprintf('\nSeed genes provided : %d\n', length(all_seeds)))
cat(sprintf('Seeds in network    : %d\n', length(valid_seeds)))
if (length(missing_seeds) > 0) {
  cat(sprintf('Seeds NOT in network: %s\n', paste(missing_seeds, collapse = ', ')))
}

# ── 6. Build seed score vector (friend's distributed approach) ────────────────
# Each seed gene gets weight tau[layer] / n_seeds in the supra-vector for
# that layer. tau is normalised by L before use, matching the paper.
tau_norm <- tau / L     # per-layer weight (sums to 1 across all layers)

n_seeds     <- length(valid_seeds)
Seeds_Score <- data.frame(
  Seeds_ID = character(n_seeds * L),
  Score    = numeric(n_seeds * L),
  stringsAsFactors = FALSE
)

for (k in seq_len(L)) {
  for (j in seq_len(n_seeds)) {
    idx                    <- (k - 1) * n_seeds + j
    Seeds_Score$Seeds_ID[idx] <- paste(valid_seeds[j], k, sep = '_')
    Seeds_Score$Score[idx]    <- tau_norm[k] / n_seeds
  }
}

cat('\nSeed score vector built (uniform weight across seeds and layers).\n')

# ── 7. Run RWR-M (iterative diffusion — friend's implementation) ──────────────
cat(sprintf('\nRunning RWR-M (r = %.2f, convergence threshold = 1e-10)...\n', r))

NetworkSize    <- ncol(SupraAdj_Norm)
prox_vector    <- matrix(0, nrow = NetworkSize, ncol = 1)
seed_positions <- which(colnames(SupraAdj_Norm) %in% Seeds_Score$Seeds_ID)
prox_vector[seed_positions] <- Seeds_Score$Score[
  match(colnames(SupraAdj_Norm)[seed_positions], Seeds_Score$Seeds_ID)
]
prox_vector    <- prox_vector / sum(prox_vector)
restart_vector <- prox_vector

threshold <- 1e-10
residue   <- 1
iter      <- 0

while (residue >= threshold) {
  old_vec     <- prox_vector
  prox_vector <- (1 - r) * (SupraAdj_Norm %*% prox_vector) + r * restart_vector
  residue     <- sqrt(sum((prox_vector - old_vec)^2))
  iter        <- iter + 1
}
cat(sprintf('Converged in %d iterations. ✅\n', iter))

# ── 8. Aggregate scores: geometric mean across layers ────────────────────────
# Take scores from layer 1 (indices 1..N), then compute geometric mean over
# all L copies of each gene — same as friend's Geometric_Mean() in C++.
# Pure-R implementation so no Rcpp dependency is needed.
geometric_mean_scores <- function(full_vector, L, N) {
  sapply(seq_len(N), function(i) {
    layer_vals <- sapply(seq_len(L), function(k) full_vector[(k - 1) * N + i])
    prod(layer_vals) ^ (1 / L)
  })
}

gene_names <- gsub('_1$', '', rownames(SupraAdj_Norm)[seq_len(N)])
geo_scores <- geometric_mean_scores(as.vector(prox_vector), L, N)

rank_global <- data.frame(
  gene_symbol = gene_names,
  score       = geo_scores,
  stringsAsFactors = FALSE
)
rank_global <- rank_global[order(-rank_global$score, rank_global$gene_symbol), ]

cat('\nTop 15 ranked candidate genes (raw RWR-M scores, before seed pinning):\n')
print(head(rank_global, 15))

# ── 9. Pin seed genes to score = 1, rebuild ranked output ────────────────────
# Keep your original step3 output format: seeds at score=1, then candidates.
candidates_df <- rank_global[!(rank_global$gene_symbol %in% valid_seeds), ]

seeds_out_df <- data.frame(
  gene_symbol = valid_seeds,
  score       = 1.0,
  stringsAsFactors = FALSE
)

results_df         <- rbind(seeds_out_df, candidates_df)
results_df$rank    <- seq_len(nrow(results_df))
results_df$method  <- 'rwr_multiplex'
results_df$disease <- 'BC'
results_df         <- results_df[, c('gene_symbol', 'score', 'rank', 'method', 'disease')]

cat(sprintf('\nSeed genes pinned (score = 1.0) : %d\n', length(valid_seeds)))
cat(sprintf('Candidate genes                  : %d\n', nrow(candidates_df)))

# ── 10. Save output ──────────────────────────────────────────────────────────
output_dir  <- file.path(project_dir, 'results', 'genes')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, 'expanded_BC_rwr_multiplex.tsv')

write.table(results_df, output_file, sep = '\t', row.names = FALSE, quote = FALSE)
cat(sprintf('\n✅ Step 3 (multiplex) complete — %d total genes (%d seeds + %d candidates)\n',
            nrow(results_df), length(valid_seeds), nrow(candidates_df)))
cat(sprintf('   Output saved to: results/genes/expanded_BC_rwr_multiplex.tsv\n'))

# ── 11. Quick verification ───────────────────────────────────────────────────
df <- read.table(output_file, sep = '\t', header = TRUE)

cat(sprintf('\nTotal genes in output : %d\n', nrow(df)))
cat(sprintf('Genes with score = 1  : %d (seed genes)\n', sum(df$score == 1.0)))
cat(sprintf('Candidate genes       : %d\n', sum(df$score < 1.0)))
cat(sprintf('Columns               : %s\n', paste(colnames(df), collapse = ', ')))
cat(sprintf('\nTop %d genes (seeds first, then candidates by RWR-M score):\n', K))
print(head(df[, c('rank', 'gene_symbol', 'score')], K))
cat('\n✅ Output verified. Continue with Steps 4–9.\n')