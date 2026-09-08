# ============================================================================
# step9_disease_propagation_rwrm.R
# Random Walk with Restart on Multiplex (RWR-M) — Disease PPI Propagation
# Adapted from teacher's RWR-M code (Valdeolivas et al., 2017)
#
# Same methodology as step3, but runs on the disease PPI network from step8
# instead of the full STRING background network.
#
# Run in RStudio:  source("step9_disease_propagation_rwrm.R")
# ============================================================================

rm(list = ls())

library(igraph)
library(Matrix)

# ── 0. Configuration ─────────────────────────────────────────────────────────
project_dir <- "C:/Users/anike/OneDrive/Documents/Prob stats project/project phase 1 final/Project"

# RWR parameters (from teacher's defaults)
r     <- 0.7    # Restart probability (0.7 = stay close to seeds)
delta <- 0.5    # Inter-layer jump probability (irrelevant for single layer)
L     <- 1      # Number of layers (monoplex = 1)
tau   <- 1      # Layer restart weight

# ============================================================================
# FUNCTIONS (from teacher's RWR-M / All_Functions.R)
# ============================================================================

# --- Check which seed genes are present in the network ---
check.seeds <- function(Seeds, All_proteins) {
  Genes_Seeds_OK <- Seeds[which(Seeds %in% All_proteins)]
  Genes_Seeds_KO <- Seeds[which(!Seeds %in% All_proteins)]

  if (length(Genes_Seeds_OK) == 0) {
    stop("Seeds not found in our network")
  } else {
    if (length(Genes_Seeds_KO) > 0) {
      cat("Some seed genes not present in the network:\n")
      for (i in seq_along(Genes_Seeds_KO)) {
        cat("  ", Genes_Seeds_KO[i], "\n")
      }
    }
    return(Genes_Seeds_OK)
  }
}

# --- Build supra-adjacency matrix for multiplex network ---
get.supra.adj.multiplex <- function(Layers, delta, N) {
  Idem_Matrix <- Diagonal(N, x = 1)
  L <- length(Layers)

  SupraAdjacencyMatrix <- Matrix(0, ncol = N * L, nrow = N * L, sparse = TRUE)

  Col_Node_Names <- character()
  Row_Node_Names <- character()

  for (i in 1:L) {
    Adjacency_Layer <- as_adjacency_matrix(Layers[[i]], sparse = TRUE)

    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),
                                       order(colnames(Adjacency_Layer))]
    Layer_Col_Names <- paste(colnames(Adjacency_Layer), i, sep = "_")
    Layer_Row_Names <- paste(rownames(Adjacency_Layer), i, sep = "_")
    Col_Node_Names  <- c(Col_Node_Names, Layer_Col_Names)
    Row_Node_Names  <- c(Row_Node_Names, Layer_Row_Names)

    pos_ini <- 1 + (i - 1) * N
    pos_end <- N + (i - 1) * N
    SupraAdjacencyMatrix[pos_ini:pos_end, pos_ini:pos_end] <- (1 - delta) * Adjacency_Layer

    for (j in 1:L) {
      if (j != i) {
        pos_ini_col <- 1 + (j - 1) * N
        pos_end_col <- N + (j - 1) * N
        SupraAdjacencyMatrix[pos_ini:pos_end, pos_ini_col:pos_end_col] <-
          (delta / (L - 1)) * Idem_Matrix
      }
    }
  }

  rownames(SupraAdjacencyMatrix) <- Row_Node_Names
  colnames(SupraAdjacencyMatrix) <- Col_Node_Names
  return(SupraAdjacencyMatrix)
}

# --- Prepare seed score vector ---
get.seed.scores <- function(Genes, Number_Layers, tau) {
  Seeds_Genes_Scores        <- numeric(length = length(Genes) * Number_Layers)
  Seed_Genes_Layer_Labeled  <- character(length = length(Genes) * Number_Layers)

  for (k in 1:Number_Layers) {
    for (j in seq_along(Genes)) {
      idx <- ((k - 1) * length(Genes)) + j
      Seed_Genes_Layer_Labeled[idx] <- paste(Genes[j], k, sep = "_")
      Seeds_Genes_Scores[idx]       <- (tau / Number_Layers) / length(Genes)
    }
  }
  Seeds_Score <- data.frame(
    Seeds_ID = Seed_Genes_Layer_Labeled,
    Score    = Seeds_Genes_Scores,
    stringsAsFactors = FALSE
  )
  return(Seeds_Score)
}

# --- Random Walk with Restart (core algorithm) ---
Random_Walk_Restart <- function(Network_Matrix, r, SeedGenes) {
  Threshold   <- 1e-10
  NetworkSize <- ncol(Network_Matrix)

  residue <- 1
  iter    <- 1

  prox_vector <- matrix(0, nrow = NetworkSize, ncol = 1)
  prox_vector[which(colnames(Network_Matrix) %in% SeedGenes[, 1])] <- SeedGenes[, 2]
  prox_vector    <- prox_vector / sum(prox_vector)
  restart_vector <- prox_vector

  while (residue >= Threshold) {
    old_prox_vector <- prox_vector
    prox_vector     <- (1 - r) * (Network_Matrix %*% prox_vector) + r * restart_vector
    residue         <- sqrt(sum((prox_vector - old_prox_vector)^2))
    iter            <- iter + 1
  }

  cat(sprintf("  RWR converged in %d iterations\n", iter))
  return(prox_vector)
}

# ============================================================================
# PIPELINE
# ============================================================================

# ── 1. Load Disease PPI network (from step8) ─────────────────────────────────
cat("Step 9 (RWR-M): Loading Disease PPI network...\n")
bg_df <- read.table(
  file.path(project_dir, "results/networks/ppi_disease_BC.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)
cat(sprintf("  Loaded %d edges\n", nrow(bg_df)))

# ── 2. Build igraph object ───────────────────────────────────────────────────
ppi_graph <- graph_from_data_frame(
  d        = bg_df[, c("nodeA", "nodeB")],
  directed = FALSE
)
ppi_graph <- simplify(ppi_graph, remove.multiple = TRUE, remove.loops = TRUE)
cat(sprintf("  igraph: %d nodes, %d edges\n", vcount(ppi_graph), ecount(ppi_graph)))

# ── 3. Get pool of nodes ─────────────────────────────────────────────────────
pool_nodes        <- V(ppi_graph)$name
pool_nodes_sorted <- sort(pool_nodes)
N                 <- length(pool_nodes_sorted)
cat(sprintf("  Pool of nodes: %d\n", N))

# ── 4. Load and validate seed genes ──────────────────────────────────────────
seeds_df <- read.table(
  file.path(project_dir, "results/genes/seeds_BC.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)
all_seeds   <- seeds_df$gene_symbol
valid_seeds <- check.seeds(all_seeds, pool_nodes_sorted)
cat(sprintf("  Seed genes provided : %d\n", length(all_seeds)))
cat(sprintf("  Seeds in network    : %d\n", length(valid_seeds)))

# ── 5. Build supra-adjacency matrix ─────────────────────────────────────────
cat("  Building supra-adjacency matrix...\n")
Layers_list <- list(ppi_graph)

Node_Names_Layer <- V(Layers_list[[1]])$name
Missing_Nodes    <- pool_nodes[which(!pool_nodes %in% Node_Names_Layer)]
if (length(Missing_Nodes) > 0) {
  Layers_list[[1]] <- add_vertices(Layers_list[[1]], length(Missing_Nodes), name = Missing_Nodes)
}

SupraAdjacencyMatrix <- get.supra.adj.multiplex(Layers_list, delta, N)

# ── 6. Column-normalize ─────────────────────────────────────────────────────
cat("  Normalizing adjacency matrix...\n")
col_sums <- Matrix::colSums(SupraAdjacencyMatrix, na.rm = FALSE, sparseResult = FALSE)
col_sums[col_sums == 0] <- 1
Supra_Adj_Matrix_Normalized <- t(t(SupraAdjacencyMatrix) / col_sums)

# ── 7. Prepare seeds and run RWR-M ──────────────────────────────────────────
tau_per_layer <- tau / L
Seeds_Score   <- get.seed.scores(valid_seeds, L, tau)

cat(sprintf("  Running RWR-M with r = %.1f ...\n", r))
Random_Walk_Results <- Random_Walk_Restart(Supra_Adj_Matrix_Normalized, r, Seeds_Score)

# ── 8. Extract ranking ──────────────────────────────────────────────────────
rank_global <- data.frame(
  gene_symbol = gsub("_1$", "", rownames(Random_Walk_Results)[1:N]),
  score       = as.numeric(Random_Walk_Results[1:N, 1]),
  stringsAsFactors = FALSE
)

rank_global <- rank_global[order(-rank_global$score, rank_global$gene_symbol), ]

cat(sprintf("\n  Top 15 genes by RWR-M score:\n"))
print(head(rank_global, 15))

# ── 9. Add seed genes back with score = 1.0, format output ──────────────────
candidates_df <- rank_global[!(rank_global$gene_symbol %in% valid_seeds), ]

seeds_out_df <- data.frame(
  gene_symbol = valid_seeds,
  score       = 1.0,
  stringsAsFactors = FALSE
)

results_df <- rbind(seeds_out_df, candidates_df)
results_df$rank    <- seq_len(nrow(results_df))
results_df$method  <- "rwr"
results_df$disease <- "BC"
results_df         <- results_df[, c("gene_symbol", "score", "rank", "method", "disease")]

cat(sprintf("\n  Seed genes (score=1.0): %d\n", length(valid_seeds)))
cat(sprintf("  Candidate genes      : %d\n", nrow(candidates_df)))
cat(sprintf("  Total genes          : %d\n", nrow(results_df)))

# ── 10. Save output ─────────────────────────────────────────────────────────
output_dir  <- file.path(project_dir, "results", "genes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "expanded_BC_disease_rwr.tsv")

write.table(
  results_df, output_file,
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat(sprintf("\n✅ Step 9 (RWR-M) complete!\n"))
cat(sprintf("   Output: %s\n", output_file))
cat(sprintf("   %d total genes (%d seeds + %d candidates)\n",
            nrow(results_df), length(valid_seeds), nrow(candidates_df)))
cat(sprintf("   Parameters: r=%.1f, delta=%.1f, L=%d\n", r, delta, L))
