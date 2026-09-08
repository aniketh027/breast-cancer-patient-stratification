# ============================================================================
# step9_disease_propagation.R
# Random Walk with Restart (RWR) on Disease PPI
# Same methodology as step3 but on the disease PPI network from step8.
#
# References:
#   [1] Valentín, A. et al. (2017). RandomWalkRestartMH: RWR on multiplex
#       and heterogeneous biological networks. Bioinformatics.
#       https://doi.org/10.1093/bioinformatics/btx043
# ============================================================================

library(igraph)
library(Matrix)
library(RandomWalkRestartMH)

# ── 0. Set working directory ─────────────────────────────────────────────────
project_dir <- "C:/Users/anike/OneDrive/Documents/Prob stats project/project phase 1 final/Project"

# ── 1. Load Disease PPI network ──────────────────────────────────────────────
bg_df <- read.table(
  file.path(project_dir, "results/networks/ppi_disease_BC.tsv"),
  sep = '\t', header = TRUE, stringsAsFactors = FALSE
)
cat(sprintf('Loaded Disease PPI network: %d edges\n', nrow(bg_df)))

# ── 2. Build igraph object with weights ───────────────────────────────────────
ppi_graph <- graph_from_data_frame(
  d        = bg_df[, c('nodeA', 'nodeB', 'weight')],
  directed = FALSE
)

ppi_graph <- simplify(
  ppi_graph,
  remove.multiple = TRUE,
  remove.loops    = TRUE,
  edge.attr.comb  = list(weight = 'max')
)
cat(sprintf('igraph: %d nodes, %d edges\n', vcount(ppi_graph), ecount(ppi_graph)))

# ── 3. Create Multiplex object (monoplex — 1 layer) ──────────────────────────
PPI_Multiplex <- create.multiplex(list(DiseasePPI = ppi_graph))
cat(sprintf('Multiplex layers: %d | Nodes: %d\n',
            PPI_Multiplex$Number_of_Layers,
            PPI_Multiplex$Number_of_Nodes_Multiplex))

# ── 4. Compute and normalize the adjacency matrix ────────────────────────────
cat('Computing and normalizing adjacency matrix...\n')
AdjMatrix     <- compute.adjacency.matrix(PPI_Multiplex)
AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
cat('✅ Adjacency matrix ready\n')

# ── 5. Load seed genes ───────────────────────────────────────────────────────
seeds_df <- read.table(
  file.path(project_dir, "results/genes/seeds_BC.tsv"),
  sep = '\t', header = TRUE, stringsAsFactors = FALSE
)

all_seeds     <- seeds_df$gene_symbol
network_nodes <- PPI_Multiplex$Pool_of_Nodes
valid_seeds   <- intersect(all_seeds, network_nodes)
missing_seeds <- setdiff(all_seeds, network_nodes)

cat(sprintf('\nSeed genes provided : %d\n', length(all_seeds)))
cat(sprintf('Seeds in network    : %d\n', length(valid_seeds)))
if (length(missing_seeds) > 0) {
  cat(sprintf('Seeds NOT in network: %s\n', paste(missing_seeds, collapse = ', ')))
}

# ── 6. Run RWR ───────────────────────────────────────────────────────────────
cat('\nRunning Random Walk with Restart (r = 0.3)...\n')
RWR_Results <- Random.Walk.Restart.Multiplex(
  x               = AdjMatrixNorm,
  MultiplexObject = PPI_Multiplex,
  Seeds           = valid_seeds,
  r               = 0.3,
  tau             = c(1)
)

cat('\nTop 15 ranked candidate genes:\n')
print(head(RWR_Results$RWRM_Results, 15))

# ── 7. Add seed genes back with max score, then save ─────────────────────────
candidates_df           <- as.data.frame(RWR_Results$RWRM_Results)
colnames(candidates_df) <- c('gene_symbol', 'score')

candidates_df <- candidates_df[!(candidates_df$gene_symbol %in% valid_seeds), ]

seeds_out_df <- data.frame(
  gene_symbol = valid_seeds,
  score       = 1.0,
  stringsAsFactors = FALSE
)

results_df <- rbind(seeds_out_df, candidates_df)

results_df$rank    <- seq_len(nrow(results_df))
results_df$method  <- 'rwr'
results_df$disease <- 'BC'
results_df         <- results_df[, c('gene_symbol', 'score', 'rank', 'method', 'disease')]

cat(sprintf('\nSeed genes added back: %d (score = 1.0)\n', length(valid_seeds)))
cat(sprintf('New candidate genes  : %d\n', nrow(candidates_df)))

output_dir <- file.path(project_dir, 'results', 'genes')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, 'expanded_BC_disease_rwr.tsv')
write.table(
  results_df,
  output_file,
  sep = '\t', row.names = FALSE, quote = FALSE
)

cat(sprintf('\n✅ Step 9 complete — %d total genes (%d seeds + %d candidates)\n',
            nrow(results_df), length(valid_seeds), nrow(candidates_df)))
cat(sprintf('   Output saved to: %s\n', output_file))
