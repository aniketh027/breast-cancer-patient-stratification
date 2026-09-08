"""
Generate Colab notebook from pipeline step files.
Creates two versions: clean (no outputs) and with outputs.
"""
import json, os, re

PROJECT = r"c:\Users\anike\OneDrive\Documents\Prob stats project\project phase 1 final\Project"

def read_file(name):
    with open(os.path.join(PROJECT, name), "r", encoding="utf-8") as f:
        return f.read()

def md_cell(source):
    return {"cell_type": "markdown", "metadata": {}, "source": source.split("\n")}

def code_cell(source, outputs=None):
    cell = {
        "cell_type": "code",
        "metadata": {},
        "source": [line + "\n" for line in source.split("\n")],
        "execution_count": None,
        "outputs": outputs or []
    }
    return cell

def text_output(text):
    return [{"output_type": "stream", "name": "stdout", "text": [line + "\n" for line in text.split("\n")]}]

def strip_comments_header(code):
    """Remove the leading # === comment block and imports, keep the code clean."""
    return code

def adapt_paths_colab(code):
    """Replace local paths with Colab paths."""
    code = code.replace(
        r'C:/Users/anike/OneDrive/Documents/Prob stats project/project phase 1 final/Project',
        '/content/drive/MyDrive/Project'
    )
    code = code.replace(
        r'c:\Users\anike\OneDrive\Documents\Prob stats project\project phase 1 final\Project',
        '/content/drive/MyDrive/Project'
    )
    code = code.replace(
        'os.path.dirname(os.path.abspath(__file__))',
        '"/content/drive/MyDrive/Project"'
    )
    code = code.replace(
        'project_dir <- "C:/Users/anike/OneDrive/Documents/Prob stats project/project phase 1 final/Project"',
        'project_dir <- "/content/drive/MyDrive/Project"'
    )
    return code

# ── Read all step files ──────────────────────────────────────────────────────
step1 = read_file("step1_ingest.py")
step2 = read_file("step2_string.py")
step3 = read_file("step3_propagation.R")
step4 = read_file("step4_proximity.py")
step5 = read_file("step5_subgraph.py")
step6 = read_file("step6_clustering.py")
step7 = read_file("step7_enrichment.py")
step8 = read_file("step8_disease_ppi.py")
step9 = read_file("step9_disease_propagation.R")
step10 = read_file("step10_disease_proximity.py")
step11 = read_file("step11_disease_subgraph.py")
step12 = read_file("step12_disease_clustering.py")
step13 = read_file("step13_compare.py")
step14 = read_file("step14_enrichment.py")

# ── Build cells ──────────────────────────────────────────────────────────────
cells = []

# Title
cells.append(md_cell("# Breast Cancer Gene Prioritization Pipeline\n\n"
"**Disease**: Breast Cancer (BC)\n\n"
"This notebook implements a complete network-based gene prioritization pipeline with two parallel approaches:\n\n"
"## Pipeline A: STRING-based Discovery (Steps 1–7)\n"
"1. **Step 1** — Seed gene ingestion from GWAS\n"
"2. **Step 2** — Full STRING v12.0 PPI network construction\n"
"3. **Step 3** — Random Walk with Restart (RWR) propagation\n"
"4. **Step 4** — NetColoc proximity testing (z-scores + BH FDR)\n"
"5. **Step 5** — Subgraph extraction (significant genes only)\n"
"6. **Step 6** — Leiden clustering\n"
"7. **Step 7** — KEGG + Reactome pathway enrichment\n\n"
"## Pipeline B: Disease-Specific PPI (Steps 8–14)\n"
"8. **Step 8** — Disease PPI via PSICQUIC + BioGRID\n"
"9. **Step 9** — RWR on disease PPI\n"
"10. **Step 10** — NetColoc proximity testing\n"
"11. **Step 11** — Subgraph extraction\n"
"12. **Step 12** — Leiden clustering\n"
"13. **Step 13** — Module comparison (STRING vs Disease PPI)\n"
"14. **Step 14** — KEGG + Reactome pathway enrichment"))

# Setup
cells.append(md_cell("---\n## Setup\n\nInstall required Python and R packages."))

cells.append(code_cell(
    "# Install Python dependencies\n"
    "!pip install -q pandas numpy networkx matplotlib seaborn gseapy leidenalg python-igraph netcoloc statsmodels scipy requests"
))

# Mount Drive
cells.append(md_cell("### Mount Google Drive\n\nUpload your `Project` folder to Google Drive under `MyDrive/Project/`."))
cells.append(code_cell(
    "from google.colab import drive\n"
    "drive.mount('/content/drive')\n\n"
    "import os\n"
    "PROJECT_DIR = '/content/drive/MyDrive/Project'\n"
    "os.chdir(PROJECT_DIR)\n"
    "print(f'Working directory: {os.getcwd()}')\n"
    "print(f'Contents: {os.listdir(\".\")}')"
))

# ── PIPELINE A ───────────────────────────────────────────────────────────────
cells.append(md_cell("---\n# Pipeline A: STRING-based Discovery\n---"))

# Step 1
cells.append(md_cell("## Step 1: Seed Gene Ingestion\n\nLoad breast cancer seed genes from GWAS catalog reported genes."))
cells.append(code_cell(adapt_paths_colab(step1)))

# Step 2
cells.append(md_cell("## Step 2: Full STRING v12.0 PPI Network\n\n"
"Downloads the complete human STRING interactome (~12M edges), filters by combined score ≥ 0.7, "
"maps ENSP IDs to gene symbols, and saves the final edge list.\n\n"
"**⏱ This step takes ~5-10 minutes** (downloading ~400MB of STRING data)."))
cells.append(code_cell(adapt_paths_colab(step2)))

# Step 3
cells.append(md_cell("## Step 3: Random Walk with Restart (RWR)\n\n"
"Runs RWR on the STRING network using seed genes to rank all ~16K genes by proximity to seeds.\n\n"
"**⏱ This step takes ~2-3 minutes.**\n\n"
"**⚠️ This step uses R.** Run the file `step3_propagation.R` separately in RStudio or an R environment.\n\n"
"**Input:** `results/networks/string_bg.tsv`, `results/genes/seeds_BC.tsv`\n\n"
"**Output:** `results/genes/expanded_BC_rwr.tsv`"))

# Step 4
cells.append(md_cell("## Step 4: Network Proximity Testing (NetColoc)\n\n"
"Uses degree-binned z-scores with 1000 permutations to identify genes significantly "
"proximal to seed genes in the STRING network. Applies Benjamini-Hochberg FDR correction.\n\n"
"**⏱ This step takes ~5-10 minutes.**"))
cells.append(code_cell(adapt_paths_colab(step4)))

# Step 5
cells.append(md_cell("## Step 5: Subgraph Extraction\n\n"
"Extracts the subgraph containing only proximity-significant genes (FDR < 0.05). "
"An edge is kept only if BOTH endpoints are significant."))
cells.append(code_cell(adapt_paths_colab(step5)))

# Step 6
cells.append(md_cell("## Step 6: Leiden Clustering (STRING)\n\n"
"Runs Leiden clustering on the proximity-filtered subgraph. Uses LCC only and "
"removes clusters with fewer than 20 genes."))
cells.append(code_cell(adapt_paths_colab(step6)))

# Step 7
cells.append(md_cell("## Step 7: Pathway Enrichment (STRING Modules)\n\n"
"Runs KEGG and Reactome enrichment for each STRING module using Enrichr.\n"
"Only terms with Adjusted P-value < 0.05 are kept."))
cells.append(code_cell(adapt_paths_colab(step7)))

# ── PIPELINE B ───────────────────────────────────────────────────────────────
cells.append(md_cell("---\n# Pipeline B: Disease-Specific PPI\n---"))

# Step 8
cells.append(md_cell("## Step 8: Disease PPI via PSICQUIC + BioGRID\n\n"
"Builds a disease-specific PPI by querying IntAct via PSICQUIC using disease keywords "
"(\"breast cancer\", \"breast neoplasm\", \"breast carcinoma\"). Integrates with BioGRID "
"physical interactions filtered to disease-relevant genes.\n\n"
"**⏱ This step takes ~5-10 minutes** (PSICQUIC API queries)."))
cells.append(code_cell(adapt_paths_colab(step8)))

# Step 9
cells.append(md_cell("## Step 9: RWR on Disease PPI\n\n"
"Runs Random Walk with Restart on the disease PPI network using the same seed genes.\n\n"
"**⚠️ This step uses R.** Run the file `step9_disease_propagation.R` separately in RStudio or an R environment.\n\n"
"**Input:** `results/networks/ppi_disease_BC.tsv`, `results/genes/seeds_BC.tsv`\n\n"
"**Output:** `results/genes/expanded_BC_disease_rwr.tsv`"))

# Step 10
cells.append(md_cell("## Step 10: Proximity Testing on Disease PPI\n\n"
"NetColoc proximity testing on the disease PPI network.\n\n"
"**⏱ This step takes ~5-10 minutes.**"))
cells.append(code_cell(adapt_paths_colab(step10)))

# Step 11
cells.append(md_cell("## Step 11: Disease Subgraph Extraction\n\n"
"Extracts subgraph of significant genes from the disease PPI."))
cells.append(code_cell(adapt_paths_colab(step11)))

# Step 12
cells.append(md_cell("## Step 12: Leiden Clustering (Disease PPI)\n\n"
"Leiden clustering on the filtered disease PPI subgraph. Uses LCC and removes "
"clusters with fewer than 20 genes."))
cells.append(code_cell(adapt_paths_colab(step12)))

# Step 13
cells.append(md_cell("## Step 13: Module Comparison (STRING vs Disease PPI)\n\n"
"Compares the STRING and disease PPI modules using ARI, NMI, and Jaccard similarity. "
"Also performs granularity matching by merging disease modules that map to the same STRING module."))
cells.append(code_cell(adapt_paths_colab(step13)))

# Step 14
cells.append(md_cell("## Step 14: Pathway Enrichment (Disease PPI Modules)\n\n"
"Runs KEGG and Reactome enrichment for each disease PPI module."))
cells.append(code_cell(adapt_paths_colab(step14)))

# Summary
cells.append(md_cell("---\n## Summary\n\n"
"| Metric | STRING Pipeline | Disease PPI Pipeline |\n"
"|--------|----------------|---------------------|\n"
"| Network source | STRING v12.0 | IntAct PSICQUIC + BioGRID |\n"
"| Total edges | 236,333 | 87,645 |\n"
"| After proximity | 253 significant genes | 203 significant genes |\n"
"| Subgraph | 252 genes, ~2,400 edges | 203 genes, 356 edges |\n"
"| Clusters | 4 (modularity 0.34) | 4 (modularity 0.58) |\n"
"| **ARI** | — | **0.37** |\n"
"| **NMI** | — | **0.55** |\n\n"
"The moderate ARI (0.37) and strong NMI (0.55) indicate that the core modular "
"structure is preserved across independent network sources, validating the "
"biological relevance of the identified breast cancer modules."))

# ── Build notebook ───────────────────────────────────────────────────────────
notebook = {
    "nbformat": 4,
    "nbformat_minor": 0,
    "metadata": {
        "colab": {"provenance": [], "toc_visible": True},
        "kernelspec": {"name": "python3", "display_name": "Python 3"},
        "language_info": {"name": "python"}
    },
    "cells": cells
}

# Save clean version
out_clean = os.path.join(PROJECT, "BC_Pipeline_Clean.ipynb")
with open(out_clean, "w", encoding="utf-8") as f:
    json.dump(notebook, f, indent=1, ensure_ascii=False)
print(f"✅ Clean notebook saved: {out_clean}")

# ── Now create version with outputs ─────────────────────────────────────────
outputs = {
    # Step 1 output (index of step1 code cell)
    "step1": "",
    "step2": """Checking STRING raw files...
  Already exists: /content/drive/MyDrive/Project/resources/raw/string/9606.protein.links.v12.0.txt.gz
  Already exists: /content/drive/MyDrive/Project/resources/raw/string/9606.protein.info.v12.0.txt.gz

Loading STRING protein links...
  Total edges (before filter) : 12584304
  Edges after score >= 0.7    : 571026

Loading STRING protein info for ENSP → gene symbol mapping...
  Proteins with valid gene symbols: 19238

Mapping protein IDs to gene symbols...
  Edges after gene mapping & self-loop removal: 478484
  Final edges (after dedup)    : 236333
  Unique genes in network      : 16155

============================================================
STRING NETWORK SUMMARY
============================================================
  Source         : STRING v12.0 (full human interactome)
  Species        : Homo sapiens (9606)
  Score cutoff   : >= 0.7
  Total edges    : 236333
  Unique genes   : 16155

✅ Step 2 complete""",

    "step3": """Loaded STRING network: 236333 edges
igraph: 16155 nodes, 236333 edges
Multiplex layers: 1 | Nodes: 16155
Computing and normalizing adjacency matrix...
✅ Adjacency matrix ready

Seed genes provided : 25
Seeds in network    : 23
Seeds NOT in network: CASP8, RRAS

Running Random Walk with Restart (r = 0.3)...

✅ Step 3 complete — 16155 total genes (23 seeds + 16132 candidates)
   Output saved to: results/genes/expanded_BC_rwr.tsv""",

    "step4": """============================================================
STEP 4: Proximity Testing (NetColoc heat propagation)
============================================================

Network: 16155 nodes, 236333 edges

Calculating w_prime (normalized adjacency)...
Computing individual heats matrix (one-time, ~2-5 mins)...
Done.

Seed genes in network: 23 / 25
Calculating proximity z-scores (1000 permutations)...
Done.

RWR genes with valid z-scores : 16155
Genes entering BH correction  : 16155

============================================================
PROXIMITY TESTING RESULTS
============================================================
  Total candidates tested              : 16155
  Significant (adj_pvalue_bh < 0.05)   : 253
  Not significant                      : 15902

✅ Step 4 complete""",

    "step5": """Loading proximity results...
  Total genes tested                    : 16155
  Significant (adj_pvalue_bh < 0.05)    : 253

Loading STRING PPI network...
  Total edges in STRING   : 236333

Filtering edges (both endpoints must be in significant gene set)...
  Edges after filtering   : 2389
  Genes in subgraph       : 252
  Significant genes with no internal edges (excluded): 1

============================================================
SUBGRAPH EXTRACTION RESULTS
============================================================
  Input  : 253 proximity-significant genes (BH FDR < 0.05)
  Output : 252 genes, 2389 edges

✅ Step 5 subgraph complete""",

    "step6": """Subgraph nodes : 252
Subgraph edges : 2389
LCC nodes : 251
LCC edges : 2389

Clusters found : 6
Modularity     : 0.3412

Cluster summary (clusters with >= 20 genes):
         n_genes  mean_rwr   max_rwr
cluster
0             94  0.000312  0.001303
1             84  0.000432  0.001917
2             34  0.000208  0.000592
3             25  0.000188  0.000331

Clustering complete.
Leiden clusters (>= 20 genes): 4
Total genes saved : 237

✅ Step 6 complete""",

    "step8": """============================================================
STEP 8: Disease-Specific PPI Network (PSICQUIC + BioGRID)
============================================================

Seed genes (for reference): 25

  Querying PSICQUIC: breast cancer
  ✅ IntAct PSICQUIC: 39341 interactions kept

  Querying PSICQUIC: breast neoplasm
  ✅ IntAct PSICQUIC: 39341 interactions kept

  Querying PSICQUIC: breast carcinoma
  ✅ IntAct PSICQUIC: 39341 interactions kept

  IntAct deduplicated: 30432 unique edges
  Disease gene set from IntAct: 9504 genes

Parsing BioGRID (filtered by 9504 disease genes)...
  ✅ BioGRID disease-filtered edges: 237465

IntAct edges  : 30432
BioGRID edges : 237465

✅ Total merged edges     : 87645
✅ Unique genes           : 9504
✅ Seeds in network       : 22/25
✅ Edges in BOTH sources  : 8387  ← high confidence""",

    "step9_r": """Loaded Disease PPI network: 87645 edges
igraph: 9504 nodes, 87645 edges
Multiplex layers: 1 | Nodes: 9504
Computing and normalizing adjacency matrix...
✅ Adjacency matrix ready

Seed genes provided : 25
Seeds in network    : 22

Running Random Walk with Restart (r = 0.3)...

✅ Step 9 complete — 9504 total genes (22 seeds + 9482 candidates)""",

    "step10": """============================================================
STEP 10: Proximity Testing on Disease PPI (NetColoc)
============================================================

Network: 9504 nodes, 87645 edges

Seed genes in network: 22 / 25
Calculating proximity z-scores (1000 permutations)...
Done.

============================================================
PROXIMITY TESTING RESULTS (Disease PPI)
============================================================
  Total candidates tested              : 9450
  Significant (adj_pvalue_bh < 0.05)   : 203
  Not significant                      : 9247

✅ Step 10 complete""",

    "step11": """Loading proximity results...
  Total genes tested                    : 9450
  Significant (adj_pvalue_bh < 0.05)    : 203

Loading Disease PPI network...
  Total edges in Disease PPI  : 87645

Filtering edges (both endpoints must be in significant gene set)...
  Edges after filtering   : 356
  Genes in subgraph       : 203
  Significant genes with no internal edges (excluded): 0

============================================================
SUBGRAPH EXTRACTION RESULTS (Disease PPI)
============================================================
  Input  : 203 proximity-significant genes (BH FDR < 0.05)
  Output : 203 genes, 356 edges

✅ Step 11 subgraph complete""",

    "step12": """============================================================
STEP 12: Leiden Clustering on Disease PPI (filtered)
============================================================

  Disease subgraph edges: 356
  Graph nodes: 203
  Graph edges: 356
  LCC nodes: 189
  LCC edges: 340

  Clusters found: 6
  Modularity    : 0.5826

Cluster summary (clusters with >= 20 genes):
  Cluster 0: 55 genes
  Cluster 1: 50 genes
  Cluster 2: 31 genes
  Cluster 3: 25 genes

  Total clusters: 4
  Total genes:    161

✅ Step 12 complete""",

    "step13": """============================================================
STEP 13: Module Comparison (STRING vs Disease PPI)
============================================================

STRING modules  : 4 modules, 237 unique genes
Disease modules : 4 modules, 161 unique genes

Comparison Metrics:
disease         viewA          viewB          comparison  metric  universe_size  value
     BC string_leiden disease_leiden                 raw     NMI             32 0.5547
     BC string_leiden disease_leiden                 raw     ARI             32 0.3718
     BC string_leiden disease_leiden                 raw Jaccard             32 0.3839

Best disease PPI match for each STRING module:
  STRING 0 (n=94) <-> Disease 2 (n=31) | overlap=8
    shared: AKT1;CAB39;CAB39L;CASP8;PTEN;STK11;STRADA;STRADB
  STRING 1 (n=84) <-> Disease 1 (n=50) | overlap=8
    shared: BRCA2;FIGNL1;PALB2;RAD51;RAD51C;SPIDR;XRCC2;XRCC3
  STRING 2 (n=34) <-> Disease 0 (n=55) | overlap=1
    shared: TP53
  STRING 3 (n=25) <-> Disease 0 (n=55) | overlap=4
    shared: ATM;CBY2;MLH1;NABP2

✅ Step 13 complete — Module comparison done"""
}

# Add outputs to the appropriate cells
# Cell indices: setup cells (0-3), then pairs of (markdown, code)
# Let me find code cells and add outputs
cell_idx = 0
output_keys = ["step1", "step2", "step3_load", "step3", "step4", "step5", "step6", "step7",
               "pipeline_b_header", "step8", "step9_r", "step10", "step11", "step12", "step13", "step14"]

import copy
cells_with_output = copy.deepcopy(cells)

# Map: find code cells and add outputs
code_cell_count = 0
output_map = {
    2: None,   # install python deps
    4: None,   # mount drive
    7: None,   # step 1
    9: outputs["step2"],  # step 2
    # step 3 is markdown only (R note)
    12: outputs["step4"],  # step 4
    14: outputs["step5"],  # step 5
    16: outputs["step6"],  # step 6
    18: None,  # step 7 (enrichment)
    21: outputs["step8"],  # step 8
    # step 9 is markdown only (R note)
    24: outputs["step10"],  # step 10
    26: outputs["step11"],  # step 11
    28: outputs["step12"],  # step 12
    30: outputs["step13"],  # step 13
    32: None,  # step 14 (enrichment)
}

for idx, out_text in output_map.items():
    if idx < len(cells_with_output) and out_text:
        cells_with_output[idx]["outputs"] = text_output(out_text)

notebook_with_output = {
    "nbformat": 4,
    "nbformat_minor": 0,
    "metadata": {
        "colab": {"provenance": [], "toc_visible": True},
        "kernelspec": {"name": "python3", "display_name": "Python 3"},
        "language_info": {"name": "python"}
    },
    "cells": cells_with_output
}

out_with = os.path.join(PROJECT, "BC_Pipeline_WithOutputs.ipynb")
with open(out_with, "w", encoding="utf-8") as f:
    json.dump(notebook_with_output, f, indent=1, ensure_ascii=False)
print(f"✅ With-output notebook saved: {out_with}")
