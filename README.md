# Breast Cancer Patient Stratification

A dual-pipeline framework for breast cancer gene prioritization and patient stratification using Protein-Protein Interaction (PPI) networks, Random Walk with Restart (RWR), and community detection.

## Overview

This project builds two independent PPI-based pipelines and compares their modular structure to identify robust, source-independent gene modules linked to breast cancer:

| | Pipeline A (STRING) | Pipeline B (Disease PPI) |
|---|---|---|
| **Network Source** | STRING v12.0 (score ≥ 0.7) | IntAct (PSICQUIC) + BioGRID |
| **Network Size** | 16,155 genes / 236,333 edges | 9,504 genes / 87,645 edges |
| **Significant Genes** | 253 (FDR < 0.05) | 203 (FDR < 0.05) |
| **Modules Found** | 4 (modularity 0.34) | 4 (modularity 0.58) |
| **Key Pathways** | mTOR, MAPK, DNA repair | Homologous recombination, Fanconi anemia |

**Cross-pipeline agreement:** ARI = 0.37, NMI = 0.55 (ARI = 0.51, NMI = 0.60 after granularity matching)

## Pipeline Steps

### Pipeline A: STRING-based Discovery
| Step | Script | Description |
|------|--------|-------------|
| 1 | `step1_ingest.py` | Load breast cancer seed genes from GWAS catalog |
| 2 | `step2_string.py` | Download & filter STRING v12.0 network |
| 3 | `step3_propagation.R` | Random Walk with Restart (RWR) on STRING |
| 4 | `step4_proximity.py` | NetColoc proximity testing (1000 permutations) |
| 5 | `step5_subgraph.py` | Extract significant subgraph (FDR < 0.05) |
| 6 | `step6_clustering.py` | Leiden clustering (resolution 0.8) |
| 7 | `step7_enrichment.py` | KEGG + Reactome enrichment (adj. p < 0.05) |

### Pipeline B: Disease-Specific PPI
| Step | Script | Description |
|------|--------|-------------|
| 8 | `step8_disease_ppi.py` | Build disease PPI from IntAct + BioGRID |
| 9 | `step9_disease_propagation.R` | RWR on disease PPI |
| 10 | `step10_disease_proximity.py` | Proximity testing on disease network |
| 11 | `step11_disease_subgraph.py` | Extract disease subgraph |
| 12 | `step12_disease_clustering.py` | Leiden clustering on disease PPI |
| 13 | `step13_compare.py` | Compare modules (ARI, NMI, Jaccard) |
| 14 | `step14_enrichment.py` | Enrichment on disease modules |

## Setup

### Python Dependencies
```bash
pip install pandas numpy networkx matplotlib seaborn gseapy leidenalg python-igraph netcoloc statsmodels scipy requests
```

### R Dependencies
```r
install.packages(c("igraph", "Matrix", "remotes"))
remotes::install_github("alberto-valdeolivas/RandomWalkRestartMH")
```

## Project Structure
```
Project/
├── step1_ingest.py          # Seed gene loading
├── step2_string.py          # STRING network construction
├── step3_propagation.R      # RWR (STRING)
├── step4_proximity.py       # Proximity testing (STRING)
├── step5_subgraph.py        # Subgraph extraction (STRING)
├── step6_clustering.py      # Leiden clustering (STRING)
├── step7_enrichment.py      # Enrichment (STRING modules)
├── step8_disease_ppi.py     # Disease PPI construction
├── step9_disease_propagation.R  # RWR (Disease PPI)
├── step10_disease_proximity.py  # Proximity testing (Disease)
├── step11_disease_subgraph.py   # Subgraph extraction (Disease)
├── step12_disease_clustering.py # Leiden clustering (Disease)
├── step13_compare.py        # Module comparison
├── step14_enrichment.py     # Enrichment (Disease modules)
├── Literature_Survey.md     # Literature review (32 references)
├── set of genes.txt         # Input seed genes
├── others/                  # Colab notebooks
│   ├── BC_Pipeline_Clean.ipynb
│   └── BC_Pipeline_WithOutputs.ipynb
├── resources/               # Raw data (downloaded at runtime)
│   ├── biogrid/
│   └── raw/
└── results/                 # Pipeline outputs
    ├── genes/               # RWR & proximity results
    ├── networks/            # PPI edge lists & subgraphs
    ├── modules/             # Clustering assignments
    ├── enrichment/          # Per-module enrichment CSVs
    └── figures/             # Generated plots
```

## Key References

- Szklarczyk et al. (2023). STRING database. *Nucleic Acids Research*
- Valdeolivas et al. (2019). RandomWalkRestartMH. *Bioinformatics*
- Rosenthal et al. (2023). NetColoc. *Nature Protocols*
- Traag et al. (2019). Leiden algorithm. *Scientific Reports*
- Barabasi et al. (2011). Network Medicine. *Nature Reviews Genetics*

See [`Literature_Survey.md`](Literature_Survey.md) for the full 32-reference literature review.
