# Breast Cancer Patient Stratification using Graph Analytics
![Python](https://img.shields.io/badge/Python-3.x-blue)
![NetworkX](https://img.shields.io/badge/NetworkX-Graph%20Analytics-green)
![License](https://img.shields.io/badge/License-MIT-yellow)

A Python-based computational pipeline for identifying breast cancer gene modules using protein-protein interaction (PPI) networks and graph analytics techniques. The project integrates biological interaction data from multiple sources and applies network propagation, community detection, and pathway enrichment to generate disease-relevant gene modules for downstream patient stratification.

---

## Overview

Patient stratification is an important step in precision medicine, helping identify biologically meaningful subgroups of patients that may respond differently to treatments.

This project builds an end-to-end graph analytics pipeline that:

- Integrates protein-protein interaction (PPI) datasets
- Expands disease-associated genes using Random Walk with Restart (RWR)
- Identifies significant disease modules through network proximity analysis
- Detects gene communities using Leiden clustering
- Performs pathway enrichment analysis
- Compares disease-specific interaction networks using clustering evaluation metrics

---

## Data Sources

- STRING v12.0 Human Protein Interaction Network
- IntAct (PSICQUIC)
- BioGRID
- GWAS Breast Cancer Reported Genes

---

## Pipeline

```
Breast Cancer Seed Genes
            │
            ▼
     STRING PPI Network
            │
            ▼
 Random Walk with Restart (RWR)
            │
            ▼
   Network Proximity Analysis
            │
            ▼
 Significant Gene Subgraph
            │
            ▼
     Leiden Clustering
            │
            ▼
    Pathway Enrichment
            │
            ▼
────────────────────────────────────
Disease-Specific PPI Construction
(IntAct + BioGRID)
            │
            ▼
Disease Network Propagation
            │
            ▼
 Disease Module Detection
            │
            ▼
Cluster Comparison
(ARI • NMI • Modularity)
            │
            ▼
Final Disease Modules
```

---

## Methodology

### Pipeline A — STRING-Based Discovery

- Seed gene ingestion
- STRING protein interaction network construction
- Random Walk with Restart (RWR)
- NetColoc proximity testing
- Significant subgraph extraction
- Leiden community detection
- Pathway enrichment analysis

### Pipeline B — Disease-Specific Network

- Disease-specific PPI construction using IntAct and BioGRID
- Network propagation
- Proximity analysis
- Disease subgraph extraction
- Leiden clustering
- Disease module comparison
- Pathway enrichment

---

## Repository Structure (To be updated soon)

```
step1_ingest.py
step2_string.py
step3_propagation.py
step4_proximity.py
step5_subgraph.py
step6_clustering.py
step7_enrichment.py

step8_disease_ppi.py
step9_disease_propagation.py
step10_disease_proximity.py
step11_disease_subgraph.py
step12_disease_clustering.py
step13_compare.py
step14_enrichment.py

Breast_Cancer_Patient_Stratification_Pipeline.ipynb
Literature_Survey.md
```

---

## Technologies

- Python
- NetworkX
- pandas
- NumPy
- SciPy
- StatsModels
- python-igraph
- LeidenAlg
- GSEApy
- Matplotlib

---

## Biological Resources

- STRING v12.0
- IntAct (PSICQUIC)
- BioGRID
- GWAS Catalog

---

## Graph Analytics Techniques

- Random Walk with Restart (RWR)
- Network Propagation
- NetColoc Proximity Analysis
- Leiden Community Detection
- Graph Subgraph Extraction
- Pathway Enrichment Analysis

---

## Evaluation

The pipeline evaluates disease modules using:

- Adjusted Rand Index (ARI)
- Normalized Mutual Information (NMI)
- Modularity

These metrics were used to iteratively refine the disease-specific graph pipeline and improve agreement between STRING-derived and disease-specific gene modules.

---

## Key Features

- Multi-source biological network integration
- End-to-end graph analytics workflow
- Disease-specific protein interaction analysis
- Automated module detection pipeline
- Structured feature generation for downstream machine learning
- Reproducible computational workflow

---

## Future Work

- Integrate transcriptomic and methylation datasets
- Explore Graph Neural Networks (GNNs)
- Extend patient stratification using multi-omics data
- Compare additional community detection algorithms

---
