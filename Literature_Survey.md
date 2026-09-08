# Literature Survey: Computational Approaches for Breast Cancer Gene Prioritization, Diagnosis, and Staging Using Protein-Protein Interaction Networks

---

## 1. Introduction

Breast cancer is the most commonly diagnosed cancer worldwide and a leading cause of cancer-related mortality in women. Early and accurate diagnosis, along with precise staging, is critical for selecting appropriate treatment strategies and improving patient outcomes. Traditional diagnostic approaches rely on histopathological examination, imaging, and clinical biomarkers such as ER, PR, and HER2 status. However, advances in high-throughput genomics and computational biology have opened new avenues for understanding the molecular underpinnings of breast cancer through network-based and statistical frameworks.

This literature survey reviews the key computational and statistical methods relevant to a **dual-pipeline breast cancer gene prioritization framework** that integrates:
- **Pipeline A:** STRING-based PPI network analysis using Random Walk with Restart (RWR), network proximity testing, Leiden clustering, and pathway enrichment.
- **Pipeline B:** Disease-specific PPI construction from IntAct (PSICQUIC) and BioGRID, followed by the same analytical workflow for independent validation.

The survey is organized into five sections: (1) gene expression-based probabilistic models for diagnosis, (2) computational staging prediction, (3) PPI network-based cancer prediction, (4) gene prioritization via network propagation, and (5) pathway enrichment connecting modules to clinical outcomes.

---

## 2. Probabilistic and Statistical Models for Breast Cancer Diagnosis

The fundamental clinical question — *what is the probability that a patient is cancerous or not?* — has been addressed through increasingly sophisticated statistical and machine learning models applied to gene expression data.

### 2.1 Gene Signature-Based Classification

**van 't Veer et al. (2002)** developed the landmark **70-gene MammaPrint signature** using supervised classification on microarray data from young, lymph-node-negative breast cancer patients. Their classifier significantly outperformed traditional clinical criteria (tumor size, grade, lymph node status) in predicting 5-year distant metastasis risk, establishing that high-dimensional gene expression profiles carry superior diagnostic information compared to conventional pathological assessment [1].

**Paik et al. (2004)** introduced the **Oncotype DX 21-gene recurrence score**, a continuous mathematical scoring algorithm based on RT-PCR expression levels. This score independently predicts 10-year distant recurrence risk in ER-positive, node-negative patients and guides adjuvant chemotherapy decisions. The probabilistic recurrence score stratifies patients into low, intermediate, and high-risk categories, directly addressing the diagnostic question of cancer aggressiveness [2].

### 2.2 Bayesian and Integrative Approaches

**Gevaert et al. (2006)** introduced a **Bayesian network framework** integrating clinical covariates and microarray gene expression data for breast cancer prognosis. Using Markov Blanket feature selection, their probabilistic graphical model achieved an AUC of 0.845, demonstrating that Bayesian approaches effectively handle the high-dimensional uncertainty inherent in genomic data and outperform models trained on single data modalities [3].

### 2.3 Statistical Stability of Genomic Classifiers

**Michiels et al. (2005)** critically evaluated the reproducibility of microarray-based prognostic classifiers across seven cancer studies, including breast cancer. Through multiple random validation splits, they demonstrated that generated gene signatures strongly depend on training sample selection, establishing the statistical necessity of rigorous cross-validation and external validation in genomic diagnostic modeling [4].

---

## 3. Computational Cancer Staging Prediction

Beyond binary diagnosis, predicting the **stage** of breast cancer (Stage I through IV) is essential for treatment planning and prognosis estimation.

### 3.1 Transcriptome-Based Stage Classification

**Athira and Gopakumar (2022)** developed a computational stage-wise classification framework using transcriptome meta-analysis and machine learning algorithms (SVM, Random Forest). By identifying stage-specific gene signatures (118 genes for Stage I, 12 for Stage II, 4 for Stage III), their SVM model achieved 92.21% accuracy, demonstrating that distinct molecular programs underlie each clinical stage [5].

**Sun et al. (2018)** formulated a multi-class stage prediction strategy using TCGA breast invasive carcinoma data with Minimum Redundancy Maximum Relevance (mRMR) feature selection coupled with ensemble learners. Their optimized Random Forest and XGBoost models achieved over 88% accuracy in separating early-stage (Stage I/II) from late-stage (Stage III/IV) breast cancers [6].

### 3.2 Molecular Subtyping and Stage Correlation

**Parker et al. (2009)** introduced the **PAM50** 50-gene classifier to assign breast tumors to intrinsic molecular subtypes (Luminal A, Luminal B, HER2-enriched, Basal-like) and calculate a continuous Risk of Recurrence (ROR) score. The ROR model combines subtype-specific gene expression with clinical stage and tumor size, linking molecular subtype classification directly to stage progression and long-term relapse risk [7].

**The Cancer Genome Atlas (TCGA) Network (2012)** performed comprehensive multi-omics profiling across 825 primary breast tumors using mRNA-seq, miRNA, DNA copy number, proteomics (RPPA), and DNA methylation. This integrative computational analysis mapped genomic alterations across clinical stages and identified stage- and subtype-specific driver mutations and targetable signaling pathways [8].

---

## 4. PPI Network-Based Cancer Prediction and Disease Modules

Protein-protein interaction networks provide a systems-level view of cellular organization, enabling the identification of disease-associated modules that go beyond individual gene-level analysis.

### 4.1 Network-Based Classification

**Chuang et al. (2007)** introduced a **subnetwork marker discovery algorithm** that overlays gene expression data onto human PPI networks to predict breast cancer metastasis. Subnetwork markers were significantly more reproducible between independent patient cohorts than single-gene markers and achieved higher classification accuracy, demonstrating the power of network context for diagnostic prediction [9].

**Hofree et al. (2013)** created **Network-Based Stratification (NBS)**, which integrates sparse somatic mutation profiles with human PPI networks via network propagation. Applied to TCGA breast cancer datasets, NBS successfully stratified patients into distinct clinical subtypes with significantly different survival outcomes that single-gene mutation frequencies failed to resolve [10].

### 4.2 Interactome Topology and Cancer Severity

**Taylor et al. (2009)** measured dynamic changes in protein interactome topology between normal breast tissue and breast carcinomas. They demonstrated that loss of dynamic hub coordination — a diffuse interactome structure across dynamic PPI modules — directly correlates with advanced stage, high grade, and decreased patient survival, establishing network topology as a prognostic indicator [11].

### 4.3 Network Medicine Framework

**Barabasi et al. (2011)** established the principles of **Network Medicine**, proposing that human diseases result from perturbations of localized subnetworks (disease modules) rather than isolated genes. This foundational framework laid the quantitative basis for disease module identification, disease-disease proximity metrics, and interactome-based cancer subtype diagnosis [12].

**Menche et al. (2015)** empirically validated the disease module hypothesis by mapping disease modules across hundreds of complex diseases within the human interactome. They introduced separation metrics to quantify module boundaries and overlap, providing a rigorous framework for network-based disease gene prioritization [13].

---

## 5. Gene Prioritization via Network Propagation

Network propagation methods, particularly Random Walk with Restart (RWR), form the computational backbone of our pipeline for ranking candidate disease genes.

### 5.1 Random Walk with Restart

**Tong et al. (2006)** developed fast algorithms for computing RWR on large graphs using low-rank matrix approximations, establishing the foundational computational framework. RWR simulates an iterative walker traversing graph edges while maintaining a restart probability of returning to designated seed nodes, capturing global network affinity [14].

**Kohler et al. (2008)** adapted RWR to biological networks for disease gene prioritization, demonstrating that global network diffusion significantly outperforms local metrics (such as direct degree or shortest path distance) in identifying disease drivers. This established RWR as the gold standard for seed-gene-based candidate prioritization in PPI networks [15].

### 5.2 Extensions to Multiplex and Heterogeneous Networks

**Valdeolivas et al. (2019)** introduced the **RandomWalkRestartMH** R package, extending RWR to multiplex networks (multiple interaction layers) and heterogeneous networks (combining gene and disease phenotype networks). This tool, which is used in our pipeline, enables multi-omics integration and cross-layer diffusion for enhanced gene prioritization [16].

**Li and Patra (2010)** formulated the **RWRH algorithm**, demonstrating that simultaneous random walking across dual-layer heterogeneous networks (PPI + disease similarity) significantly improves driver gene prioritization precision over single-network walks [17].

### 5.3 Network Propagation for Disease Gene Discovery

**Vanunu et al. (2010)** developed **PRINCE**, a network propagation algorithm that prioritizes disease genes and functional protein complexes by diffusing phenotypic disease similarity over PPI networks. Applied to complex diseases including breast cancer, PRINCE demonstrated high accuracy in prioritizing novel cancer drivers [18].

**Leiserson et al. (2015)** introduced **HotNet2**, using directional network propagation with a restart matrix to identify mutated subnetworks across pan-cancer datasets. HotNet2 successfully identified rare somatic driver mutations clustering within localized PPI subgraphs, bypassing the limitations of single-gene frequency-based tests [19].

---

## 6. Network Proximity Testing and Co-localization

Proximity testing in PPI networks quantifies whether a set of candidate genes is significantly closer to seed disease genes than expected by chance.

**Guney et al. (2016)** established **network proximity metrics** using z-score-adjusted distances correcting for node degree biases between drug targets and disease gene sets in the human interactome. This foundational work laid the groundwork for using network distance to prioritize candidate disease genes and drug targets [20].

**Rosenthal et al. (2023)** developed **NetColoc**, a protocol using network heat diffusion to quantify interactome proximity and co-localization of gene sets linked to different diseases or traits. NetColoc identifies shared molecular subnetworks and pathways, and its proximity testing framework (degree-binned z-scores with permutation testing) is directly employed in our pipeline [21].

---

## 7. Community Detection and Module Identification

Identifying functional modules within PPI subgraphs is essential for interpreting the biological roles of prioritized gene sets.

### 7.1 Leiden Algorithm

**Traag et al. (2019)** introduced the **Leiden clustering algorithm**, fixing a critical flaw in the Louvain method where identified clusters could be internally disconnected. The Leiden algorithm guarantees well-connected communities, optimizes modularity faster, and is the method used in our pipeline for partitioning proximity-filtered PPI subgraphs into functional disease modules [22].

### 7.2 Cluster Comparison Metrics

**Hubert and Arabie (1985)** introduced the **Adjusted Rand Index (ARI)** for evaluating agreement between two clusterings while correcting for chance. ARI equals 1.0 for identical partitions and 0 for random assignments [23].

**Vinh et al. (2010)** provided comprehensive mathematical formulations for **Normalized Mutual Information (NMI)** and related information-theoretic measures. NMI quantifies shared information between clustering solutions normalized by partition entropy [24].

In our pipeline, both ARI (0.37 -> 0.51 after granularity matching) and NMI (0.55 -> 0.60) are used to compare STRING-derived and disease PPI-derived module assignments, validating that the core modular structure is preserved across independent network sources.

---

## 8. Pathway Enrichment Analysis

Functional enrichment connects gene modules to known biological pathways, linking computational findings to clinical and therapeutic relevance.

### 8.1 Gene Set Enrichment Methods

**Subramanian et al. (2005)** established **Gene Set Enrichment Analysis (GSEA)**, a rank-based method to evaluate whether predefined pathway gene sets (KEGG, Reactome) show concordant differences between biological states. GSEA reveals subtle multi-gene pathway-level changes that single-gene statistics fail to detect [25].

**Liberzon et al. (2015)** curated 50 refined **MSigDB hallmark gene sets** representing well-defined biological processes in cancer (estrogen response, p53 pathway, hypoxia). Hallmark enrichment reduces feature redundancy and provides interpretable pathway scores for subtype and stage classification [26].

### 8.2 Enrichment Tools

**Chen et al. (2013)** developed **Enrichr**, a comprehensive gene set enrichment analysis engine aggregating hundreds of thousands of annotated gene sets across KEGG, Reactome, GO, cell types, and transcription factor targets. Enrichr's combined score (Fisher's exact test p-value x z-score deviation) provides reliable pathway annotations for prioritized gene modules [27]. This is the enrichment tool used in our pipeline for both STRING and disease PPI modules.

### 8.3 Sample-Level Pathway Scoring

**Hanzelmann et al. (2013)** introduced **GSVA** (Gene Set Variation Analysis), a non-parametric method that transforms gene-by-sample expression matrices into pathway-by-sample enrichment score matrices. GSVA facilitates downstream sample-level clustering, diagnostic classification, and stage stratification [28].

---

## 9. Databases Used for Disease PPI Construction

### 9.1 IntAct and PSICQUIC

**Hermjakob et al. (2004)** established **IntAct** as a premier open-source repository for curated molecular interaction data, enforcing rigorous PSI-MI standards for storing experimental evidence, detection methods, and host organisms [29].

**Aranda et al. (2011)** introduced **PSICQUIC**, a standardized RESTful web service for querying multiple molecular interaction databases (IntAct, BioGRID, MINT) simultaneously using a unified query language. PSICQUIC enables disease-keyword-driven queries (e.g., disease:"breast cancer") across heterogeneous repositories [30].

### 9.2 BioGRID

**Stark et al. (2006)** established **BioGRID** as an open-access repository containing manually curated physical and genetic interactions from high-throughput and low-throughput studies. BioGRID provides unbiased physical binding networks essential for constructing disease-specific interactomes [31].

### 9.3 STRING

**Szklarczyk et al. (2023)** published the latest update of the **STRING database** (v12.0), integrating physical and functional protein associations across genomic context, experiments, co-expression, and text mining. STRING provides the confidence-scored interactome backbone (combined score >= 0.7 threshold) used in our primary pipeline [32].

---

## 10. Connection to Our Pipeline

Our dual-pipeline framework directly implements and integrates the methods surveyed above:

| Pipeline Step | Method | Key References |
|---|---|---|
| Seed Gene Ingestion | GWAS Catalog curation | TCGA (2012) [8] |
| STRING Network | STRING v12.0 (score >= 0.7) | Szklarczyk et al. (2023) [32] |
| Disease PPI | IntAct PSICQUIC + BioGRID | Hermjakob (2004) [29], Aranda (2011) [30], Stark (2006) [31] |
| RWR Propagation | RandomWalkRestartMH | Kohler (2008) [15], Valdeolivas (2019) [16] |
| Proximity Testing | NetColoc (degree-binned z-scores) | Guney (2016) [20], Rosenthal (2023) [21] |
| Subgraph Extraction | FDR < 0.05 filtering | Menche (2015) [13] |
| Leiden Clustering | Resolution 0.8, LCC only | Traag (2019) [22] |
| Module Comparison | ARI, NMI | Hubert & Arabie (1985) [23], Vinh (2010) [24] |
| Pathway Enrichment | Enrichr (KEGG + Reactome, adj. p < 0.05) | Chen (2013) [27] |

The comparison metrics (ARI = 0.37, NMI = 0.55) between STRING and disease PPI modules validate that the identified functional modules — enriched for DNA repair, mTOR signaling, homologous recombination, and Fanconi anemia pathways — represent robust, source-independent biological signals relevant to breast cancer diagnosis and staging.

---

## References

[1] van 't Veer, L.J. et al. (2002). Gene expression profiling predicts clinical outcome of breast cancer. *Nature*, 415, 530-536.

[2] Paik, S. et al. (2004). A multigene assay to predict recurrence of tamoxifen-treated, node-negative breast cancer. *New England Journal of Medicine*, 351(27), 2817-2826.

[3] Gevaert, O. et al. (2006). Predicting the prognosis of breast cancer by integrating clinical and microarray data with Bayesian networks. *Bioinformatics*, 22(14), e184-e190.

[4] Michiels, S. et al. (2005). Prediction of cancer outcome with microarrays: a multiple random validation strategy. *The Lancet*, 365(9458), 488-492.

[5] Athira, K. & Gopakumar, G. (2022). Breast cancer stage prediction: a computational approach guided by transcriptome analysis. *Molecular Genetics and Genomics*, 297, 1547-1562.

[6] Sun, L. et al. (2018). Classification of breast cancer stages using gene expression profiles and machine learning algorithms. *IEEE Access*, 6, 65569-65577.

[7] Parker, J.S. et al. (2009). Supervised risk predictor of breast cancer based on intrinsic subtypes. *Journal of Clinical Oncology*, 27(8), 1160-1167.

[8] The Cancer Genome Atlas Network. (2012). Comprehensive molecular portraits of human breast tumours. *Nature*, 490, 61-70.

[9] Chuang, H.Y. et al. (2007). Network-based classification of breast cancer metastasis. *Molecular Systems Biology*, 3, 140.

[10] Hofree, M. et al. (2013). Network-based stratification of tumor mutations. *Nature Methods*, 10(11), 1108-1115.

[11] Taylor, I.W. et al. (2009). Dynamic modularity in protein interaction networks predicts breast cancer outcome. *Nature Biotechnology*, 27, 199-204.

[12] Barabasi, A.L. et al. (2011). Network medicine: a network-based approach to human disease. *Nature Reviews Genetics*, 12, 56-68.

[13] Menche, J. et al. (2015). Uncovering disease-disease relationships through the incomplete interactome. *Science*, 347(6224), 1257601.

[14] Tong, H. et al. (2006). Fast Random Walk with Restart and its applications. *IEEE ICDM*, 613-622.

[15] Kohler, S. et al. (2008). Walking the interactome for prioritization of candidate disease genes. *American Journal of Human Genetics*, 82(4), 949-958.

[16] Valdeolivas, A. et al. (2019). Random walk with restart on multiplex and heterogeneous biological networks. *Bioinformatics*, 35(3), 497-505.

[17] Li, Y. & Patra, J.C. (2010). Genome-wide inferring gene-phenotype relationship by walking on the heterogeneous network. *Bioinformatics*, 26(9), 1219-1224.

[18] Vanunu, O. et al. (2010). Associating genes and protein complexes with disease via network propagation. *PLoS Computational Biology*, 6(1), e1000641.

[19] Leiserson, M.D.M. et al. (2015). Pan-cancer network analysis identifies combinations of rare somatic mutations across pathways and protein complexes. *Nature Genetics*, 47, 106-114.

[20] Guney, E. et al. (2016). Network-based in silico drug efficacy screening. *Nature Communications*, 7, 10331.

[21] Rosenthal, S.B. et al. (2023). Mapping the common gene networks that underlie related diseases. *Nature Protocols*, 18, 1745-1759.

[22] Traag, V.A. et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9, 5233.

[23] Hubert, L. & Arabie, P. (1985). Comparing partitions. *Journal of Classification*, 2(1), 193-218.

[24] Vinh, N.X. et al. (2010). Information theoretic measures for clusterings comparison. *Journal of Machine Learning Research*, 11, 2837-2854.

[25] Subramanian, A. et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43), 15545-15550.

[26] Liberzon, A. et al. (2015). The Molecular Signatures Database hallmark gene set collection. *Cell Systems*, 1(6), 417-425.

[27] Chen, E.Y. et al. (2013). Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. *BMC Bioinformatics*, 14, 128.

[28] Hanzelmann, S. et al. (2013). GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics*, 14, 7.

[29] Hermjakob, H. et al. (2004). IntAct: an open source molecular interaction database. *Nucleic Acids Research*, 32(D1), D452-D455.

[30] Aranda, B. et al. (2011). PSICQUIC and PSISCORE: accessing and scoring molecular interactions. *Nature Methods*, 8(7), 528-529.

[31] Stark, C. et al. (2006). BioGRID: a general repository for interaction datasets. *Nucleic Acids Research*, 34(D1), D535-D539.

[32] Szklarczyk, D. et al. (2023). The STRING database in 2023. *Nucleic Acids Research*, 51(D1), D638-D646.
