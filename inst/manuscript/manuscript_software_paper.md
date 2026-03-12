# CellDeathAnalysis: An R Package for Comprehensive Analysis of Regulated Cell Death Pathways in Transcriptomic Data

## Authors

Keran Sun^1,*

^1 [Your Institution], [City], [Country]

*Corresponding author: Email: s1214844197@163.com

---

## Abstract

**Summary:** Regulated cell death (RCD) encompasses diverse molecular mechanisms including ferroptosis, cuproptosis, pyroptosis, and other programmed death pathways that play critical roles in cancer biology and therapeutic response. Despite growing interest in these pathways, no unified computational tool exists for their systematic analysis. Here we present CellDeathAnalysis, an R package providing curated gene sets for 14 cell death types comprising over 500 genes, multiple scoring algorithms (ssGSEA, GSVA, AUCell, Z-score), and integrated analysis modules for survival analysis, pathway enrichment, single-cell RNA-seq, and machine learning-based prediction. The package includes an interactive Shiny application enabling researchers without programming expertise to perform comprehensive analyses. CellDeathAnalysis streamlines cell death pathway analysis and facilitates biomarker discovery in cancer research.

**Availability:** https://github.com/keransun/CellDeathAnalysis

**Contact:** s1214844197@163.com

**Supplementary information:** Supplementary data are available online.

---

## 1. Introduction

Regulated cell death (RCD) represents genetically controlled cellular demise executed through specific molecular machinery (Galluzzi et al., 2018). Beyond classical apoptosis, numerous RCD modalities have been characterized, including ferroptosis (iron-dependent lipid peroxidation), cuproptosis (copper-induced protein aggregation), pyroptosis (inflammatory cell death), necroptosis (programmed necrosis), and others (Tang et al., 2019).

Recent discoveries have expanded the RCD landscape. Ferroptosis, first described in 2012, has emerged as a promising therapeutic target in drug-resistant cancers (Dixon et al., 2012; Lei et al., 2022). Cuproptosis was identified in 2022 as a distinct copper-dependent death mechanism (Tsvetkov et al., 2022). Most recently, disulfidptosis was characterized in 2023 as a death form triggered by disulfide stress (Liu et al., 2023). Understanding these pathways' expression patterns and clinical associations is crucial for developing novel therapeutic strategies.

However, analyzing cell death pathways faces several challenges: (1) gene sets are scattered across multiple databases (FerrDb, MSigDB, KEGG) and primary literature; (2) no unified tool provides integrated analysis across all major death types; (3) existing resources lack support for single-cell data and machine learning applications; and (4) most tools require programming expertise, limiting accessibility.

To address these gaps, we developed CellDeathAnalysis, a comprehensive R package that integrates curated gene sets, multiple scoring methods, and analysis modules into a unified framework with an interactive interface.

---

## 2. Implementation

### 2.1 Gene Set Curation

We compiled gene sets for 14 RCD types from authoritative sources (Table 1). Ferroptosis genes were obtained from FerrDb V3 (Zhou et al., 2023), cuproptosis genes from Tsvetkov et al. (2022), and disulfidptosis genes from Liu et al. (2023). Apoptosis, autophagy, and necroptosis gene sets were derived from MSigDB hallmark collections and KEGG pathways (Liberzon et al., 2015; Kanehisa & Goto, 2000). Additional pathways (pyroptosis, PANoptosis, NETosis, parthanatos, entosis, oxeiptosis, alkaliptosis, lysosome-dependent cell death) were curated from primary literature. All gene symbols were standardized to HGNC nomenclature. The final collection comprises 523 unique genes.

### 2.2 Scoring Methods

CellDeathAnalysis implements six pathway scoring algorithms:

- **ssGSEA**: Single-sample Gene Set Enrichment Analysis using ranked gene expression (Barbie et al., 2009)
- **GSVA**: Gene Set Variation Analysis for pathway activity estimation (Hänzelmann et al., 2013)
- **AUCell**: Area Under the Curve-based scoring optimized for single-cell data (Aibar et al., 2017)
- **Z-score**: Mean of z-score normalized expression values (fast computation)
- **Mean/Median**: Simple aggregation methods for quick analysis

### 2.3 Analysis Modules

The package provides five integrated analysis modules:

**Survival Analysis Module** integrates Kaplan-Meier analysis, Cox proportional hazards regression, batch analysis across pathways, time-dependent ROC curves, and forest plot visualization.

**Enrichment Analysis Module** implements over-representation analysis (ORA) using Fisher's exact test and Gene Set Enrichment Analysis (GSEA) via the fgsea package.

**Single-Cell Analysis Module** supports Seurat and SingleCellExperiment objects, providing pathway scoring, UMAP/tSNE visualization, cluster comparison, and pseudotime analysis.

**Machine Learning Module** includes classification models (Random Forest, LASSO, Ridge, Elastic Net, SVM, XGBoost, ensemble) and survival prediction models (Cox-LASSO, Random Survival Forest).

**Visualization Module** offers publication-ready plots including heatmaps, boxplots, radar charts, violin plots, and correlation matrices.

### 2.4 Interactive Application

An interactive Shiny application (`launch_death_app()`) provides a graphical interface for all major functions, enabling researchers without R programming experience to perform comprehensive analyses through point-and-click operations.

---

## 3. Features and Usage

### 3.1 Quick Start

```r
# Install
devtools::install_github("keransun/CellDeathAnalysis")
library(CellDeathAnalysis)

# Calculate pathway scores
scores <- calculate_death_score(expr_matrix, method = "zscore")

# Visualize
plot_death_heatmap(scores, group = sample_groups)
plot_death_boxplot(scores, group = sample_groups)

# Survival analysis
surv_result <- death_survival(scores, time, status, pathway = "ferroptosis")
plot_survival(surv_result)

# Launch interactive app
launch_death_app()
```

### 3.2 Gene Set Access

```r
# List all pathways
list_death_pathways(detailed = TRUE)

# Get specific gene set
ferro_genes <- get_death_geneset("ferroptosis", type = "all")

# Get all death-related genes
all_genes <- get_all_death_genes()
```

### 3.3 Single-Cell Analysis

```r
# Score Seurat object
seurat_obj <- sc_death_score(seurat_obj, method = "aucell")

# Visualize on UMAP
sc_plot_death_feature(seurat_obj, pathways = c("ferroptosis", "pyroptosis"))

# Compare clusters
sc_compare_clusters(seurat_obj, group_by = "seurat_clusters")
```

### 3.4 Machine Learning

```r
# Build classification model
model <- death_build_model(scores, outcome, method = "rf")
print(model)  # Shows accuracy, AUC, etc.

# Cross-validation
cv_result <- death_cv(scores, outcome, method = "lasso", nfolds = 5)
```

---

## 4. Demonstration

To illustrate package functionality, we provide example datasets and demonstrate key analyses.

### 4.1 Package Contents

The package includes example expression data (`example_expr`) containing 5,000 genes across 100 samples (50 tumor, 50 normal) and corresponding clinical information (`example_clinical`) with survival data.

### 4.2 Pathway Scoring Demonstration

Using example data, we calculated death pathway scores and compared tumor versus normal samples (Figure 2). Multiple pathways showed differential activity between groups, demonstrating the package's ability to detect biological differences.

### 4.3 Survival Analysis Demonstration

We performed survival analysis on tumor samples using the example clinical data (Figure 3). The forest plot displays hazard ratios for all pathways, and Kaplan-Meier curves illustrate stratification by pathway scores. These demonstrations show the package's integrated survival analysis capabilities.

### 4.4 Correlation Analysis

Pathway correlation analysis revealed relationships between death mechanisms (Figure 4). Strong correlations between certain pathways (e.g., ferroptosis-autophagy, pyroptosis-necroptosis) reflect shared molecular components and biological crosstalk.

### 4.5 Enrichment Analysis

Over-representation analysis of example differentially expressed genes identified significantly enriched death pathways (Figure 5), demonstrating the enrichment module's functionality.

---

## 5. Comparison with Existing Tools

Table 2 compares CellDeathAnalysis with existing resources. CellDeathAnalysis is the only tool that: (1) covers all 14 major RCD types in one package; (2) provides integrated survival and enrichment analysis; (3) supports single-cell RNA-seq data; (4) includes machine learning models; and (5) offers an interactive graphical interface.

---

## 6. Application Scenarios

CellDeathAnalysis supports diverse research applications:

1. **Pan-cancer analysis**: Calculate death pathway scores across TCGA or other cancer cohorts to identify cancer-specific patterns.

2. **Biomarker discovery**: Use survival analysis to identify pathways associated with patient outcomes.

3. **Drug response prediction**: Correlate death pathway activities with drug sensitivity data.

4. **Single-cell characterization**: Analyze death pathway heterogeneity across cell populations in tumor microenvironments.

5. **Functional annotation**: Use enrichment analysis to interpret differentially expressed gene lists.

---

## 7. Conclusion

CellDeathAnalysis provides a comprehensive, user-friendly platform for analyzing cell death pathways in transcriptomic data. By integrating curated gene sets, multiple scoring methods, and diverse analysis modules, the package streamlines research workflows and enables systematic exploration of cell death mechanisms. The interactive Shiny application further extends accessibility to researchers without programming expertise. We anticipate CellDeathAnalysis will facilitate cell death research and contribute to therapeutic target discovery.

---

## Funding

[Add funding information if applicable]

## Conflict of Interest

None declared.

---

## References

Aibar S, González-Blas CB, Moerman T, et al. (2017) SCENIC: single-cell regulatory network inference and clustering. *Nat Methods*, 14:1083-1086.

Barbie DA, Tamayo P, Boehm JS, et al. (2009) Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*, 462:108-112.

Dixon SJ, Lemberg KM, Lamprecht MR, et al. (2012) Ferroptosis: an iron-dependent form of nonapoptotic cell death. *Cell*, 149:1060-1072.

Galluzzi L, Vitale I, Aaronson SA, et al. (2018) Molecular mechanisms of cell death: recommendations of the Nomenclature Committee on Cell Death 2018. *Cell Death Differ*, 25:486-541.

Hänzelmann S, Castelo R, Guinney J. (2013) GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics*, 14:7.

Kanehisa M, Goto S. (2000) KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Res*, 28:27-30.

Lei G, Zhuang L, Gan B. (2022) Targeting ferroptosis as a vulnerability in cancer. *Nat Rev Cancer*, 22:381-396.

Liberzon A, Birger C, Thorvaldsdóttir H, et al. (2015) The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst*, 1:417-425.

Liu X, Nie L, Zhang Y, et al. (2023) Actin cytoskeleton vulnerability to disulfide stress mediates disulfidptosis. *Nat Cell Biol*, 25:404-414.

Tang D, Kang R, Berghe TV, et al. (2019) The molecular machinery of regulated cell death. *Cell Res*, 29:347-364.

Tsvetkov P, Coy S, Petrova B, et al. (2022) Copper induces cell death by targeting lipoylated TCA cycle proteins. *Science*, 375:1254-1261.

Zhou N, Yuan X, Du Q, et al. (2023) FerrDb V2: update of the manually curated database of ferroptosis regulators and ferroptosis-disease associations. *Nucleic Acids Res*, 51:D483-D492.

---

## Tables

### Table 1. Cell Death Gene Sets in CellDeathAnalysis

| Cell Death Type | Year | Genes | Key Regulators | Source |
|-----------------|------|-------|----------------|--------|
| Ferroptosis | 2012 | 78 | GPX4, SLC7A11, ACSL4 | FerrDb V3 |
| Cuproptosis | 2022 | 23 | FDX1, LIAS, DLAT | Tsvetkov et al. |
| Disulfidptosis | 2023 | 20 | SLC7A11, NDUFS1 | Liu et al. |
| Pyroptosis | 2001 | 52 | NLRP3, CASP1, GSDMD | Literature |
| Necroptosis | 2005 | 48 | RIPK1, RIPK3, MLKL | MSigDB/KEGG |
| Apoptosis | 1972 | 89 | CASP3, BAX, BCL2 | MSigDB/KEGG |
| Autophagy | 1963 | 112 | BECN1, ATG5, LC3B | MSigDB/KEGG |
| PANoptosis | 2019 | 28 | ZBP1, CASP8 | Literature |
| NETosis | 2004 | 35 | MPO, ELANE, PADI4 | Literature |
| Parthanatos | 2009 | 18 | PARP1, AIFM1 | Literature |
| Entosis | 2007 | 15 | RHOA, ROCK1 | Literature |
| Oxeiptosis | 2018 | 8 | KEAP1, PGAM5 | Literature |
| Alkaliptosis | 2018 | 10 | CA9, SLC4A7 | Literature |
| LDCD | 2000 | 22 | LAMP1, CTSL | Literature |
| **Total** | - | **523** | - | - |

### Table 2. Comparison with Existing Tools

| Feature | CellDeathAnalysis | FerrDb | GSVA | AUCell |
|---------|-------------------|--------|------|--------|
| Death types covered | 14 | 1 | User-defined | User-defined |
| Curated gene sets | ✓ | ✓ | ✗ | ✗ |
| Multiple scoring methods | ✓ | ✗ | ✓ | ✗ |
| Survival analysis | ✓ | ✗ | ✗ | ✗ |
| Enrichment analysis | ✓ | ✗ | ✗ | ✗ |
| Single-cell support | ✓ | ✗ | Limited | ✓ |
| Machine learning | ✓ | ✗ | ✗ | ✗ |
| Interactive GUI | ✓ | Web only | ✗ | ✗ |
| R package | ✓ | ✗ | ✓ | ✓ |

---

## Figure Legends

**Figure 1.** CellDeathAnalysis package architecture and workflow. The package accepts expression data (bulk RNA-seq, single-cell RNA-seq, or microarray) and provides curated gene sets for 14 cell death types. Multiple scoring methods calculate pathway activity scores, which feed into analysis modules for survival analysis, enrichment analysis, single-cell analysis, and machine learning. An interactive Shiny application provides a graphical interface for all functions.

**Figure 2.** Gene set statistics and pathway scoring demonstration. (A) Bar plot showing the number of genes in each cell death pathway. (B) UpSet plot displaying gene overlap between pathways. (C) Heatmap of pathway scores across samples. (D) Boxplot comparing pathway scores between tumor and normal samples using example data.

**Figure 3.** Survival analysis demonstration. (A) Forest plot showing hazard ratios for all pathways. (B) Kaplan-Meier survival curves stratified by ferroptosis, cuproptosis, pyroptosis, and apoptosis scores. (C) Heatmap summarizing survival associations across pathways.

**Figure 4.** Pathway correlation analysis. (A) Correlation heatmap showing Spearman correlations between pathway scores. (B) Scatter plots for top correlated pathway pairs.

**Figure 5.** Enrichment analysis demonstration. (A) Bar plot showing over-representation analysis results. (B) Dot plot displaying fold enrichment and gene counts. (C) Volcano plot of pathway differential analysis between groups.

**Figure 6.** Interactive Shiny application interface. Screenshots showing (A) data upload page, (B) pathway scoring interface, (C) visualization options, and (D) survival analysis results.

---

## Supplementary Materials

**Table S1.** Complete gene list for all 14 cell death pathways with annotations.

**Table S2.** Scoring method comparison across different data types.

**Figure S1.** Comparison of scoring methods (ssGSEA, GSVA, AUCell, Z-score) using example data.

**Figure S2.** Shiny application detailed workflow and user guide.

---

## Software Availability

- **Package name:** CellDeathAnalysis
- **Version:** 0.3.0
- **Repository:** https://github.com/keransun/CellDeathAnalysis
- **Language:** R (≥ 4.0.0)
- **License:** MIT
- **Documentation:** Package vignettes and function help pages
- **Dependencies:** ggplot2, dplyr, tidyr (required); GSVA, AUCell, survival, Seurat, shiny (suggested)
