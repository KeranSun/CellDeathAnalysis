# CellDeathAnalysis

<!-- badges: start -->
[![R-CMD-check](https://github.com/keransun/CellDeathAnalysis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/keransun/CellDeathAnalysis/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

**CellDeathAnalysis** is a comprehensive R package for analyzing cell death pathway-related gene expression in transcriptomic data. 

**Author:** Keran Sun  
**Email:** s1214844197@163.com

### Key Features

- **14 Cell Death Types**: Ferroptosis, Cuproptosis, Disulfidptosis, Pyroptosis, Necroptosis, Apoptosis, Autophagy, PANoptosis, NETosis, Parthanatos, Entosis, Oxeiptosis, Alkaliptosis, LDCD
- **Multiple Scoring Methods**: ssGSEA, GSVA, AUCell, Z-score, Mean
- **Survival Analysis**: Kaplan-Meier, Cox regression, Time-dependent ROC
- **Enrichment Analysis**: ORA, GSEA, Group comparison
- **Single-Cell Analysis**: Seurat/SCE integration, Cluster comparison
- **Machine Learning**: Random Forest, LASSO, XGBoost, Ensemble models
- **Interactive Shiny App**: No coding required

## Installation

```r
# Install from GitHub
devtools::install_github("keransun/CellDeathAnalysis")

# After first install, generate data files:
setwd("path/to/CellDeathAnalysis")
source("data-raw/build_data.R")
devtools::install()
```

## Quick Start

### Option 1: Interactive Shiny App (No Coding)

```r
library(CellDeathAnalysis)
launch_death_app()
```

### Option 2: R Script

```r
library(CellDeathAnalysis)

# Load example data
data(example_expr)
data(example_clinical)

# Calculate pathway scores
scores <- calculate_death_score(example_expr, method = "zscore")

# Visualize
plot_death_boxplot(scores, group = example_clinical$group)
plot_death_radar(scores, group = example_clinical$group)
```

## Main Modules

### 1. Gene Set Management
```r
list_death_pathways()                    # List all pathways
get_death_geneset("ferroptosis")         # Get genes
get_pathway_info("ferroptosis")          # Get details
```

### 2. Pathway Scoring
```r
scores <- calculate_death_score(expr, method = "zscore")
groups <- classify_by_score(scores$ferroptosis, method = "median")
```

### 3. Visualization
```r
plot_death_heatmap(scores)               # Heatmap
plot_death_boxplot(scores, group)        # Boxplot
plot_death_radar(scores, group)          # Radar chart
plot_pathway_correlation(scores)         # Correlation
```

### 4. Survival Analysis
```r
result <- death_survival(scores, time, status, pathway = "ferroptosis")
plot_survival(result)

batch_result <- batch_survival(scores, time, status)
plot_forest(batch_result)
```

### 5. Enrichment Analysis
```r
# Over-representation analysis
ora <- death_enrich_ora(gene_list)
plot_enrichment_bar(ora)

# Group comparison
compare <- death_compare_groups(expr, group)
plot_compare_volcano(compare)
```

### 6. Single-Cell Analysis
```r
# Score Seurat object
seurat_obj <- sc_death_score(seurat_obj, method = "aucell")

# Visualize
sc_plot_death_feature(seurat_obj, pathways = c("ferroptosis", "pyroptosis"))
sc_plot_cluster_heatmap(seurat_obj)
sc_plot_cluster_violin(seurat_obj)

# Compare clusters
sc_compare_clusters(seurat_obj, group_by = "seurat_clusters")
```

### 7. Machine Learning
```r
# Build prediction model
model <- death_build_model(scores, outcome, method = "rf")
print(model)
plot_model(model, type = "importance")

# Cross-validation
cv_result <- death_cv(scores, outcome, method = "rf", nfolds = 5)

# Survival prediction
surv_model <- death_surv_model(scores, time, status, method = "cox_lasso")
```

### 8. Shiny App
```r
launch_death_app()        # Full-featured app
launch_death_app_mini()   # Quick scoring app
```

## Function Reference

| Module | Function | Description |
|--------|----------|-------------|
| **Gene Sets** | `get_death_geneset()` | Get pathway genes |
| | `list_death_pathways()` | List all pathways |
| **Scoring** | `calculate_death_score()` | Calculate scores |
| | `classify_by_score()` | Classify samples |
| **Visualization** | `plot_death_heatmap()` | Heatmap |
| | `plot_death_boxplot()` | Boxplot comparison |
| | `plot_death_radar()` | Radar chart |
| **Survival** | `death_survival()` | KM + Cox analysis |
| | `batch_survival()` | Batch analysis |
| | `plot_survival()` | Survival curves |
| | `plot_forest()` | Forest plot |
| **Enrichment** | `death_enrich_ora()` | ORA analysis |
| | `death_enrich_gsea()` | GSEA analysis |
| | `death_compare_groups()` | Group comparison |
| **Single-Cell** | `sc_death_score()` | Score scRNA-seq |
| | `sc_plot_death_feature()` | UMAP visualization |
| | `sc_compare_clusters()` | Cluster comparison |
| **ML** | `death_build_model()` | Build model |
| | `death_predict()` | Make predictions |
| | `death_cv()` | Cross-validation |
| **Shiny** | `launch_death_app()` | Interactive app |

## Required Packages

### Core (auto-installed)
- ggplot2, dplyr, tidyr

### Optional (install as needed)
```r
# Scoring
BiocManager::install(c("GSVA", "AUCell", "fgsea"))

# Survival
install.packages(c("survival", "survminer", "timeROC"))

# Single-cell
BiocManager::install(c("Seurat", "SingleCellExperiment"))

# Machine learning
install.packages(c("randomForest", "glmnet", "xgboost", "e1071"))

# Shiny app
install.packages(c("shiny", "shinydashboard", "DT"))

# Visualization
BiocManager::install("ComplexHeatmap")
install.packages(c("pheatmap", "patchwork", "ggrepel"))
```

## Citation

If you use CellDeathAnalysis in your research, please cite:

```
Sun K (2024). CellDeathAnalysis: Comprehensive Analysis of Cell Death 
Pathways in Gene Expression Data. R package version 0.3.0.
https://github.com/keransun/CellDeathAnalysis
```

## Data Sources

- [FerrDb V3](http://www.zhounan.org/ferrdb/)
- [MSigDB](https://www.gsea-msigdb.org)
- [KEGG](https://www.kegg.jp)
- Primary literature (2012-2024)

## License

MIT License

## Contact

- **Author:** Keran Sun
- **Email:** s1214844197@163.com
- **GitHub:** https://github.com/keransun/CellDeathAnalysis
- **Issues:** https://github.com/keransun/CellDeathAnalysis/issues
