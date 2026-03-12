# CellDeathAnalysis 0.3.0

## New Features

### Single-Cell Analysis Module
* `sc_death_score()`: Calculate death scores for Seurat/SCE objects
  - Supports AUCell (recommended), Seurat-style, mean, and Z-score methods
  - Automatic integration with object metadata
* `sc_plot_death_feature()`: Visualize scores on UMAP/tSNE
* `sc_compare_clusters()`: Compare scores across cell clusters
* `sc_plot_cluster_heatmap()`: Heatmap of mean scores by cluster
* `sc_plot_cluster_violin()`: Violin plots by cluster
* `sc_identify_high_death()`: Identify cells with high pathway activity
* `sc_diff_death()`: Differential analysis between conditions
* `sc_plot_trajectory()`: Visualize scores along pseudotime

### Machine Learning Module
* `death_build_model()`: Build predictive models
  - Random Forest, LASSO, Ridge, Elastic Net, SVM, XGBoost, Ensemble
  - Automatic train/test split and evaluation
* `death_predict()`: Make predictions on new data
* `death_cv()`: K-fold cross-validation
* `death_surv_model()`: Survival prediction (Cox-LASSO, RSF)
* `plot_model()`: Visualize model performance and importance

### Shiny Interactive Application
* `launch_death_app()`: Full-featured interactive dashboard
  - Data upload (CSV/TSV/RDS)
  - Gene set exploration
  - Pathway scoring
  - Visualization (heatmap, boxplot, radar, correlation)
  - Survival analysis
  - Enrichment analysis
  - Machine learning
  - Result export
* `launch_death_app_mini()`: Simplified quick-scoring app

---

# CellDeathAnalysis 0.2.0

## New Features

### Survival Analysis Module
* `death_survival()`: Perform Kaplan-Meier and Cox regression analysis
* `plot_survival()`: Create publication-ready survival curves with risk tables
* `batch_survival()`: Analyze all pathways at once with multiple testing correction
* `plot_forest()`: Forest plot for hazard ratios across pathways
* `death_cox_multivariate()`: Multivariate Cox regression with clinical variables
* `death_time_roc()`: Time-dependent ROC analysis
* `plot_time_roc()`: Visualize time-dependent ROC curves

### Enrichment Analysis Module
* `death_enrich_ora()`: Over-representation analysis (Fisher's exact test)
* `death_enrich_gsea()`: Gene set enrichment analysis (requires fgsea)
* `plot_enrichment_bar()`: Bar plot for enrichment results
* `plot_enrichment_dot()`: Dot plot for enrichment results
* `plot_gsea_curve()`: Classic GSEA enrichment plot
* `death_compare_groups()`: Compare pathway enrichment between groups
* `plot_compare_volcano()`: Volcano plot for group comparisons

### Other Improvements
* Updated author information
* Added more example data features
* Improved documentation

---

# CellDeathAnalysis 0.1.0

## Initial Release

### New Features

* **Gene Sets**: Curated gene sets for 14 types of regulated cell death:
  - Ferroptosis
  - Cuproptosis
  - Disulfidptosis
  - Pyroptosis
  - Necroptosis
  - Apoptosis
  - Autophagy
  - PANoptosis
  - NETosis
  - Parthanatos
  - Entosis
  - Oxeiptosis
  - Alkaliptosis
  - LDCD (Lysosome-dependent cell death)

* **Gene Set Management**:
  - `get_death_geneset()`: Retrieve gene sets for specific pathways
  - `list_death_pathways()`: List all available pathways
  - `get_pathway_info()`: Get detailed pathway information
  - `get_all_death_genes()`: Get union of all genes
  - `check_pathway_overlap()`: Check overlap between pathways
  - `add_custom_geneset()`: Add user-defined gene sets

* **Scoring Methods**:
  - `calculate_death_score()`: Calculate pathway enrichment scores
    - Supports ssGSEA, GSVA, AUCell, Z-score, mean, and median methods
  - `score_pathway()`: Quick scoring for single pathway
  - `classify_by_score()`: Classify samples into high/low groups
  - `calculate_death_index()`: Calculate composite death index

* **Visualization**:
  - `plot_death_heatmap()`: Heatmap of pathway scores
  - `plot_death_boxplot()`: Boxplot comparison between groups
  - `plot_death_radar()`: Radar/spider chart
  - `plot_pathway_correlation()`: Correlation heatmap
  - `plot_pathway_genes()`: Gene expression heatmap
  - `plot_score_distribution()`: Score distribution plots
  - `plot_pathway_network()`: Pathway network visualization

### Data

* `death_genesets`: Complete gene sets with sub-categories
* `death_genesets_simple`: Simplified core gene sets
* `death_pathway_info`: Pathway metadata and descriptions
* `example_expr`: Example expression matrix (200 samples)
* `example_clinical`: Example clinical data with survival information

### Documentation

* Comprehensive function documentation
* Quick start vignette
* README with examples

## Data Sources

Gene sets were curated from:
- FerrDb V3 (http://www.zhounan.org/ferrdb/)
- MSigDB (https://www.gsea-msigdb.org/)
- KEGG (https://www.kegg.jp/)
- Primary literature (2012-2024)
