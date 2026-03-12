# =============================================================================
# CellDeathAnalysis Complete Example Analysis v0.2.0
# =============================================================================
#
# This script demonstrates all main functions of the CellDeathAnalysis package
# including scoring, visualization, survival analysis, and enrichment analysis.
#
# Author: Keran Sun
# Email: s1214844197@163.com
# =============================================================================

# Load package
library(CellDeathAnalysis)

# =============================================================================
# 1. Load Example Data
# =============================================================================

cat("========== 1. Load Data ==========\n")

# Load example expression data
data(example_expr)
cat("Expression matrix:", dim(example_expr)[1], "genes x", 
    dim(example_expr)[2], "samples\n")

# Load example clinical data
data(example_clinical)
cat("Clinical data:", nrow(example_clinical), "samples x", 
    ncol(example_clinical), "variables\n")

# View data structure
cat("\nFirst 5 rows and columns of expression matrix:\n")
print(example_expr[1:5, 1:5])

cat("\nClinical data structure:\n")
print(head(example_clinical, 5))

cat("\nSample distribution:\n")
print(table(example_clinical$group))

# =============================================================================
# 2. Explore Gene Sets
# =============================================================================

cat("\n========== 2. Explore Gene Sets ==========\n")

# List all available pathways
cat("\nAvailable cell death types:\n")
print(list_death_pathways())

# View detailed information
cat("\nDetailed pathway information:\n")
info <- list_death_pathways(detailed = TRUE)
print(info[, c("pathway", "chinese_name", "total_genes", "core_genes")])

# Get specific pathway info
cat("\nFerroptosis information:\n")
get_pathway_info("ferroptosis")

# Get gene sets
cat("\nGet ferroptosis genes:\n")
ferro_genes <- get_death_geneset("ferroptosis", type = "all")
cat("Number of ferroptosis genes:", length(ferro_genes), "\n")
cat("First 10 genes:", paste(head(ferro_genes, 10), collapse = ", "), "\n")

# Check pathway overlap
cat("\nPathway gene overlap:\n")
overlap <- check_pathway_overlap(c("ferroptosis", "apoptosis", "autophagy"))
print(overlap)

# =============================================================================
# 3. Calculate Pathway Scores
# =============================================================================

cat("\n========== 3. Calculate Pathway Scores ==========\n")

# Method 1: Z-score (fastest)
cat("\nCalculating scores using Z-score method...\n")
scores <- calculate_death_score(
  example_expr,
  pathways = "all",
  method = "zscore",
  verbose = TRUE
)

cat("\nScore matrix:", dim(scores)[1], "samples x", 
    dim(scores)[2], "pathways\n")
cat("\nScore summary (first 4 pathways):\n")
print(summary(scores[, 1:4]))

# =============================================================================
# 4. Visualization
# =============================================================================

cat("\n========== 4. Visualization ==========\n")

# Prepare group info
group <- example_clinical$group

# 4.1 Boxplot comparison
cat("\nCreating boxplot...\n")
p1 <- plot_death_boxplot(
  scores, 
  group = group,
  pathways = c("ferroptosis", "pyroptosis", "cuproptosis", "apoptosis"),
  ncol = 2
)
print(p1)

# 4.2 Radar chart
cat("\nCreating radar chart...\n")
p2 <- plot_death_radar(scores, group = group)
print(p2)

# 4.3 Score distribution
cat("\nCreating distribution plot...\n")
p3 <- plot_score_distribution(
  scores, 
  pathways = c("ferroptosis", "pyroptosis"),
  type = "density",
  group = group
)
print(p3)

# 4.4 Heatmap (requires ComplexHeatmap or pheatmap)
cat("\nCreating heatmap...\n")
tryCatch({
  p4 <- plot_death_heatmap(scores)
  cat("Heatmap created successfully!\n")
}, error = function(e) {
  cat("Heatmap requires ComplexHeatmap or pheatmap package\n")
})

# 4.5 Pathway correlation
cat("\nCreating correlation plot...\n")
tryCatch({
  p5 <- plot_pathway_correlation(scores[, 1:6], method = "spearman")
  cat("Correlation plot created successfully!\n")
}, error = function(e) {
  cat("Correlation heatmap requires ComplexHeatmap package\n")
})

# =============================================================================
# 5. Sample Classification
# =============================================================================

cat("\n========== 5. Sample Classification ==========\n")

# Classify by ferroptosis score
ferro_groups <- classify_by_score(scores$ferroptosis, method = "median")
cat("\nFerroptosis score groups (median method):\n")
print(table(ferro_groups))

# Classify by quantile
ferro_groups_q <- classify_by_score(scores$ferroptosis, 
                                     method = "quantile", quantile = 0.75)
cat("\nFerroptosis score groups (75th percentile):\n")
print(table(ferro_groups_q))

# Calculate composite death index
death_index <- calculate_death_index(scores, method = "pca")
cat("\nDeath index summary:\n")
print(summary(death_index))

# =============================================================================
# 6. Survival Analysis (NEW in v0.2.0)
# =============================================================================

cat("\n========== 6. Survival Analysis ==========\n")

if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)
  
  # Analyze only tumor samples
  tumor_idx <- example_clinical$group == "Tumor"
  tumor_scores <- scores[tumor_idx, ]
  tumor_clinical <- example_clinical[tumor_idx, ]
  
  # 6.1 Single pathway survival analysis
  cat("\nSingle pathway survival analysis (ferroptosis):\n")
  surv_result <- death_survival(
    tumor_scores,
    time = tumor_clinical$OS_time,
    status = tumor_clinical$OS_status,
    pathway = "ferroptosis",
    method = "median"
  )
  print(surv_result)
  
  # 6.2 Plot survival curves
  if (requireNamespace("survminer", quietly = TRUE)) {
    cat("\nCreating survival curves...\n")
    p_surv <- plot_survival(surv_result)
    print(p_surv)
  }
  
  # 6.3 Batch survival analysis
  cat("\nBatch survival analysis for all pathways:\n")
  batch_result <- batch_survival(
    tumor_scores,
    time = tumor_clinical$OS_time,
    status = tumor_clinical$OS_status,
    method = "median"
  )
  print(batch_result)
  
  # 6.4 Forest plot
  cat("\nCreating forest plot...\n")
  p_forest <- plot_forest(batch_result)
  print(p_forest)
  
  # 6.5 Multivariate Cox regression
  cat("\nMultivariate Cox regression:\n")
  clinical_vars <- data.frame(
    age = tumor_clinical$age,
    stage = tumor_clinical$stage
  )
  # Remove samples with NA stage
  complete_idx <- complete.cases(clinical_vars)
  
  if (sum(complete_idx) > 20) {
    multi_cox <- death_cox_multivariate(
      tumor_scores[complete_idx, ],
      time = tumor_clinical$OS_time[complete_idx],
      status = tumor_clinical$OS_status[complete_idx],
      clinical_vars = clinical_vars[complete_idx, ],
      pathway = "ferroptosis"
    )
    print(multi_cox)
  }
  
  # 6.6 Time-dependent ROC
  if (requireNamespace("timeROC", quietly = TRUE)) {
    cat("\nTime-dependent ROC analysis:\n")
    roc_result <- death_time_roc(
      scores = tumor_scores$ferroptosis,
      time = tumor_clinical$OS_time,
      status = tumor_clinical$OS_status,
      times = c(12, 24, 36),
      pathway = "Ferroptosis"
    )
    cat("AUC at 12 months:", round(roc_result$auc[1], 3), "\n")
    cat("AUC at 24 months:", round(roc_result$auc[2], 3), "\n")
    cat("AUC at 36 months:", round(roc_result$auc[3], 3), "\n")
    
    # Plot ROC curves
    p_roc <- plot_time_roc(roc_result)
    print(p_roc)
  }
  
} else {
  cat("Package 'survival' is required for survival analysis\n")
  cat("Install with: install.packages('survival')\n")
}

# =============================================================================
# 7. Enrichment Analysis (NEW in v0.2.0)
# =============================================================================

cat("\n========== 7. Enrichment Analysis ==========\n")

# 7.1 Over-representation analysis (ORA)
cat("\nOver-representation analysis:\n")

# Create example gene list (e.g., differentially expressed genes)
example_degs <- c(
  # Ferroptosis genes
  "GPX4", "SLC7A11", "ACSL4", "TFRC", "NCOA4", "PTGS2",
  # Pyroptosis genes  
  "NLRP3", "CASP1", "GSDMD", "IL1B", "IL18", "PYCARD",
  # Cuproptosis genes
  "FDX1", "LIAS", "DLAT",
  # Some random genes
  "TP53", "MYC", "EGFR", "KRAS"
)

ora_result <- death_enrich_ora(
  gene_list = example_degs,
  pathways = "all",
  min_genes = 2
)
print(ora_result)

# Plot enrichment results
if (!is.null(ora_result) && nrow(ora_result) > 0) {
  cat("\nCreating enrichment bar plot...\n")
  p_enrich <- plot_enrichment_bar(ora_result, top_n = 8)
  print(p_enrich)
  
  cat("\nCreating enrichment dot plot...\n")
  p_dot <- plot_enrichment_dot(ora_result, top_n = 8)
  print(p_dot)
}

# 7.2 Compare groups
cat("\nComparing death pathways between Normal and Tumor:\n")
compare_result <- death_compare_groups(
  example_expr,
  group = example_clinical$group,
  pathways = "all",
  method = "wilcox"
)
print(compare_result)

# Volcano plot
cat("\nCreating volcano plot...\n")
p_volcano <- plot_compare_volcano(compare_result, 
                                   fc_cutoff = 0.3, 
                                   p_cutoff = 0.05)
print(p_volcano)

# 7.3 GSEA (requires fgsea package)
if (requireNamespace("fgsea", quietly = TRUE)) {
  cat("\nGene Set Enrichment Analysis:\n")
  
  # Create ranked gene list (simulated log2FC)
  set.seed(123)
  ranked_genes <- setNames(
    rnorm(nrow(example_expr)),
    rownames(example_expr)
  )
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  gsea_result <- death_enrich_gsea(
    ranked_genes,
    pathways = "all",
    nperm = 1000,
    min_size = 5
  )
  print(gsea_result)
  
  # Plot GSEA curve for top pathway
  if (nrow(gsea_result) > 0) {
    top_pathway <- gsea_result$pathway[1]
    cat("\nCreating GSEA curve for", top_pathway, "...\n")
    p_gsea <- plot_gsea_curve(ranked_genes, top_pathway)
    print(p_gsea)
  }
} else {
  cat("Package 'fgsea' is required for GSEA\n")
  cat("Install with: BiocManager::install('fgsea')\n")
}

# =============================================================================
# 8. Export Results
# =============================================================================

cat("\n========== 8. Export Results ==========\n")

# Combine scores with clinical data
output_data <- cbind(
  sample_id = rownames(scores),
  group = example_clinical$group,
  scores
)

cat("\nOutput data preview:\n")
print(head(output_data[, 1:5]))

cat("\nTo save results, use:\n")
cat("  write.csv(output_data, 'death_scores.csv', row.names = FALSE)\n")

# =============================================================================
# 9. Summary
# =============================================================================

cat("\n")
cat("=============================================================\n")
cat("                   Analysis Complete!                        \n")
cat("=============================================================\n")
cat("\n")
cat("Summary:\n")
cat("  - Analyzed", ncol(scores), "cell death pathways\n")
cat("  - Compared", sum(group == "Normal"), "normal vs", 
    sum(group == "Tumor"), "tumor samples\n")
cat("  - Ferroptosis High/Low:", table(ferro_groups)["High"], "vs", 
    table(ferro_groups)["Low"], "\n")
cat("\n")
cat("New in v0.2.0:\n")
cat("  - Survival analysis: death_survival(), batch_survival()\n")
cat("  - Time-dependent ROC: death_time_roc()\n")
cat("  - Enrichment analysis: death_enrich_ora(), death_enrich_gsea()\n")
cat("  - Group comparison: death_compare_groups()\n")
cat("\n")
cat("Next steps:\n")
cat("  - Use ssGSEA or GSVA for more accurate scores\n")
cat("  - Perform differential expression analysis\n")
cat("  - Integrate with other clinical features\n")
cat("  - Build predictive models\n")
cat("\n")
cat("=============================================================\n")
