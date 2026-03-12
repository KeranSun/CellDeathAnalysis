# =============================================================================
# GEO Data Analysis using CellDeathAnalysis (Simplified Alternative)
# =============================================================================
#
# This script provides a simpler analysis using GEO data
# Can be run more quickly than full TCGA analysis
#
# Author: Keran Sun
# =============================================================================

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

# Install required packages
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

library(GEOquery)
library(tidyverse)
library(survival)
library(survminer)
library(CellDeathAnalysis)

# Create output directory
dir.create("results_geo", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Example 1: Liver Cancer (GSE14520)
# -----------------------------------------------------------------------------

cat("=== Downloading GSE14520 (Liver Cancer) ===\n")

# Download GEO data
gse14520 <- getGEO("GSE14520", GSEMatrix = TRUE)[[1]]

# Get expression matrix
expr <- exprs(gse14520)

# Get phenotype data
pdata <- pData(gse14520)

# Check available columns
cat("Available clinical variables:\n")
print(colnames(pdata))

# Get gene symbols from annotation
fdata <- fData(gse14520)
gene_symbols <- fdata$`Gene Symbol`

# Handle probe-to-gene mapping (keep probe with highest mean expression for each gene)
expr_df <- as.data.frame(expr)
expr_df$gene <- gene_symbols

# Remove probes without gene symbols
expr_df <- expr_df[expr_df$gene != "" & !is.na(expr_df$gene), ]

# For duplicated genes, keep the one with highest mean expression
expr_df$mean_expr <- rowMeans(expr_df[, -ncol(expr_df)])
expr_df <- expr_df %>%
  group_by(gene) %>%
  slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

# Create final expression matrix
rownames_new <- expr_df$gene
expr_matrix <- as.matrix(expr_df[, 1:(ncol(expr_df)-2)])
rownames(expr_matrix) <- rownames_new

cat("Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")

# -----------------------------------------------------------------------------
# Calculate Cell Death Scores
# -----------------------------------------------------------------------------

cat("\n=== Calculating Cell Death Scores ===\n")

scores <- calculate_death_score(
  expr_matrix,
  pathways = "all",
  method = "zscore",
  verbose = TRUE
)

cat("Scores calculated for", ncol(scores), "pathways\n")

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------

cat("\n=== Creating Visualizations ===\n")

# Get sample groups (tumor vs normal if available)
if ("source_name_ch1" %in% colnames(pdata)) {
  group <- pdata$source_name_ch1
} else if ("tissue:ch1" %in% colnames(pdata)) {
  group <- pdata$`tissue:ch1`
} else {
  group <- rep("Sample", nrow(scores))
}

# Boxplot
pdf("results_geo/boxplot_gse14520.pdf", width = 12, height = 8)
p1 <- plot_death_boxplot(
  scores,
  group = group,
  pathways = c("ferroptosis", "pyroptosis", "cuproptosis", "apoptosis")
)
print(p1)
dev.off()

# Radar chart
pdf("results_geo/radar_gse14520.pdf", width = 10, height = 10)
p2 <- plot_death_radar(scores, group = group)
print(p2)
dev.off()

# Correlation heatmap
pdf("results_geo/correlation_gse14520.pdf", width = 10, height = 10)
tryCatch({
  p3 <- plot_pathway_correlation(scores)
  print(p3)
}, error = function(e) {
  # Fallback if ComplexHeatmap not available
  cor_mat <- cor(scores, use = "pairwise.complete.obs")
  heatmap(cor_mat, main = "Pathway Correlation")
})
dev.off()

cat("Visualizations saved to results_geo/\n")

# -----------------------------------------------------------------------------
# Survival Analysis (if survival data available)
# -----------------------------------------------------------------------------

cat("\n=== Survival Analysis ===\n")

# Check for survival data
has_survival <- FALSE

# Common survival column names in GEO
surv_time_cols <- c("os:ch1", "survival_months:ch1", "follow_up_time:ch1", 
                    "OS_MONTHS", "survival_time")
surv_status_cols <- c("status:ch1", "event:ch1", "OS_STATUS", "vital_status")

for (tc in surv_time_cols) {
  if (tc %in% colnames(pdata)) {
    time <- as.numeric(pdata[[tc]])
    for (sc in surv_status_cols) {
      if (sc %in% colnames(pdata)) {
        status <- pdata[[sc]]
        # Convert to 0/1 if needed
        if (is.character(status)) {
          status <- ifelse(tolower(status) %in% c("dead", "1", "yes", "deceased"), 1, 0)
        }
        has_survival <- TRUE
        break
      }
    }
    break
  }
}

if (has_survival && sum(!is.na(time)) > 20) {
  cat("Survival data found! Performing analysis...\n")
  
  # Ferroptosis survival analysis
  surv_result <- death_survival(
    scores,
    time = time,
    status = status,
    pathway = "ferroptosis",
    method = "median"
  )
  
  print(surv_result)
  
  # KM plot
  pdf("results_geo/survival_ferroptosis.pdf", width = 8, height = 6)
  tryCatch({
    p_surv <- plot_survival(surv_result)
    print(p_surv)
  }, error = function(e) {
    cat("survminer not available for KM plot\n")
  })
  dev.off()
  
  # Batch survival analysis
  batch_result <- batch_survival(scores, time, status)
  write.csv(batch_result, "results_geo/survival_all_pathways.csv", row.names = FALSE)
  print(batch_result)
  
} else {
  cat("No survival data available in this dataset\n")
}

# -----------------------------------------------------------------------------
# Enrichment Analysis Example
# -----------------------------------------------------------------------------

cat("\n=== Enrichment Analysis Example ===\n")

# Simulate differentially expressed genes (top variable genes)
gene_vars <- apply(expr_matrix, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:200]

# Intersect with death genes to create a realistic gene list
all_death_genes <- get_all_death_genes()
example_genes <- intersect(top_genes, all_death_genes)

if (length(example_genes) > 5) {
  ora_result <- death_enrich_ora(example_genes, min_genes = 2)
  
  if (!is.null(ora_result)) {
    print(ora_result)
    
    pdf("results_geo/enrichment_bar.pdf", width = 8, height = 6)
    p_enrich <- plot_enrichment_bar(ora_result, show_all = TRUE)
    print(p_enrich)
    dev.off()
  }
}

# -----------------------------------------------------------------------------
# Machine Learning Example
# -----------------------------------------------------------------------------

cat("\n=== Machine Learning Example ===\n")

# Create binary outcome (high vs low ferroptosis for demo)
outcome <- factor(ifelse(scores$ferroptosis > median(scores$ferroptosis), "High", "Low"))

# Build model
tryCatch({
  model <- death_build_model(
    scores,
    outcome = outcome,
    method = "rf",
    train_ratio = 0.7
  )
  
  print(model)
  
  # Plot importance
  pdf("results_geo/feature_importance.pdf", width = 8, height = 6)
  p_imp <- plot_model(model, type = "importance")
  print(p_imp)
  dev.off()
  
}, error = function(e) {
  cat("ML analysis failed:", e$message, "\n")
  cat("Make sure randomForest package is installed\n")
})

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

cat("\n=== Saving Results ===\n")

# Save scores
write.csv(scores, "results_geo/death_scores_gse14520.csv")

# Save summary
sink("results_geo/analysis_summary.txt")
cat("CellDeathAnalysis - GEO Data Analysis Summary\n")
cat("==============================================\n\n")
cat("Dataset: GSE14520 (Liver Cancer)\n")
cat("Samples:", ncol(expr_matrix), "\n")
cat("Genes:", nrow(expr_matrix), "\n")
cat("Pathways analyzed:", ncol(scores), "\n\n")
cat("Score ranges:\n")
print(summary(scores[, 1:5]))
sink()

cat("\n")
cat("=============================================================\n")
cat("                    Analysis Complete!                       \n")
cat("=============================================================\n")
cat("\nOutput files:\n")
cat("  - results_geo/death_scores_gse14520.csv\n")
cat("  - results_geo/boxplot_gse14520.pdf\n")
cat("  - results_geo/radar_gse14520.pdf\n")
cat("  - results_geo/correlation_gse14520.pdf\n")
if (has_survival) {
  cat("  - results_geo/survival_ferroptosis.pdf\n")
  cat("  - results_geo/survival_all_pathways.csv\n")
}
cat("  - results_geo/enrichment_bar.pdf\n")
cat("  - results_geo/feature_importance.pdf\n")
cat("  - results_geo/analysis_summary.txt\n")
cat("\n")

# -----------------------------------------------------------------------------
# Example 2: Quick analysis template for any GEO dataset
# -----------------------------------------------------------------------------

# Function for quick GEO analysis
quick_geo_analysis <- function(gse_id, output_dir = NULL) {
  
  if (is.null(output_dir)) {
    output_dir <- paste0("results_", gse_id)
  }
  dir.create(output_dir, showWarnings = FALSE)
  
  cat("Downloading", gse_id, "...\n")
  gse <- getGEO(gse_id, GSEMatrix = TRUE)[[1]]
  
  expr <- exprs(gse)
  fdata <- fData(gse)
  
  # Try to find gene symbol column
  symbol_cols <- c("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_symbol")
  gene_col <- NULL
  for (col in symbol_cols) {
    if (col %in% colnames(fdata)) {
      gene_col <- col
      break
    }
  }
  
  if (is.null(gene_col)) {
    stop("Cannot find gene symbol column in annotation")
  }
  
  # Prepare expression matrix
  genes <- fdata[[gene_col]]
  expr_df <- as.data.frame(expr)
  expr_df$gene <- genes
  expr_df <- expr_df[expr_df$gene != "" & !is.na(expr_df$gene), ]
  expr_df$mean_expr <- rowMeans(expr_df[, -ncol(expr_df)])
  expr_df <- expr_df %>%
    group_by(gene) %>%
    slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  expr_matrix <- as.matrix(expr_df[, 1:(ncol(expr_df)-2)])
  rownames(expr_matrix) <- expr_df$gene
  
  # Calculate scores
  scores <- calculate_death_score(expr_matrix, method = "zscore")
  
  # Save
  write.csv(scores, file.path(output_dir, "death_scores.csv"))
  
  # Visualize
  pdf(file.path(output_dir, "analysis_plots.pdf"), width = 12, height = 10)
  
  # Boxplot
  print(plot_death_boxplot(scores, 
                           pathways = c("ferroptosis", "pyroptosis", "apoptosis", "autophagy")))
  
  # Radar
  print(plot_death_radar(scores))
  
  dev.off()
  
  cat("Analysis complete! Results saved to", output_dir, "\n")
  
  return(list(
    expression = expr_matrix,
    scores = scores,
    phenotype = pData(gse)
  ))
}

# Usage example:
# result <- quick_geo_analysis("GSE62944")
