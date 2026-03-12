# =============================================================================
# TCGA Pan-Cancer Analysis using CellDeathAnalysis
# =============================================================================
#
# This script performs comprehensive cell death pathway analysis on TCGA data
# for the CellDeathAnalysis paper
#
# Author: Keran Sun
# Email: s1214844197@163.com
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Environment Setup
# -----------------------------------------------------------------------------

# Install required packages if not available
packages_cran <- c("tidyverse", "survival", "survminer", "ggpubr", 
                   "pheatmap", "RColorBrewer", "ggrepel", "patchwork",
                   "corrplot", "reshape2")

packages_bioc <- c("TCGAbiolinks", "SummarizedExperiment", "DESeq2",
                   "ComplexHeatmap", "circlize")

# Install CRAN packages
for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load packages
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(patchwork)
library(CellDeathAnalysis)

# Set working directory and create output folders
# setwd("your/working/directory")
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Download TCGA Pan-Cancer Data
# -----------------------------------------------------------------------------

# Define cancer types to analyze (33 TCGA cancer types)
cancer_types <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)

# For demonstration, use a subset (you can expand to all)
# demo_cancers <- c("TCGA-BRCA", "TCGA-LUAD", "TCGA-STAD", "TCGA-LIHC", "TCGA-COAD")
demo_cancers <- cancer_types  # Use all for publication

# Function to download and process TCGA data
download_tcga_data <- function(project, save_dir = "data/tcga") {
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  save_file <- file.path(save_dir, paste0(project, "_data.rds"))
  
  # Check if already downloaded
  if (file.exists(save_file)) {
    message("Loading cached data for ", project)
    return(readRDS(save_file))
  }
  
  message("Downloading ", project, "...")
  
  tryCatch({
    # Query gene expression data
    query <- GDCquery(
      project = project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    
    # Download
    GDCdownload(query, method = "api", files.per.chunk = 50)
    
    # Prepare data - use legacy = FALSE to avoid some issues
    data <- GDCprepare(query, summarizedExperiment = TRUE)
    
    # Extract expression matrix (TPM)
    if ("tpm_unstrand" %in% assayNames(data)) {
      expr <- assay(data, "tpm_unstrand")
    } else if ("unstranded" %in% assayNames(data)) {
      expr <- assay(data, "unstranded")
    } else {
      expr <- assay(data, 1)
    }
    
    # Get gene symbols
    gene_info <- rowData(data)
    if ("gene_name" %in% colnames(gene_info)) {
      rownames(expr) <- gene_info$gene_name
    }
    
    # Remove duplicates (keep highest mean expression)
    dup_genes <- duplicated(rownames(expr))
    expr <- expr[!dup_genes, ]
    
    # Get sample info from colData
    sample_info <- as.data.frame(colData(data))
    
    # Determine sample type
    if ("sample_type" %in% colnames(sample_info)) {
      sample_info$sample_type_simple <- ifelse(
        grepl("Normal", sample_info$sample_type, ignore.case = TRUE), 
        "Normal", "Tumor"
      )
    } else if ("shortLetterCode" %in% colnames(sample_info)) {
      sample_info$sample_type_simple <- ifelse(
        sample_info$shortLetterCode == "NT", "Normal", "Tumor"
      )
    } else {
      # Use barcode to determine (positions 14-15: 01-09 = tumor, 10-19 = normal)
      sample_code <- substr(colnames(expr), 14, 15)
      sample_info$sample_type_simple <- ifelse(
        as.numeric(sample_code) < 10, "Tumor", "Normal"
      )
    }
    
    # Get clinical data separately to avoid merge issues
    clinical <- tryCatch({
      GDCquery_clinic(project, type = "clinical")
    }, error = function(e) {
      message("  Warning: Could not fetch clinical data: ", e$message)
      NULL
    })
    
    # Create result object
    result <- list(
      expr = expr,
      sample_info = sample_info,
      clinical = clinical,
      project = project
    )
    
    # Save
    saveRDS(result, save_file)
    
    return(result)
    
  }, error = function(e) {
    message("Error downloading ", project, ": ", e$message)
    return(NULL)
  })
}

# Download data for each cancer type
tcga_data_list <- list()

for (cancer in demo_cancers) {
  result <- download_tcga_data(cancer)
  if (!is.null(result)) {
    tcga_data_list[[cancer]] <- result
    message("Successfully processed ", cancer)
  } else {
    message("Skipped ", cancer, " due to errors")
  }
}

# Remove NULL entries
tcga_data_list <- tcga_data_list[!sapply(tcga_data_list, is.null)]

if (length(tcga_data_list) == 0) {
  stop("No data successfully downloaded!")
}

message("\nSuccessfully downloaded ", length(tcga_data_list), " cancer types")

# Save combined data
saveRDS(tcga_data_list, "results/tcga_data_list.rds")

# -----------------------------------------------------------------------------
# 2. Calculate Cell Death Scores for All Cancers
# -----------------------------------------------------------------------------

# Function to calculate death scores for one cancer
calculate_cancer_scores <- function(tcga_data) {
  
  if (is.null(tcga_data)) return(NULL)
  
  # Log2 transform TPM
  expr <- log2(tcga_data$expr + 1)
  
  # Calculate scores
  scores <- calculate_death_score(
    expr,
    pathways = "all",
    method = "zscore",
    min_genes = 5,
    verbose = FALSE
  )
  
  # Add sample info
  scores$sample_id <- rownames(scores)
  
  # Match sample type
  if ("barcode" %in% colnames(tcga_data$sample_info)) {
    match_idx <- match(scores$sample_id, tcga_data$sample_info$barcode)
  } else {
    match_idx <- match(scores$sample_id, rownames(tcga_data$sample_info))
  }
  
  scores$sample_type <- tcga_data$sample_info$sample_type_simple[match_idx]
  scores$project <- tcga_data$project
  
  return(scores)
}

# Calculate scores for all cancers
all_scores <- list()

for (cancer in names(tcga_data_list)) {
  message("Calculating scores for ", cancer)
  result <- calculate_cancer_scores(tcga_data_list[[cancer]])
  if (!is.null(result)) {
    all_scores[[cancer]] <- result
  }
}

# Remove NULL entries
all_scores <- all_scores[!sapply(all_scores, is.null)]

# Combine all scores
pan_cancer_scores <- do.call(rbind, all_scores)
rownames(pan_cancer_scores) <- NULL

# Remove rows with NA sample_type
pan_cancer_scores <- pan_cancer_scores[!is.na(pan_cancer_scores$sample_type), ]

message("Total samples with scores: ", nrow(pan_cancer_scores))

# Save
saveRDS(pan_cancer_scores, "results/pan_cancer_scores.rds")
write.csv(pan_cancer_scores, "results/tables/pan_cancer_scores.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Pan-Cancer Heatmap (Figure 3A)
# -----------------------------------------------------------------------------

# Calculate mean scores per cancer type and sample type
pathway_cols <- names(get_death_geneset("all", type = "simple"))

mean_scores <- pan_cancer_scores %>%
  group_by(project, sample_type) %>%
  summarise(across(all_of(pathway_cols), mean, na.rm = TRUE), .groups = "drop")

# Create matrix for heatmap (Tumor samples)
tumor_matrix <- mean_scores %>%
  filter(sample_type == "Tumor") %>%
  select(-sample_type) %>%
  column_to_rownames("project") %>%
  as.matrix()

# Scale by pathway
tumor_matrix_scaled <- scale(tumor_matrix)

# Color palette
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

# Create annotation
cancer_abbr <- gsub("TCGA-", "", rownames(tumor_matrix_scaled))

# Plot heatmap
pdf("results/figures/Figure3A_pancancer_heatmap.pdf", width = 12, height = 10)

ht <- Heatmap(
  tumor_matrix_scaled,
  name = "Z-score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 45,
  row_labels = cancer_abbr,
  column_title = "Cell Death Pathways",
  row_title = "TCGA Cancer Types",
  heatmap_legend_param = list(
    title = "Z-score",
    legend_direction = "vertical"
  )
)

draw(ht)
dev.off()

# -----------------------------------------------------------------------------
# 4. Tumor vs Normal Comparison (Figure 3B)
# -----------------------------------------------------------------------------

# Function to compare tumor vs normal for each pathway
compare_tumor_normal <- function(scores_df, pathway) {
  
  tumor <- scores_df[[pathway]][scores_df$sample_type == "Tumor"]
  normal <- scores_df[[pathway]][scores_df$sample_type == "Normal"]
  
  if (length(normal) < 3) {
    return(data.frame(
      pathway = pathway,
      mean_tumor = mean(tumor, na.rm = TRUE),
      mean_normal = NA,
      log2FC = NA,
      p_value = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  test <- wilcox.test(tumor, normal)
  
  data.frame(
    pathway = pathway,
    mean_tumor = mean(tumor, na.rm = TRUE),
    mean_normal = mean(normal, na.rm = TRUE),
    log2FC = log2(mean(tumor, na.rm = TRUE) / mean(normal, na.rm = TRUE)),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

# Compare for each cancer
comparison_results <- list()

for (cancer in names(tcga_data_list)) {
  cancer_scores <- all_scores[[cancer]]
  
  if (sum(cancer_scores$sample_type == "Normal") >= 3) {
    results <- lapply(pathway_cols, function(p) {
      compare_tumor_normal(cancer_scores, p)
    })
    results_df <- do.call(rbind, results)
    results_df$cancer <- cancer
    comparison_results[[cancer]] <- results_df
  }
}

comparison_df <- do.call(rbind, comparison_results)
comparison_df$p_adjust <- p.adjust(comparison_df$p_value, method = "BH")
comparison_df$significance <- ifelse(comparison_df$p_adjust < 0.001, "***",
                                      ifelse(comparison_df$p_adjust < 0.01, "**",
                                             ifelse(comparison_df$p_adjust < 0.05, "*", "")))

# Save
write.csv(comparison_df, "results/tables/tumor_vs_normal_comparison.csv", row.names = FALSE)

# Plot heatmap of log2FC
fc_matrix <- comparison_df %>%
  select(cancer, pathway, log2FC) %>%
  pivot_wider(names_from = pathway, values_from = log2FC) %>%
  column_to_rownames("cancer") %>%
  as.matrix()

fc_matrix[is.na(fc_matrix)] <- 0
fc_matrix[is.infinite(fc_matrix)] <- 0

pdf("results/figures/Figure3B_tumor_normal_FC.pdf", width = 12, height = 8)

col_fun_fc <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(
  fc_matrix,
  name = "log2FC",
  col = col_fun_fc,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_labels = gsub("TCGA-", "", rownames(fc_matrix)),
  column_names_rot = 45,
  column_title = "Tumor vs Normal (log2 Fold Change)",
  row_title = "Cancer Type"
)

dev.off()

# -----------------------------------------------------------------------------
# 5. Boxplot for Key Pathways (Figure 3C)
# -----------------------------------------------------------------------------

# Select key pathways
key_pathways <- c("ferroptosis", "cuproptosis", "pyroptosis", "apoptosis")

# Prepare data
boxplot_data <- pan_cancer_scores %>%
  select(sample_type, project, all_of(key_pathways)) %>%
  pivot_longer(cols = all_of(key_pathways), 
               names_to = "pathway", 
               values_to = "score") %>%
  mutate(
    cancer = gsub("TCGA-", "", project),
    pathway = tools::toTitleCase(pathway)
  )

# Plot
pdf("results/figures/Figure3C_pathway_boxplots.pdf", width = 14, height = 10)

p_box <- ggplot(boxplot_data, aes(x = cancer, y = score, fill = sample_type)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~pathway, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Normal" = "#4DAF4A", "Tumor" = "#E41A1C")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(
    x = "Cancer Type",
    y = "Death Score (Z-score)",
    fill = "Sample Type"
  )

print(p_box)
dev.off()

# -----------------------------------------------------------------------------
# 6. Survival Analysis (Figure 4)
# -----------------------------------------------------------------------------

# Function to perform survival analysis for one cancer
perform_survival_analysis <- function(tcga_data, scores) {
  
  if (is.null(tcga_data) || is.null(tcga_data$clinical)) {
    return(NULL)
  }
  
  # Get clinical data
  clinical <- tcga_data$clinical
  
  # Check required columns
  if (!all(c("submitter_id") %in% colnames(clinical))) {
    message("  Missing submitter_id column")
    return(NULL)
  }
  
  # Match samples
  scores_tumor <- scores %>%
    filter(sample_type == "Tumor")
  
  if (nrow(scores_tumor) < 20) {
    message("  Too few tumor samples")
    return(NULL)
  }
  
  # Get patient IDs (first 12 characters of barcode)
  scores_tumor$patient_id <- substr(scores_tumor$sample_id, 1, 12)
  
  # Prepare survival columns
  clinical$OS_time <- NA
  clinical$OS_status <- NA
  
  # Try different column names for survival data
  if ("days_to_death" %in% colnames(clinical)) {
    clinical$OS_time <- ifelse(
      !is.na(clinical$days_to_death),
      as.numeric(clinical$days_to_death),
      as.numeric(clinical$days_to_last_follow_up)
    )
  } else if ("days_to_last_follow_up" %in% colnames(clinical)) {
    clinical$OS_time <- as.numeric(clinical$days_to_last_follow_up)
  }
  
  if ("vital_status" %in% colnames(clinical)) {
    clinical$OS_status <- ifelse(
      tolower(clinical$vital_status) %in% c("dead", "deceased"), 1, 0
    )
  }
  
  # Match with clinical
  matched <- merge(
    scores_tumor,
    clinical[, c("submitter_id", "OS_time", "OS_status")],
    by.x = "patient_id",
    by.y = "submitter_id",
    all.x = FALSE
  )
  
  # Remove NA and invalid values
  matched <- matched %>%
    filter(!is.na(OS_time) & !is.na(OS_status) & OS_time > 0)
  
  if (nrow(matched) < 30) {
    message("  Too few samples with survival data: ", nrow(matched))
    return(NULL)
  }
  
  return(matched)
}

# Perform survival analysis for each cancer
survival_results <- list()

for (cancer in names(tcga_data_list)) {
  message("Survival analysis for ", cancer)
  
  if (!cancer %in% names(all_scores)) {
    message("  No scores available, skipping")
    next
  }
  
  matched <- perform_survival_analysis(
    tcga_data_list[[cancer]],
    all_scores[[cancer]]
  )
  
  if (is.null(matched)) {
    message("  No valid survival data")
    next
  }
  
  message("  Samples with survival data: ", nrow(matched))
  
  # Analyze each pathway
  for (pathway in pathway_cols) {
    
    if (!pathway %in% colnames(matched)) next
    if (all(is.na(matched[[pathway]]))) next
    
    tryCatch({
      # Classify by median
      pathway_values <- matched[[pathway]]
      med_val <- median(pathway_values, na.rm = TRUE)
      
      matched$group <- ifelse(pathway_values > med_val, "High", "Low")
      matched$group <- factor(matched$group, levels = c("Low", "High"))
      
      # Skip if not enough in each group
      group_table <- table(matched$group)
      if (any(group_table < 10)) next
      
      # Cox regression
      cox_data <- data.frame(
        time = matched$OS_time,
        status = matched$OS_status,
        score = pathway_values
      )
      cox_data <- cox_data[complete.cases(cox_data), ]
      
      cox_fit <- coxph(Surv(time, status) ~ score, data = cox_data)
      cox_sum <- summary(cox_fit)
      
      # Log-rank test
      surv_data <- data.frame(
        time = matched$OS_time,
        status = matched$OS_status,
        group = matched$group
      )
      surv_data <- surv_data[complete.cases(surv_data), ]
      
      surv_diff <- survdiff(Surv(time, status) ~ group, data = surv_data)
      log_rank_p <- 1 - pchisq(surv_diff$chisq, df = 1)
      
      survival_results[[paste(cancer, pathway, sep = "_")]] <- data.frame(
        cancer = cancer,
        pathway = pathway,
        n = nrow(cox_data),
        hr = cox_sum$conf.int[1, 1],
        hr_lower = cox_sum$conf.int[1, 3],
        hr_upper = cox_sum$conf.int[1, 4],
        cox_p = cox_sum$coefficients[1, 5],
        log_rank_p = log_rank_p,
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      # Silently skip errors for individual pathways
    })
  }
}

# Check if we have results
if (length(survival_results) == 0) {
  message("\nWarning: No survival analysis results generated")
  survival_df <- data.frame()
} else {
  survival_df <- do.call(rbind, survival_results)
  rownames(survival_df) <- NULL
  
  # Adjust p-values
  survival_df$cox_p_adj <- p.adjust(survival_df$cox_p, method = "BH")
  survival_df$log_rank_p_adj <- p.adjust(survival_df$log_rank_p, method = "BH")
  
  # Add significance
  survival_df$significance <- ifelse(survival_df$cox_p_adj < 0.001, "***",
                                      ifelse(survival_df$cox_p_adj < 0.01, "**",
                                             ifelse(survival_df$cox_p_adj < 0.05, "*", "")))
  
  message("\nSurvival analysis completed for ", 
          length(unique(survival_df$cancer)), " cancer types")
}

# Save
write.csv(survival_df, "results/tables/survival_analysis_results.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Forest Plot (Figure 4A)
# -----------------------------------------------------------------------------

if (nrow(survival_df) > 0 && "ferroptosis" %in% survival_df$pathway) {
  
  # Select ferroptosis results
  ferro_survival <- survival_df %>%
    filter(pathway == "ferroptosis") %>%
    arrange(hr)
  
  if (nrow(ferro_survival) > 0) {
    ferro_survival$cancer_abbr <- gsub("TCGA-", "", ferro_survival$cancer)
    ferro_survival$cancer_abbr <- factor(ferro_survival$cancer_abbr, 
                                          levels = ferro_survival$cancer_abbr)
    
    pdf("results/figures/Figure4A_ferroptosis_forest.pdf", width = 10, height = 12)
    
    p_forest <- ggplot(ferro_survival, aes(x = hr, y = cancer_abbr)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      geom_errorbarh(aes(xmin = hr_lower, xmax = hr_upper), height = 0.2) +
      geom_point(aes(color = cox_p_adj < 0.05), size = 3) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                         labels = c("TRUE" = "p < 0.05", "FALSE" = "NS"),
                         name = "Significance") +
      scale_x_log10() +
      theme_bw() +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 10)
      ) +
      labs(
        x = "Hazard Ratio (95% CI)",
        y = NULL,
        title = "Ferroptosis Score and Overall Survival"
      )
    
    print(p_forest)
    dev.off()
    
    message("Forest plot saved")
  } else {
    message("No ferroptosis survival data for forest plot")
  }
} else {
  message("Skipping forest plot - no survival results")
}

# -----------------------------------------------------------------------------
# 8. Kaplan-Meier Curves for Selected Cancers (Figure 4B)
# -----------------------------------------------------------------------------

if (nrow(survival_df) > 0) {
  
  # Select cancers with significant ferroptosis association
  ferro_results <- survival_df %>% filter(pathway == "ferroptosis")
  
  if (nrow(ferro_results) > 0) {
    sig_cancers <- ferro_results %>%
      filter(cox_p_adj < 0.05) %>%
      arrange(hr) %>%
      head(4) %>%
      pull(cancer)
    
    # If no significant, use top 4 by p-value
    if (length(sig_cancers) == 0) {
      sig_cancers <- ferro_results %>%
        arrange(cox_p) %>%
        head(4) %>%
        pull(cancer)
    }
    
    # Create KM plots
    km_plots <- list()
    
    for (cancer in sig_cancers) {
      
      if (!cancer %in% names(tcga_data_list) || !cancer %in% names(all_scores)) next
      
      matched <- perform_survival_analysis(
        tcga_data_list[[cancer]],
        all_scores[[cancer]]
      )
      
      if (is.null(matched) || nrow(matched) < 20) next
      
      tryCatch({
        matched$group <- ifelse(
          matched$ferroptosis > median(matched$ferroptosis, na.rm = TRUE),
          "High", "Low"
        )
        matched$group <- factor(matched$group, levels = c("Low", "High"))
        
        # Convert days to months
        matched$OS_months <- matched$OS_time / 30
        
        fit <- survfit(Surv(OS_months, OS_status) ~ group, data = matched)
        
        km_plots[[cancer]] <- ggsurvplot(
          fit,
          data = matched,
          pval = TRUE,
          risk.table = TRUE,
          palette = c("#2E9FDF", "#E7B800"),
          title = gsub("TCGA-", "", cancer),
          xlab = "Time (months)",
          legend.title = "Ferroptosis",
          legend.labs = c("Low", "High"),
          risk.table.height = 0.25,
          ggtheme = theme_bw()
        )
        
      }, error = function(e) {
        message("  Error creating KM plot for ", cancer, ": ", e$message)
      })
    }
    
    # Combine plots if we have any
    if (length(km_plots) > 0) {
      pdf("results/figures/Figure4B_KM_curves.pdf", width = 12, height = 12)
      tryCatch({
        arrange_ggsurvplots(km_plots, ncol = 2, nrow = 2)
      }, error = function(e) {
        # If arrange fails, plot individually
        for (i in seq_along(km_plots)) {
          print(km_plots[[i]])
        }
      })
      dev.off()
      message("KM curves saved")
    } else {
      message("No KM plots generated")
    }
  }
} else {
  message("Skipping KM plots - no survival results")
}

# -----------------------------------------------------------------------------
# 9. Pathway Correlation Analysis (Figure 5)
# -----------------------------------------------------------------------------

# Calculate correlation between pathways
cor_matrix <- cor(
  pan_cancer_scores[, pathway_cols],
  use = "pairwise.complete.obs",
  method = "spearman"
)

pdf("results/figures/Figure5_pathway_correlation.pdf", width = 10, height = 10)

col_fun_cor <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(
  cor_matrix,
  name = "Correlation",
  col = col_fun_cor,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y, gp = gpar(fontsize = 8))
  },
  column_title = "Cell Death Pathway Correlation (Spearman)",
  heatmap_legend_param = list(
    title = "Spearman\nCorrelation"
  )
)

dev.off()

# -----------------------------------------------------------------------------
# 10. Summary Statistics (Table 1)
# -----------------------------------------------------------------------------

# Gene set statistics
geneset_stats <- data.frame(
  pathway = pathway_cols,
  n_genes = sapply(pathway_cols, function(p) {
    length(get_death_geneset(p, type = "all"))
  }),
  stringsAsFactors = FALSE
)

# Add pathway info
pathway_info <- list_death_pathways(detailed = TRUE)
geneset_stats <- merge(geneset_stats, pathway_info, by = "pathway")

write.csv(geneset_stats, "results/tables/Table1_geneset_statistics.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 11. Generate Summary Report
# -----------------------------------------------------------------------------

# Summary statistics
summary_stats <- list(
  n_cancers = length(unique(pan_cancer_scores$project)),
  n_samples_total = nrow(pan_cancer_scores),
  n_tumor = sum(pan_cancer_scores$sample_type == "Tumor", na.rm = TRUE),
  n_normal = sum(pan_cancer_scores$sample_type == "Normal", na.rm = TRUE),
  n_pathways = length(pathway_cols),
  n_genes_total = length(get_all_death_genes()),
  n_significant_survival = if(nrow(survival_df) > 0) sum(survival_df$cox_p_adj < 0.05, na.rm = TRUE) else 0
)

# Print summary
cat("\n")
cat("=============================================================\n")
cat("           TCGA Pan-Cancer Analysis Summary                  \n")
cat("=============================================================\n")
cat("\n")
cat("Data Summary:\n")
cat("  Cancer types analyzed:", summary_stats$n_cancers, "\n")
cat("  Total samples:", summary_stats$n_samples_total, "\n")
cat("  Tumor samples:", summary_stats$n_tumor, "\n")
cat("  Normal samples:", summary_stats$n_normal, "\n")
cat("\n")
cat("Gene Sets:\n")
cat("  Cell death pathways:", summary_stats$n_pathways, "\n")
cat("  Total genes:", summary_stats$n_genes_total, "\n")
cat("\n")
cat("Key Findings:\n")
cat("  Significant survival associations:", summary_stats$n_significant_survival, "\n")
cat("\n")
cat("Output Files:\n")
cat("  - results/tables/pan_cancer_scores.csv\n")
cat("  - results/tables/tumor_vs_normal_comparison.csv\n")
cat("  - results/tables/survival_analysis_results.csv\n")
cat("  - results/figures/Figure3A_pancancer_heatmap.pdf\n")
cat("  - results/figures/Figure3B_tumor_normal_FC.pdf\n")
cat("  - results/figures/Figure3C_pathway_boxplots.pdf\n")
cat("  - results/figures/Figure4A_ferroptosis_forest.pdf\n")
cat("  - results/figures/Figure4B_KM_curves.pdf\n")
cat("  - results/figures/Figure5_pathway_correlation.pdf\n")
cat("\n")
cat("=============================================================\n")

# Save summary
saveRDS(summary_stats, "results/summary_statistics.rds")

cat("\nAnalysis completed successfully!\n")
