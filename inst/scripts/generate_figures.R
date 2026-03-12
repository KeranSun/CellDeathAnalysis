# =============================================================================
# CellDeathAnalysis - Publication-Ready Figures
# =============================================================================
#
# This script generates all figures for the manuscript using example data
# Can be adapted for TCGA or other real datasets
#
# Author: Keran Sun
# Email: s1214844197@163.com
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

# Install required packages
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "patchwork", "ggrepel",
  "RColorBrewer", "viridis", "scales", "grid", "gridExtra",
  "pheatmap", "corrplot", "ggpubr", "survival", "survminer",
  "UpSetR", "circlize"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(scales)
library(pheatmap)
library(corrplot)
library(ggpubr)
library(survival)
library(survminer)
library(CellDeathAnalysis)

# Set theme
theme_set(theme_bw(base_size = 12))

# Create output directory
dir.create("figures", showWarnings = FALSE)

# Color palettes
death_colors <- c(
  "ferroptosis" = "#E41A1C",
  "cuproptosis" = "#FF7F00",
  "disulfidptosis" = "#FFFF33",
  "pyroptosis" = "#4DAF4A",
  "necroptosis" = "#377EB8",
  "apoptosis" = "#984EA3",
  "autophagy" = "#A65628",
  "panoptosis" = "#F781BF",
  "netosis" = "#999999",
  "parthanatos" = "#66C2A5",
  "entosis" = "#FC8D62",
  "oxeiptosis" = "#8DA0CB",
  "alkaliptosis" = "#E78AC3",
  "ldcd" = "#A6D854"
)

group_colors <- c("Normal" = "#4DAF4A", "Tumor" = "#E41A1C")

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------

cat("Loading example data...\n")
data(example_expr)
data(example_clinical)

# Calculate scores
scores <- calculate_death_score(example_expr, method = "zscore", verbose = FALSE)
pathway_cols <- colnames(scores)

cat("Data loaded:", nrow(example_expr), "genes,", ncol(example_expr), "samples\n")

# =============================================================================
# FIGURE 1: Package Overview (Workflow Diagram)
# =============================================================================

cat("\n=== Creating Figure 1: Package Overview ===\n")

# Create workflow diagram using ggplot
fig1_data <- data.frame(
  x = c(1, 2, 3, 4, 5, 2, 3, 4),
  y = c(3, 3, 3, 3, 3, 1.5, 1.5, 1.5),
  label = c(
    "Input Data\n(Expression Matrix)",
    "Gene Sets\n(14 Death Types)",
    "Scoring\n(6 Methods)",
    "Analysis\n(Survival/Enrichment)",
    "Output\n(Results/Figures)"
  ,
    "Single-Cell\nModule",
    "Machine Learning\nModule",
    "Shiny App\nModule"
  ),
  type = c("input", "process", "process", "process", "output",
           "module", "module", "module")
)

fig1_arrows <- data.frame(
  x = c(1.3, 2.3, 3.3, 4.3),
  xend = c(1.7, 2.7, 3.7, 4.7),
  y = c(3, 3, 3, 3),
  yend = c(3, 3, 3, 3)
)

fig1 <- ggplot() +
  # Main flow boxes
  geom_tile(data = fig1_data[1:5,], 
            aes(x = x, y = y, fill = type),
            width = 0.8, height = 0.8, color = "black", size = 1) +
  # Module boxes
  geom_tile(data = fig1_data[6:8,],
            aes(x = x, y = y, fill = type),
            width = 0.8, height = 0.6, color = "black", size = 0.5) +
  # Labels
  geom_text(data = fig1_data, aes(x = x, y = y, label = label),
            size = 3.5, fontface = "bold") +
  # Arrows
 geom_segment(data = fig1_arrows,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               size = 1, color = "gray30") +
  # Vertical connections
  geom_segment(aes(x = 3, xend = 3, y = 2.6, yend = 1.8),
               linetype = "dashed", color = "gray50") +
  geom_segment(aes(x = 2, xend = 2, y = 2.6, yend = 1.8),
               linetype = "dashed", color = "gray50") +
  geom_segment(aes(x = 4, xend = 4, y = 2.6, yend = 1.8),
               linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c(
    "input" = "#E8F4F8",
    "process" = "#FFF2CC",
    "output" = "#D5E8D4",
    "module" = "#F8CECC"
  )) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "CellDeathAnalysis Package Workflow") +
  coord_fixed(ratio = 0.8) +
  xlim(0.3, 5.7) + ylim(1, 3.6)

ggsave("figures/Figure1_workflow.pdf", fig1, width = 12, height = 6)
ggsave("figures/Figure1_workflow.png", fig1, width = 12, height = 6, dpi = 300)
cat("Figure 1 saved.\n")

# =============================================================================
# FIGURE 2: Gene Set Statistics
# =============================================================================

cat("\n=== Creating Figure 2: Gene Set Statistics ===\n")

# 2A: Bar plot of gene counts
pathway_info <- list_death_pathways(detailed = TRUE)
pathway_info$pathway <- factor(pathway_info$pathway, 
                                levels = pathway_info$pathway[order(pathway_info$total_genes)])

fig2a <- ggplot(pathway_info, aes(x = pathway, y = total_genes, fill = pathway)) +
  geom_col(color = "black", size = 0.3) +
  geom_text(aes(label = total_genes), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = death_colors) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 11),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Number of Genes",
    title = "A. Gene Count per Pathway"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# 2B: Gene overlap heatmap
genesets <- get_death_geneset("all", type = "simple")
pathways <- names(genesets)

# Calculate overlap matrix
overlap_matrix <- matrix(0, length(pathways), length(pathways))
rownames(overlap_matrix) <- colnames(overlap_matrix) <- pathways

for (i in seq_along(pathways)) {
  for (j in seq_along(pathways)) {
    overlap_matrix[i, j] <- length(intersect(genesets[[i]], genesets[[j]]))
  }
}

# Convert to percentage of smaller set
overlap_pct <- overlap_matrix
for (i in seq_along(pathways)) {
  for (j in seq_along(pathways)) {
    min_size <- min(length(genesets[[i]]), length(genesets[[j]]))
    overlap_pct[i, j] <- round(overlap_matrix[i, j] / min_size * 100, 1)
  }
}

# Create heatmap data
overlap_long <- as.data.frame(as.table(overlap_matrix))
names(overlap_long) <- c("Pathway1", "Pathway2", "Overlap")

fig2b <- ggplot(overlap_long, aes(x = Pathway1, y = Pathway2, fill = Overlap)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(Overlap > 0, Overlap, "")), size = 2.5) +
  scale_fill_gradient2(low = "white", mid = "#FFEDA0", high = "#E41A1C",
                       midpoint = max(overlap_long$Overlap)/2,
                       name = "Shared\nGenes") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "B. Gene Overlap Between Pathways"
  ) +
  coord_fixed()

# Combine
fig2 <- fig2a + fig2b + plot_layout(widths = c(1, 1.2))

ggsave("figures/Figure2_geneset_stats.pdf", fig2, width = 14, height = 6)
ggsave("figures/Figure2_geneset_stats.png", fig2, width = 14, height = 6, dpi = 300)
cat("Figure 2 saved.\n")

# =============================================================================
# FIGURE 3: Expression Analysis
# =============================================================================

cat("\n=== Creating Figure 3: Expression Analysis ===\n")

group <- example_clinical$group

# 3A: Heatmap of pathway scores
score_matrix <- as.matrix(scores)
annotation_col <- data.frame(
  Group = group,
  row.names = rownames(scores)
)

# Order by group
order_idx <- order(group)
score_matrix_ordered <- score_matrix[order_idx, ]
annotation_col_ordered <- data.frame(
  Group = group[order_idx],
  row.names = rownames(score_matrix_ordered)
)

# Scale for visualization
score_scaled <- scale(score_matrix_ordered)

pdf("figures/Figure3A_heatmap.pdf", width = 10, height = 8)
pheatmap(
  t(score_scaled),
  annotation_col = annotation_col_ordered,
  annotation_colors = list(Group = group_colors),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 10,
  main = "A. Cell Death Pathway Scores Across Samples",
  breaks = seq(-2, 2, length.out = 101)
)
dev.off()

# 3B: Boxplot comparison
score_long <- scores %>%
  mutate(sample_id = rownames(scores), group = group) %>%
  pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "score")

# Calculate p-values
p_values <- score_long %>%
  group_by(pathway) %>%
  summarise(
    p = wilcox.test(score[group == "Tumor"], score[group == "Normal"])$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))),
    y_pos = max(score_long$score) * 1.1
  )

fig3b <- ggplot(score_long, aes(x = pathway, y = score, fill = group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.7) +
  geom_text(data = p_values, aes(x = pathway, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 4) +
  scale_fill_manual(values = group_colors, name = "Group") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top"
  ) +
  labs(
    x = NULL,
    y = "Pathway Score (Z-score)",
    title = "B. Pathway Score Comparison: Tumor vs Normal"
  )

ggsave("figures/Figure3B_boxplot.pdf", fig3b, width = 12, height = 6)
ggsave("figures/Figure3B_boxplot.png", fig3b, width = 12, height = 6, dpi = 300)

# 3C: Radar chart
score_summary <- score_long %>%
  group_by(pathway, group) %>%
  summarise(mean_score = mean(score), .groups = "drop")

# Normalize for radar
score_summary <- score_summary %>%
  group_by(pathway) %>%
  mutate(norm_score = (mean_score - min(mean_score)) / (max(mean_score) - min(mean_score) + 0.01))

# Create radar plot data
radar_data <- score_summary %>%
  select(pathway, group, norm_score) %>%
  pivot_wider(names_from = group, values_from = norm_score)

# For ggplot radar
n_pathway <- length(unique(score_summary$pathway))
angles <- seq(0, 2*pi, length.out = n_pathway + 1)[1:n_pathway]

score_summary$angle <- rep(angles, 2)
score_summary$x <- score_summary$norm_score * cos(score_summary$angle)
score_summary$y <- score_summary$norm_score * sin(score_summary$angle)

# Labels
label_data <- data.frame(
  pathway = unique(score_summary$pathway),
  angle = angles,
  x = 1.15 * cos(angles),
  y = 1.15 * sin(angles)
)

fig3c <- ggplot() +
  # Grid circles
  geom_path(data = data.frame(
    x = c(cos(seq(0, 2*pi, length.out = 100))),
    y = c(sin(seq(0, 2*pi, length.out = 100)))
  ), aes(x = x * 0.5, y = y * 0.5), color = "gray80", linetype = "dashed") +
  geom_path(data = data.frame(
    x = c(cos(seq(0, 2*pi, length.out = 100))),
    y = c(sin(seq(0, 2*pi, length.out = 100)))
  ), aes(x = x, y = y), color = "gray80", linetype = "dashed") +
  # Grid lines
  geom_segment(data = label_data, aes(x = 0, y = 0, xend = x/1.15, yend = y/1.15),
               color = "gray80") +
  # Polygon
  geom_polygon(data = score_summary, aes(x = x, y = y, fill = group, group = group),
               alpha = 0.3) +
  geom_path(data = score_summary, aes(x = x, y = y, color = group, group = group),
            size = 1) +
  geom_point(data = score_summary, aes(x = x, y = y, color = group), size = 2) +
  # Labels
  geom_text(data = label_data, aes(x = x, y = y, label = pathway), size = 3) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "C. Radar Plot of Mean Pathway Scores", fill = "Group", color = "Group")

ggsave("figures/Figure3C_radar.pdf", fig3c, width = 8, height = 8)
ggsave("figures/Figure3C_radar.png", fig3c, width = 8, height = 8, dpi = 300)

# Combine Figure 3
png("figures/Figure3_combined.png", width = 16, height = 12, units = "in", res = 300)
grid.arrange(
  grid::rasterGrob(png::readPNG("figures/Figure3A_heatmap.pdf")),
  fig3b, fig3c,
  layout_matrix = rbind(c(1, 1), c(2, 3)),
  heights = c(1, 1)
)
dev.off()

cat("Figure 3 saved.\n")

# =============================================================================
# FIGURE 4: Survival Analysis
# =============================================================================

cat("\n=== Creating Figure 4: Survival Analysis ===\n")

# Use tumor samples only
tumor_idx <- example_clinical$group == "Tumor"
tumor_scores <- scores[tumor_idx, ]
tumor_clinical <- example_clinical[tumor_idx, ]

# 4A: Forest plot for all pathways
surv_results <- data.frame()

for (pathway in pathway_cols) {
  tryCatch({
    score_vec <- tumor_scores[[pathway]]
    
    cox_data <- data.frame(
      time = tumor_clinical$OS_time,
      status = tumor_clinical$OS_status,
      score = scale(score_vec)[,1]
    )
    cox_data <- cox_data[complete.cases(cox_data), ]
    
    fit <- coxph(Surv(time, status) ~ score, data = cox_data)
    s <- summary(fit)
    
    surv_results <- rbind(surv_results, data.frame(
      pathway = pathway,
      hr = s$conf.int[1, 1],
      hr_lower = s$conf.int[1, 3],
      hr_upper = s$conf.int[1, 4],
      p = s$coefficients[1, 5]
    ))
  }, error = function(e) {})
}

surv_results$p_adj <- p.adjust(surv_results$p, method = "BH")
surv_results$sig <- ifelse(surv_results$p_adj < 0.05, "Significant", "Not Significant")
surv_results <- surv_results[order(surv_results$hr), ]
surv_results$pathway <- factor(surv_results$pathway, levels = surv_results$pathway)

fig4a <- ggplot(surv_results, aes(x = hr, y = pathway)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbarh(aes(xmin = hr_lower, xmax = hr_upper, color = sig),
                 height = 0.25, size = 0.8) +
  geom_point(aes(color = sig, size = -log10(p_adj)), shape = 18) +
  scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "gray50"),
                     name = "Significance") +
  scale_size_continuous(range = c(3, 8), name = "-log10(p.adj)") +
  scale_x_log10(breaks = c(0.5, 1, 2, 4)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Hazard Ratio (95% CI, log scale)",
    y = NULL,
    title = "A. Forest Plot: Cell Death Pathways and Overall Survival"
  )

ggsave("figures/Figure4A_forest.pdf", fig4a, width = 10, height = 8)
ggsave("figures/Figure4A_forest.png", fig4a, width = 10, height = 8, dpi = 300)

# 4B: Kaplan-Meier curves for top pathways
top_pathways <- c("ferroptosis", "pyroptosis", "cuproptosis", "apoptosis")
km_plots <- list()

for (pathway in top_pathways) {
  score_vec <- tumor_scores[[pathway]]
  
  km_data <- data.frame(
    time = tumor_clinical$OS_time,
    status = tumor_clinical$OS_status,
    group = ifelse(score_vec > median(score_vec), "High", "Low")
  )
  km_data <- km_data[complete.cases(km_data), ]
  km_data$group <- factor(km_data$group, levels = c("Low", "High"))
  
  fit <- survfit(Surv(time, status) ~ group, data = km_data)
  
  km_plots[[pathway]] <- ggsurvplot(
    fit, data = km_data,
    pval = TRUE, pval.size = 4,
    risk.table = FALSE,
    palette = c("#2E9FDF", "#E7B800"),
    legend.title = paste0(tools::toTitleCase(pathway), " Score"),
    legend.labs = c("Low", "High"),
    xlab = "Time (months)",
    ylab = "Survival Probability",
    title = tools::toTitleCase(pathway),
    ggtheme = theme_bw(base_size = 11),
    font.title = c(12, "bold"),
    conf.int = TRUE,
    conf.int.alpha = 0.1
  )
}

# Combine KM plots
pdf("figures/Figure4B_KM_curves.pdf", width = 12, height = 10)
arrange_ggsurvplots(km_plots, ncol = 2, nrow = 2, 
                    title = "B. Kaplan-Meier Survival Curves")
dev.off()

# Also save as PNG
png("figures/Figure4B_KM_curves.png", width = 12, height = 10, units = "in", res = 300)
arrange_ggsurvplots(km_plots, ncol = 2, nrow = 2,
                    title = "B. Kaplan-Meier Survival Curves")
dev.off()

cat("Figure 4 saved.\n")

# =============================================================================
# FIGURE 5: Pathway Correlation
# =============================================================================

cat("\n=== Creating Figure 5: Pathway Correlation ===\n")

cor_matrix <- cor(scores, method = "spearman")

# Create correlation plot
pdf("figures/Figure5_correlation.pdf", width = 10, height = 10)
corrplot(
  cor_matrix,
  method = "color",
  type = "lower",
  order = "hclust",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.9,
  addCoef.col = "black",
  number.cex = 0.7,
  col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  title = "Pathway Correlation (Spearman)",
  mar = c(0, 0, 2, 0)
)
dev.off()

# Also as ggplot version
cor_long <- as.data.frame(as.table(cor_matrix))
names(cor_long) <- c("Pathway1", "Pathway2", "Correlation")

fig5 <- ggplot(cor_long, aes(x = Pathway1, y = Pathway2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman\nCorrelation") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL, title = "Cell Death Pathway Correlation") +
  coord_fixed()

ggsave("figures/Figure5_correlation_ggplot.pdf", fig5, width = 10, height = 9)
ggsave("figures/Figure5_correlation_ggplot.png", fig5, width = 10, height = 9, dpi = 300)

cat("Figure 5 saved.\n")

# =============================================================================
# FIGURE 6: Score Distribution and Classification
# =============================================================================

cat("\n=== Creating Figure 6: Score Distribution ===\n")

# 6A: Density plots
fig6a <- ggplot(score_long %>% filter(pathway %in% top_pathways),
                aes(x = score, fill = group, color = group)) +
  geom_density(alpha = 0.4, size = 0.8) +
  facet_wrap(~pathway, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Pathway Score",
    y = "Density",
    title = "A. Score Distribution by Group",
    fill = "Group", color = "Group"
  )

# 6B: Violin plots
fig6b <- ggplot(score_long %>% filter(pathway %in% top_pathways),
                aes(x = pathway, y = score, fill = group)) +
  geom_violin(position = position_dodge(0.8), alpha = 0.7) +
  geom_boxplot(position = position_dodge(0.8), width = 0.15, 
               outlier.size = 0.5, alpha = 0.9) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(aes(group = group), method = "wilcox.test",
                     label = "p.signif", size = 4) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 11)
  ) +
  labs(
    x = NULL,
    y = "Pathway Score",
    title = "B. Violin Plot with Statistical Comparison",
    fill = "Group"
  )

fig6 <- fig6a / fig6b
ggsave("figures/Figure6_distribution.pdf", fig6, width = 10, height = 12)
ggsave("figures/Figure6_distribution.png", fig6, width = 10, height = 12, dpi = 300)

cat("Figure 6 saved.\n")

# =============================================================================
# FIGURE 7: Enrichment Analysis Demo
# =============================================================================

cat("\n=== Creating Figure 7: Enrichment Analysis ===\n")

# Simulate DEGs (genes with high variance + death pathway genes)
set.seed(123)
all_death_genes <- get_all_death_genes()
simulated_degs <- c(
  sample(rownames(example_expr)[rownames(example_expr) %in% all_death_genes], 30),
  sample(rownames(example_expr), 70)
)
simulated_degs <- unique(simulated_degs)

# ORA analysis
ora_result <- death_enrich_ora(simulated_degs, min_genes = 2)

if (!is.null(ora_result) && nrow(ora_result) > 0) {
  ora_result <- ora_result[order(ora_result$p_value), ]
  ora_result$pathway <- factor(ora_result$pathway, 
                                levels = rev(ora_result$pathway))
  
  # 7A: Bar plot
  fig7a <- ggplot(ora_result, aes(x = -log10(p_adjust), y = pathway)) +
    geom_col(aes(fill = fold_enrichment), color = "black", size = 0.3) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_fill_gradient(low = "#FEE0D2", high = "#DE2D26", name = "Fold\nEnrichment") +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 11)) +
    labs(
      x = "-log10(Adjusted p-value)",
      y = NULL,
      title = "A. Over-Representation Analysis"
    )
  
  # 7B: Dot plot
  fig7b <- ggplot(ora_result, aes(x = fold_enrichment, y = pathway)) +
    geom_point(aes(size = n_overlap, color = -log10(p_adjust))) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10\n(p.adj)") +
    scale_size_continuous(range = c(3, 10), name = "Gene\nCount") +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 11)) +
    labs(
      x = "Fold Enrichment",
      y = NULL,
      title = "B. Enrichment Dot Plot"
    )
  
  fig7 <- fig7a + fig7b
  ggsave("figures/Figure7_enrichment.pdf", fig7, width = 14, height = 6)
  ggsave("figures/Figure7_enrichment.png", fig7, width = 14, height = 6, dpi = 300)
}

cat("Figure 7 saved.\n")

# =============================================================================
# SUPPLEMENTARY FIGURES
# =============================================================================

cat("\n=== Creating Supplementary Figures ===\n")

# Figure S1: Scoring method comparison
methods <- c("zscore", "mean", "median")
method_scores <- list()

for (m in methods) {
  method_scores[[m]] <- calculate_death_score(example_expr, method = m, verbose = FALSE)
}

# Compare ferroptosis scores across methods
method_comparison <- data.frame(
  sample = rownames(method_scores[["zscore"]]),
  zscore = method_scores[["zscore"]]$ferroptosis,
  mean = method_scores[["mean"]]$ferroptosis,
  median = method_scores[["median"]]$ferroptosis
)

cor_methods <- cor(method_comparison[, -1])

figS1 <- ggplot(method_comparison, aes(x = zscore, y = mean)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  annotate("text", x = min(method_comparison$zscore), 
           y = max(method_comparison$mean),
           label = paste0("r = ", round(cor(method_comparison$zscore, 
                                            method_comparison$mean), 3)),
           hjust = 0, size = 5) +
  theme_bw() +
  labs(
    x = "Z-score Method",
    y = "Mean Method",
    title = "Figure S1: Scoring Method Comparison (Ferroptosis)"
  )

ggsave("figures/FigureS1_method_comparison.pdf", figS1, width = 8, height = 7)

# =============================================================================
# Summary Table
# =============================================================================

cat("\n=== Creating Summary Tables ===\n")

# Table 1: Gene set summary
table1 <- pathway_info %>%
  select(pathway, total_genes, core_genes, year_discovered) %>%
  arrange(desc(total_genes))

write.csv(table1, "figures/Table1_geneset_summary.csv", row.names = FALSE)

# Table 2: Survival analysis results
write.csv(surv_results, "figures/Table2_survival_results.csv", row.names = FALSE)

# =============================================================================
# Final Summary
# =============================================================================

cat("\n")
cat("=============================================================\n")
cat("              Figure Generation Complete!                    \n")
cat("=============================================================\n")
cat("\n")
cat("Generated files:\n")
cat("  Main Figures:\n")
cat("    - Figure1_workflow.pdf/png\n")
cat("    - Figure2_geneset_stats.pdf/png\n")
cat("    - Figure3A_heatmap.pdf\n")
cat("    - Figure3B_boxplot.pdf/png\n")
cat("    - Figure3C_radar.pdf/png\n")
cat("    - Figure4A_forest.pdf/png\n")
cat("    - Figure4B_KM_curves.pdf/png\n")
cat("    - Figure5_correlation.pdf/png\n")
cat("    - Figure6_distribution.pdf/png\n")
cat("    - Figure7_enrichment.pdf/png\n")
cat("\n")
cat("  Supplementary:\n")
cat("    - FigureS1_method_comparison.pdf\n")
cat("\n")
cat("  Tables:\n")
cat("    - Table1_geneset_summary.csv\n")
cat("    - Table2_survival_results.csv\n")
cat("\n")
cat("All files saved in ./figures/ directory\n")
cat("=============================================================\n")
