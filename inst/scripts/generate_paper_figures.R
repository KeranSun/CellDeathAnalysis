# =============================================================================
# CellDeathAnalysis - Publication-Ready Figures Generation
# =============================================================================
#
# This script generates all figures for the manuscript
# Uses example data + simulated TCGA-like data to ensure all figures work
#
# Author: Keran Sun
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

# Install and load packages
packages <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "scales",
              "gridExtra", "grid", "ggpubr", "patchwork", "reshape2",
              "survival", "survminer", "pheatmap", "corrplot", "ggrepel",
              "UpSetR", "circlize")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Load CellDeathAnalysis
library(CellDeathAnalysis)

# Create output directory
dir.create("figures", showWarnings = FALSE)

# Set theme for all plots
theme_publication <- theme_bw() +
  theme(
    text = element_text(size = 12, family = "sans"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 11, face = "bold")
  )

theme_set(theme_publication)

# Color palettes
cancer_colors <- c(
  "BRCA" = "#E41A1C", "LUAD" = "#377EB8", "LIHC" = "#4DAF4A",
  "STAD" = "#984EA3", "COAD" = "#FF7F00", "KIRC" = "#FFFF33",
  "HNSC" = "#A65628", "THCA" = "#F781BF", "PRAD" = "#999999",
  "BLCA" = "#66C2A5", "LUSC" = "#FC8D62", "KIRP" = "#8DA0CB"
)

pathway_colors <- c(
  "ferroptosis" = "#E41A1C", "cuproptosis" = "#FF7F00", 
  "disulfidptosis" = "#FFFF33", "pyroptosis" = "#984EA3",
  "necroptosis" = "#4DAF4A", "apoptosis" = "#377EB8",
  "autophagy" = "#A65628", "panoptosis" = "#F781BF",
  "netosis" = "#999999", "parthanatos" = "#66C2A5",
  "entosis" = "#FC8D62", "oxeiptosis" = "#8DA0CB",
  "alkaliptosis" = "#E5C494", "ldcd" = "#B3B3B3"
)

cat("Setup complete!\n\n")

# -----------------------------------------------------------------------------
# 1. Prepare Data
# -----------------------------------------------------------------------------

cat("=== Preparing Data ===\n")

# Load example data
data(example_expr)
data(example_clinical)

# Get pathway information
pathway_info <- list_death_pathways(detailed = TRUE)
pathway_cols <- pathway_info$pathway

# Calculate scores
scores <- calculate_death_score(example_expr, method = "zscore", verbose = FALSE)

cat("Example data loaded\n")
cat("  Samples:", nrow(scores), "\n")
cat("  Pathways:", ncol(scores), "\n\n")

# -----------------------------------------------------------------------------
# Simulate TCGA-like pan-cancer data for demonstration
# -----------------------------------------------------------------------------

cat("=== Simulating Pan-Cancer Data ===\n")

set.seed(42)

# Cancer types
cancers <- c("BRCA", "LUAD", "LIHC", "STAD", "COAD", "KIRC", 
             "HNSC", "THCA", "PRAD", "BLCA", "LUSC", "KIRP")

# Simulate scores for each cancer
pancancer_scores <- data.frame()
pancancer_clinical <- data.frame()

for (cancer in cancers) {
  n_tumor <- sample(150:300, 1)
  n_normal <- sample(20:50, 1)
  
  # Simulate pathway scores with cancer-specific patterns
  for (type in c("Tumor", "Normal")) {
    n <- if(type == "Tumor") n_tumor else n_normal
    
    temp_scores <- data.frame(
      sample_id = paste0(cancer, "_", type, "_", 1:n),
      cancer = cancer,
      sample_type = type
    )
    
    for (pathway in pathway_cols) {
      # Different patterns for different pathways
      base_mean <- rnorm(1, 0, 0.3)
      tumor_effect <- if(type == "Tumor") {
        switch(pathway,
               "ferroptosis" = rnorm(1, 0.5, 0.2),
               "pyroptosis" = rnorm(1, 0.4, 0.2),
               "cuproptosis" = rnorm(1, 0.3, 0.2),
               "apoptosis" = rnorm(1, -0.2, 0.2),
               rnorm(1, 0.1, 0.2))
      } else 0
      
      temp_scores[[pathway]] <- rnorm(n, base_mean + tumor_effect, 0.8)
    }
    
    pancancer_scores <- rbind(pancancer_scores, temp_scores)
  }
  
  # Simulate survival data for tumor samples
  n_tumor_surv <- n_tumor
  temp_clinical <- data.frame(
    sample_id = paste0(cancer, "_Tumor_", 1:n_tumor_surv),
    cancer = cancer,
    OS_time = rexp(n_tumor_surv, rate = 0.02) * 30,  # months
    OS_status = rbinom(n_tumor_surv, 1, 0.4),
    age = rnorm(n_tumor_surv, 60, 12),
    stage = sample(c("I", "II", "III", "IV"), n_tumor_surv, 
                   replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))
  )
  
  pancancer_clinical <- rbind(pancancer_clinical, temp_clinical)
}

# Merge scores with clinical for tumor samples
tumor_scores <- pancancer_scores %>% filter(sample_type == "Tumor")
tumor_data <- merge(tumor_scores, pancancer_clinical, by = c("sample_id", "cancer"))

cat("Simulated pan-cancer data:\n")
cat("  Total samples:", nrow(pancancer_scores), "\n")
cat("  Tumor samples:", sum(pancancer_scores$sample_type == "Tumor"), "\n")
cat("  Normal samples:", sum(pancancer_scores$sample_type == "Normal"), "\n")
cat("  Cancer types:", length(cancers), "\n\n")

# =============================================================================
# FIGURE 2: Gene Set Statistics
# =============================================================================

cat("=== Figure 2: Gene Set Statistics ===\n")

# Figure 2A: Bar plot of gene counts
gene_counts <- data.frame(
  pathway = pathway_cols,
  n_genes = sapply(pathway_cols, function(p) length(get_death_geneset(p, type = "all")))
)
gene_counts$pathway <- factor(gene_counts$pathway, 
                               levels = gene_counts$pathway[order(gene_counts$n_genes, decreasing = TRUE)])

fig2a <- ggplot(gene_counts, aes(x = pathway, y = n_genes, fill = pathway)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = n_genes), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = pathway_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = NULL,
    y = "Number of Genes",
    title = "A. Gene Set Sizes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  )

# Figure 2B: UpSet plot for pathway overlap
# Get gene lists
gene_lists <- lapply(pathway_cols[1:8], function(p) {
  get_death_geneset(p, type = "all")
})
names(gene_lists) <- pathway_cols[1:8]

# Create binary matrix for UpSet
all_genes <- unique(unlist(gene_lists))
upset_matrix <- data.frame(
  gene = all_genes,
  stringsAsFactors = FALSE
)
for (pathway in names(gene_lists)) {
  upset_matrix[[pathway]] <- as.integer(upset_matrix$gene %in% gene_lists[[pathway]])
}

pdf("figures/Figure2B_upset.pdf", width = 10, height = 6)
upset(upset_matrix[, -1], 
      sets = names(gene_lists),
      order.by = "freq",
      sets.bar.color = unname(pathway_colors[names(gene_lists)]),
      main.bar.color = "steelblue",
      matrix.color = "steelblue",
      text.scale = 1.3)
dev.off()

# Figure 2C: Pathway information table visualization
fig2c_data <- pathway_info %>%
  select(pathway, year_discovered, total_genes, core_genes) %>%
  mutate(pathway = factor(pathway, levels = rev(pathway_cols)))

fig2c <- ggplot(fig2c_data, aes(x = total_genes, y = pathway, fill = pathway)) +
  geom_col(show.legend = FALSE) +
  geom_point(aes(x = core_genes), color = "red", size = 3) +
  scale_fill_manual(values = pathway_colors) +
  labs(
    x = "Number of Genes",
    y = NULL,
    title = "C. Total Genes (bars) vs Core Genes (red dots)"
  )

# Combine 2A and 2C
fig2_combined <- fig2a + fig2c + plot_layout(ncol = 2, widths = c(1.2, 1))

ggsave("figures/Figure2AC_genesets.pdf", fig2_combined, width = 14, height = 6)

cat("Figure 2 saved!\n\n")

# =============================================================================
# FIGURE 3: Pan-Cancer Analysis
# =============================================================================

cat("=== Figure 3: Pan-Cancer Analysis ===\n")

# Figure 3A: Pan-cancer heatmap
mean_scores_tumor <- pancancer_scores %>%
  filter(sample_type == "Tumor") %>%
  group_by(cancer) %>%
  summarise(across(all_of(pathway_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("cancer") %>%
  as.matrix()

# Scale by column (pathway)
mean_scores_scaled <- scale(mean_scores_tumor)

# Custom color palette
heatmap_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

pdf("figures/Figure3A_pancancer_heatmap.pdf", width = 12, height = 8)
pheatmap(
  mean_scores_scaled,
  color = heatmap_colors,
  clustering_method = "ward.D2",
  border_color = "white",
  cellwidth = 30,
  cellheight = 25,
  fontsize = 11,
  fontsize_row = 12,
  fontsize_col = 11,
  angle_col = 45,
  main = "Cell Death Pathway Scores Across Cancer Types"
)
dev.off()

# Figure 3B: Tumor vs Normal fold change heatmap
fc_results <- pancancer_scores %>%
  group_by(cancer, sample_type) %>%
  summarise(across(all_of(pathway_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "score") %>%
  pivot_wider(names_from = sample_type, values_from = score) %>%
  mutate(log2FC = log2((Tumor + 0.01) / (Normal + 0.01)))

fc_matrix <- fc_results %>%
  select(cancer, pathway, log2FC) %>%
  pivot_wider(names_from = pathway, values_from = log2FC) %>%
  column_to_rownames("cancer") %>%
  as.matrix()

fc_colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

pdf("figures/Figure3B_tumor_normal_FC.pdf", width = 12, height = 8)
pheatmap(
  fc_matrix,
  color = fc_colors,
  clustering_method = "ward.D2",
  border_color = "white",
  cellwidth = 30,
  cellheight = 25,
  fontsize = 11,
  breaks = seq(-1.5, 1.5, length.out = 101),
  main = "Tumor vs Normal (log2 Fold Change)"
)
dev.off()

# Figure 3C: Boxplots for key pathways
key_pathways <- c("ferroptosis", "cuproptosis", "pyroptosis", "apoptosis")

boxplot_data <- pancancer_scores %>%
  select(cancer, sample_type, all_of(key_pathways)) %>%
  pivot_longer(cols = all_of(key_pathways), names_to = "pathway", values_to = "score") %>%
  mutate(pathway = factor(pathway, levels = key_pathways))

fig3c <- ggplot(boxplot_data, aes(x = cancer, y = score, fill = sample_type)) +
  geom_boxplot(outlier.size = 0.5, width = 0.7) +
  facet_wrap(~pathway, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Normal" = "#4DAF4A", "Tumor" = "#E41A1C"),
                    name = "Sample Type") +
  labs(x = NULL, y = "Death Score (Z-score)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "top"
  )

ggsave("figures/Figure3C_boxplots.pdf", fig3c, width = 14, height = 10)

# Figure 3D: Summary bar chart showing significant changes
sig_summary <- fc_results %>%
  group_by(pathway) %>%
  summarise(
    n_up = sum(log2FC > 0.3, na.rm = TRUE),
    n_down = sum(log2FC < -0.3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(n_up, n_down), names_to = "direction", values_to = "count") %>%
  mutate(
    count = ifelse(direction == "n_down", -count, count),
    direction = factor(direction, levels = c("n_up", "n_down"),
                       labels = c("Upregulated in Tumor", "Downregulated in Tumor"))
  )

fig3d <- ggplot(sig_summary, aes(x = reorder(pathway, -count), y = count, fill = direction)) +
  geom_col(position = "identity") +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = c("Upregulated in Tumor" = "#E41A1C", 
                                "Downregulated in Tumor" = "#377EB8")) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of Cancer Types",
    fill = NULL,
    title = "D. Pathway Dysregulation Summary"
  ) +
  theme(legend.position = "bottom")

ggsave("figures/Figure3D_summary.pdf", fig3d, width = 8, height = 6)

cat("Figure 3 saved!\n\n")

# =============================================================================
# FIGURE 4: Survival Analysis
# =============================================================================

cat("=== Figure 4: Survival Analysis ===\n")

# Perform survival analysis for all cancer-pathway combinations
survival_results <- data.frame()

for (cancer in cancers) {
  cancer_data <- tumor_data %>% filter(cancer == !!cancer)
  
  if (nrow(cancer_data) < 50) next
  
  for (pathway in pathway_cols) {
    tryCatch({
      # Fit Cox model
      cancer_data$pathway_score <- cancer_data[[pathway]]
      cox_fit <- coxph(Surv(OS_time, OS_status) ~ pathway_score, data = cancer_data)
      cox_sum <- summary(cox_fit)
      
      survival_results <- rbind(survival_results, data.frame(
        cancer = cancer,
        pathway = pathway,
        n = nrow(cancer_data),
        hr = cox_sum$conf.int[1, 1],
        hr_lower = cox_sum$conf.int[1, 3],
        hr_upper = cox_sum$conf.int[1, 4],
        p_value = cox_sum$coefficients[1, 5],
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {})
  }
}

survival_results$p_adjust <- p.adjust(survival_results$p_value, method = "BH")
survival_results$significance <- ifelse(survival_results$p_adjust < 0.05, "Significant", "NS")

# Figure 4A: Forest plot for ferroptosis
ferro_surv <- survival_results %>%
  filter(pathway == "ferroptosis") %>%
  arrange(hr)

ferro_surv$cancer <- factor(ferro_surv$cancer, levels = ferro_surv$cancer)

fig4a <- ggplot(ferro_surv, aes(x = hr, y = cancer)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = hr_lower, xmax = hr_upper), height = 0.25, linewidth = 0.7) +
  geom_point(aes(color = significance, size = n)) +
  scale_color_manual(values = c("Significant" = "#E41A1C", "NS" = "grey40")) +
  scale_size_continuous(range = c(2, 5), name = "Sample Size") +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2)) +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = NULL,
    color = "Significance",
    title = "A. Ferroptosis Score and Overall Survival"
  ) +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey90")
  )

ggsave("figures/Figure4A_forest_ferroptosis.pdf", fig4a, width = 10, height = 8)

# Figure 4B: Forest plot for all pathways (selected cancers)
selected_cancers <- c("LIHC", "KIRC", "LUAD", "BRCA")
surv_selected <- survival_results %>%
  filter(cancer %in% selected_cancers)

fig4b <- ggplot(surv_selected, aes(x = hr, y = pathway, color = cancer)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = hr_lower, xmax = hr_upper), 
                 height = 0.3, position = position_dodge(0.6)) +
  geom_point(position = position_dodge(0.6), size = 2.5) +
  scale_color_manual(values = cancer_colors[selected_cancers]) +
  scale_x_log10() +
  facet_wrap(~cancer, ncol = 2) +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = NULL,
    title = "B. All Pathways Across Selected Cancers"
  ) +
  theme(legend.position = "none")

ggsave("figures/Figure4B_forest_multipathway.pdf", fig4b, width = 12, height = 10)

# Figure 4C: Kaplan-Meier curves
km_plots <- list()

for (cancer in c("LIHC", "KIRC", "LUAD", "STAD")) {
  cancer_data <- tumor_data %>% filter(cancer == !!cancer)
  cancer_data$ferro_group <- ifelse(
    cancer_data$ferroptosis > median(cancer_data$ferroptosis),
    "High", "Low"
  )
  cancer_data$ferro_group <- factor(cancer_data$ferro_group, levels = c("Low", "High"))
  
  fit <- survfit(Surv(OS_time, OS_status) ~ ferro_group, data = cancer_data)
  
  km_plots[[cancer]] <- ggsurvplot(
    fit,
    data = cancer_data,
    pval = TRUE,
    pval.size = 4,
    risk.table = TRUE,
    risk.table.height = 0.3,
    palette = c("#377EB8", "#E41A1C"),
    title = cancer,
    xlab = "Time (months)",
    ylab = "Overall Survival",
    legend.title = "Ferroptosis",
    legend.labs = c("Low", "High"),
    font.main = c(14, "bold"),
    font.x = c(12, "bold"),
    font.y = c(12, "bold"),
    font.tickslab = 10,
    ggtheme = theme_publication
  )
}

pdf("figures/Figure4C_KM_curves.pdf", width = 12, height = 12)
arrange_ggsurvplots(km_plots, ncol = 2, nrow = 2)
dev.off()

# Figure 4D: Heatmap of survival significance
surv_heatmap <- survival_results %>%
  mutate(neg_log_p = -log10(p_adjust + 0.001),
         neg_log_p = ifelse(hr > 1, neg_log_p, -neg_log_p)) %>%
  select(cancer, pathway, neg_log_p) %>%
  pivot_wider(names_from = pathway, values_from = neg_log_p) %>%
  column_to_rownames("cancer") %>%
  as.matrix()

surv_colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

pdf("figures/Figure4D_survival_heatmap.pdf", width = 14, height = 8)
pheatmap(
  surv_heatmap,
  color = surv_colors,
  breaks = seq(-3, 3, length.out = 101),
  clustering_method = "ward.D2",
  border_color = "white",
  cellwidth = 28,
  cellheight = 22,
  fontsize = 10,
  main = "Survival Association (-log10 p × direction)"
)
dev.off()

cat("Figure 4 saved!\n\n")

# =============================================================================
# FIGURE 5: Pathway Correlation
# =============================================================================

cat("=== Figure 5: Pathway Correlation ===\n")

# Calculate correlation matrix
score_matrix <- pancancer_scores %>%
  filter(sample_type == "Tumor") %>%
  select(all_of(pathway_cols)) %>%
  as.matrix()

cor_matrix <- cor(score_matrix, use = "pairwise.complete.obs", method = "spearman")

# Figure 5A: Correlation heatmap with values
pdf("figures/Figure5A_correlation_heatmap.pdf", width = 12, height = 10)
corrplot(
  cor_matrix,
  method = "color",
  type = "lower",
  order = "hclust",
  hclust.method = "ward.D2",
  addCoef.col = "black",
  number.cex = 0.7,
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.9,
  col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
  title = "Cell Death Pathway Correlation (Spearman)",
  mar = c(0, 0, 2, 0)
)
dev.off()

# Figure 5B: Correlation network (simplified)
cor_long <- as.data.frame(as.table(cor_matrix))
names(cor_long) <- c("pathway1", "pathway2", "correlation")
cor_long <- cor_long %>%
  filter(as.character(pathway1) < as.character(pathway2),
         abs(correlation) > 0.3)

fig5b <- ggplot(cor_long, aes(x = pathway1, y = pathway2)) +
  geom_tile(aes(fill = correlation), color = "white") +
  geom_text(aes(label = round(correlation, 2)), size = 3) +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(
    x = NULL, y = NULL,
    fill = "Spearman\nCorrelation",
    title = "B. Strong Pathway Correlations (|r| > 0.3)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("figures/Figure5B_correlation_filtered.pdf", fig5b, width = 10, height = 8)

# Figure 5C: Scatter plots for top correlations
top_cors <- cor_long %>%
  arrange(desc(abs(correlation))) %>%
  head(4)

scatter_plots <- list()
for (i in 1:nrow(top_cors)) {
  p1 <- as.character(top_cors$pathway1[i])
  p2 <- as.character(top_cors$pathway2[i])
  r <- round(top_cors$correlation[i], 2)
  
  plot_data <- pancancer_scores %>%
    filter(sample_type == "Tumor") %>%
    select(all_of(c(p1, p2)))
  
  scatter_plots[[i]] <- ggplot(plot_data, aes_string(x = p1, y = p2)) +
    geom_point(alpha = 0.3, size = 1, color = "steelblue") +
    geom_smooth(method = "lm", color = "#E41A1C", se = TRUE) +
    annotate("text", x = Inf, y = Inf, label = paste0("r = ", r),
             hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold") +
    labs(
      x = tools::toTitleCase(p1),
      y = tools::toTitleCase(p2)
    )
}

fig5c <- wrap_plots(scatter_plots, ncol = 2) +
  plot_annotation(title = "C. Top Correlated Pathway Pairs",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave("figures/Figure5C_scatter_plots.pdf", fig5c, width = 10, height = 10)

cat("Figure 5 saved!\n\n")

# =============================================================================
# FIGURE 6: Example Analysis with Package Data
# =============================================================================

cat("=== Figure 6: Package Demo Analysis ===\n")

# Use example data from package
group <- example_clinical$group

# Figure 6A: Radar chart
fig6a <- plot_death_radar(scores, group = group) +
  labs(title = "A. Death Pathway Profile: Tumor vs Normal") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggsave("figures/Figure6A_radar.pdf", fig6a, width = 8, height = 8)

# Figure 6B: Boxplot comparison
boxplot_pathways <- c("ferroptosis", "cuproptosis", "pyroptosis", "apoptosis",
                      "autophagy", "necroptosis")

fig6b <- plot_death_boxplot(scores, group = group, pathways = boxplot_pathways, ncol = 3) +
  plot_annotation(title = "B. Pathway Score Comparison",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave("figures/Figure6B_boxplot.pdf", fig6b, width = 12, height = 8)

# Figure 6C: Score distribution
score_long <- scores %>%
  mutate(group = group) %>%
  pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "score")

fig6c <- ggplot(score_long %>% filter(pathway %in% c("ferroptosis", "pyroptosis")),
                aes(x = score, fill = group, color = group)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  facet_wrap(~pathway, scales = "free") +
  scale_fill_manual(values = c("Normal" = "#4DAF4A", "Tumor" = "#E41A1C")) +
  scale_color_manual(values = c("Normal" = "#4DAF4A", "Tumor" = "#E41A1C")) +
  labs(
    x = "Death Score",
    y = "Density",
    fill = "Group",
    color = "Group",
    title = "C. Score Distribution"
  ) +
  theme(legend.position = "top")

ggsave("figures/Figure6C_distribution.pdf", fig6c, width = 10, height = 5)

# Figure 6D: Sample classification
ferro_class <- classify_by_score(scores$ferroptosis, method = "median")
class_data <- data.frame(
  sample = rownames(scores),
  ferroptosis = scores$ferroptosis,
  group = group,
  class = ferro_class
)

fig6d <- ggplot(class_data, aes(x = ferroptosis, y = group, color = class)) +
  geom_jitter(height = 0.2, size = 2, alpha = 0.7) +
  geom_vline(xintercept = median(scores$ferroptosis), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("Low" = "#377EB8", "High" = "#E41A1C")) +
  labs(
    x = "Ferroptosis Score",
    y = NULL,
    color = "Classification",
    title = "D. Sample Classification by Ferroptosis Score"
  ) +
  theme(legend.position = "top")

ggsave("figures/Figure6D_classification.pdf", fig6d, width = 8, height = 5)

cat("Figure 6 saved!\n\n")

# =============================================================================
# FIGURE 7: Enrichment Analysis
# =============================================================================

cat("=== Figure 7: Enrichment Analysis ===\n")

# Simulate DEGs enriched in death pathways
set.seed(123)
all_death_genes <- get_all_death_genes()
example_degs <- c(
  sample(get_death_geneset("ferroptosis"), 15),
  sample(get_death_geneset("pyroptosis"), 10),
  sample(get_death_geneset("cuproptosis"), 8),
  sample(get_death_geneset("apoptosis"), 12),
  sample(all_death_genes, 20)
)
example_degs <- unique(example_degs)

# ORA analysis
ora_result <- death_enrich_ora(example_degs, min_genes = 3)

# Figure 7A: Enrichment bar plot
ora_plot_data <- ora_result %>%
  mutate(
    pathway = factor(pathway, levels = pathway[order(p_adjust)]),
    neg_log_p = -log10(p_adjust)
  ) %>%
  head(10)

fig7a <- ggplot(ora_plot_data, aes(x = neg_log_p, y = reorder(pathway, neg_log_p))) +
  geom_col(aes(fill = fold_enrichment), width = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Fold\nEnrichment") +
  labs(
    x = "-log10(Adjusted p-value)",
    y = NULL,
    title = "A. Over-Representation Analysis"
  )

ggsave("figures/Figure7A_ORA.pdf", fig7a, width = 10, height = 6)

# Figure 7B: Enrichment dot plot
fig7b <- ggplot(ora_plot_data, aes(x = fold_enrichment, y = reorder(pathway, fold_enrichment))) +
  geom_segment(aes(xend = 0, yend = pathway), color = "grey70") +
  geom_point(aes(size = n_overlap, color = neg_log_p)) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
  scale_size_continuous(range = c(3, 10), name = "Gene Count") +
  labs(
    x = "Fold Enrichment",
    y = NULL,
    title = "B. Enrichment Overview"
  )

ggsave("figures/Figure7B_dotplot.pdf", fig7b, width = 10, height = 6)

# Figure 7C: Group comparison volcano
compare_result <- death_compare_groups(example_expr, group = group, method = "wilcox")

compare_result <- compare_result %>%
  mutate(
    neg_log_p = -log10(p_adjust),
    sig_label = ifelse(p_adjust < 0.05 & abs(log2FC) > 0.3, pathway, "")
  )

fig7c <- ggplot(compare_result, aes(x = log2FC, y = neg_log_p)) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_point(aes(color = p_adjust < 0.05 & abs(log2FC) > 0.3), size = 4) +
  geom_text_repel(aes(label = sig_label), size = 4, max.overlaps = 20) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "#E41A1C"),
                     guide = "none") +
  labs(
    x = "log2(Fold Change) Tumor vs Normal",
    y = "-log10(Adjusted p-value)",
    title = "C. Pathway Differential Analysis"
  )

ggsave("figures/Figure7C_volcano.pdf", fig7c, width = 8, height = 7)

cat("Figure 7 saved!\n\n")

# =============================================================================
# SUPPLEMENTARY: Gene Table
# =============================================================================

cat("=== Generating Supplementary Tables ===\n")

# Table S1: Complete gene list
gene_table <- data.frame()
for (pathway in pathway_cols) {
  genes <- get_death_geneset(pathway, type = "all")
  temp <- data.frame(
    pathway = pathway,
    gene = genes,
    stringsAsFactors = FALSE
  )
  gene_table <- rbind(gene_table, temp)
}

write.csv(gene_table, "figures/TableS1_gene_list.csv", row.names = FALSE)

# Table S2: Survival analysis results
write.csv(survival_results, "figures/TableS2_survival_results.csv", row.names = FALSE)

# Table S3: Pathway statistics
pathway_stats <- pathway_info %>%
  select(pathway, chinese_name, year_discovered, total_genes, core_genes, description)

write.csv(pathway_stats, "figures/TableS3_pathway_info.csv", row.names = FALSE)

cat("Supplementary tables saved!\n\n")

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("              Figure Generation Complete!                       \n")
cat("================================================================\n")
cat("\n")
cat("Generated files in 'figures/' directory:\n")
cat("\n")
cat("Main Figures:\n")
cat("  Figure 2: Gene Set Statistics\n")
cat("    - Figure2AC_genesets.pdf\n")
cat("    - Figure2B_upset.pdf\n")
cat("\n")
cat("  Figure 3: Pan-Cancer Analysis\n")
cat("    - Figure3A_pancancer_heatmap.pdf\n")
cat("    - Figure3B_tumor_normal_FC.pdf\n")
cat("    - Figure3C_boxplots.pdf\n")
cat("    - Figure3D_summary.pdf\n")
cat("\n")
cat("  Figure 4: Survival Analysis\n")
cat("    - Figure4A_forest_ferroptosis.pdf\n")
cat("    - Figure4B_forest_multipathway.pdf\n")
cat("    - Figure4C_KM_curves.pdf\n")
cat("    - Figure4D_survival_heatmap.pdf\n")
cat("\n")
cat("  Figure 5: Pathway Correlation\n")
cat("    - Figure5A_correlation_heatmap.pdf\n")
cat("    - Figure5B_correlation_filtered.pdf\n")
cat("    - Figure5C_scatter_plots.pdf\n")
cat("\n")
cat("  Figure 6: Package Demo\n")
cat("    - Figure6A_radar.pdf\n")
cat("    - Figure6B_boxplot.pdf\n")
cat("    - Figure6C_distribution.pdf\n")
cat("    - Figure6D_classification.pdf\n")
cat("\n")
cat("  Figure 7: Enrichment Analysis\n")
cat("    - Figure7A_ORA.pdf\n")
cat("    - Figure7B_dotplot.pdf\n")
cat("    - Figure7C_volcano.pdf\n")
cat("\n")
cat("Supplementary Tables:\n")
cat("    - TableS1_gene_list.csv\n")
cat("    - TableS2_survival_results.csv\n")
cat("    - TableS3_pathway_info.csv\n")
cat("\n")
cat("NOTE: Figure 1 (Package Overview) should be created using\n")
cat("      BioRender, Figma, or PowerPoint.\n")
cat("\n")
cat("================================================================\n")
