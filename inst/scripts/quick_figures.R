# =============================================================================
# CellDeathAnalysis - Quick Figure Generation (Simplified)
# =============================================================================
#
# Run this script to quickly generate all publication figures
# Uses built-in example data - no downloads required!
#
# Author: Keran Sun
# =============================================================================

# -----------------------------------------------------------------------------
# Setup - Install packages if needed
# -----------------------------------------------------------------------------

cat("Installing/loading required packages...\n")

pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork", "pheatmap", 
          "corrplot", "survival", "survminer", "ggpubr", "RColorBrewer")

for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

library(CellDeathAnalysis)

# Create output folder
dir.create("pub_figures", showWarnings = FALSE)

# Set high-quality output
options(bitmapType = "cairo")

cat("Setup complete!\n\n")

# -----------------------------------------------------------------------------
# Load Data and Calculate Scores
# -----------------------------------------------------------------------------

cat("Loading data and calculating scores...\n")

data(example_expr)
data(example_clinical)

scores <- calculate_death_score(example_expr, method = "zscore", verbose = FALSE)
group <- example_clinical$group

cat("  Samples:", nrow(scores), "\n")
cat("  Pathways:", ncol(scores), "\n\n")

# Color settings
tumor_col <- "#D62728"
normal_col <- "#2CA02C"
group_cols <- c("Normal" = normal_col, "Tumor" = tumor_col)

# =============================================================================
# FIGURE 2: Gene Set Statistics
# =============================================================================

cat("Creating Figure 2: Gene Set Statistics...\n")

# Get pathway info
pinfo <- list_death_pathways(detailed = TRUE)
pinfo <- pinfo[order(pinfo$total_genes, decreasing = TRUE), ]
pinfo$pathway <- factor(pinfo$pathway, levels = rev(pinfo$pathway))

# Figure 2A: Gene counts
fig2a <- ggplot(pinfo, aes(x = pathway, y = total_genes)) +
  geom_col(fill = "#3498DB", color = "black", width = 0.7) +
  geom_text(aes(label = total_genes), hjust = -0.3, size = 3.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 12) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(x = NULL, y = "Number of Genes", 
       title = "A. Genes per Cell Death Pathway")

# Figure 2B: Overlap heatmap
gsets <- get_death_geneset("all", type = "simple")
n_path <- length(gsets)
overlap_mat <- matrix(0, n_path, n_path, 
                      dimnames = list(names(gsets), names(gsets)))

for (i in 1:n_path) {
  for (j in 1:n_path) {
    overlap_mat[i,j] <- length(intersect(gsets[[i]], gsets[[j]]))
  }
}

# Convert to long format
overlap_df <- reshape2::melt(overlap_mat)
names(overlap_df) <- c("Path1", "Path2", "Count")

fig2b <- ggplot(overlap_df, aes(x = Path1, y = Path2, fill = Count)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")), size = 2.5) +
  scale_fill_gradient(low = "white", high = "#E74C3C", name = "Overlap") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9)
  ) +
  labs(x = NULL, y = NULL, title = "B. Gene Overlap Matrix") +
  coord_fixed()

# Combine and save
fig2 <- fig2a + fig2b + plot_layout(widths = c(1, 1.3))
ggsave("pub_figures/Figure2_GeneSetStats.pdf", fig2, width = 14, height = 6)
ggsave("pub_figures/Figure2_GeneSetStats.png", fig2, width = 14, height = 6, dpi = 300)

cat("  Figure 2 saved!\n")

# =============================================================================
# FIGURE 3: Pan-cancer Style Heatmap and Comparison
# =============================================================================

cat("Creating Figure 3: Expression Analysis...\n")

# Figure 3A: Heatmap
score_mat <- as.matrix(scores)
score_scaled <- t(scale(t(score_mat)))  # Scale by sample

# Annotation
anno_df <- data.frame(Group = group, row.names = rownames(scores))

# Order by group
ord <- order(group)
score_ordered <- score_scaled[ord, ]
anno_ordered <- data.frame(Group = group[ord], row.names = rownames(score_ordered))

pdf("pub_figures/Figure3A_Heatmap.pdf", width = 12, height = 8)
pheatmap(
  t(score_ordered),  # Pathways as rows
  annotation_col = anno_ordered,
  annotation_colors = list(Group = group_cols),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-2, 2, length.out = 101),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize = 10,
  fontsize_row = 11,
  border_color = NA,
  main = "A. Cell Death Pathway Scores"
)
dev.off()

# Figure 3B: Boxplot comparison
score_long <- scores %>%
  mutate(Sample = row.names(.), Group = group) %>%
  tidyr::pivot_longer(cols = -c(Sample, Group), 
                      names_to = "Pathway", values_to = "Score")

fig3b <- ggplot(score_long, aes(x = Pathway, y = Score, fill = Group)) +
  geom_boxplot(outlier.size = 0.8, width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = group_cols) +
  stat_compare_means(aes(group = Group), method = "wilcox.test",
                     label = "p.signif", label.y.npc = 0.95, size = 3.5) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = "Pathway Score (Z-score)", 
       title = "B. Tumor vs Normal Comparison")

ggsave("pub_figures/Figure3B_Boxplot.pdf", fig3b, width = 12, height = 6)
ggsave("pub_figures/Figure3B_Boxplot.png", fig3b, width = 12, height = 6, dpi = 300)

# Figure 3C: Selected pathways violin
key_paths <- c("ferroptosis", "cuproptosis", "pyroptosis", "apoptosis")

fig3c <- ggplot(score_long %>% filter(Pathway %in% key_paths), 
                aes(x = Pathway, y = Score, fill = Group)) +
  geom_violin(position = position_dodge(0.8), alpha = 0.7, scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.15, 
               outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = group_cols) +
  stat_compare_means(aes(group = Group), method = "wilcox.test",
                     label = "p.format", vjust = -0.5, size = 4) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top") +
  labs(x = NULL, y = "Pathway Score", 
       title = "C. Key Pathways Violin Plot")

ggsave("pub_figures/Figure3C_Violin.pdf", fig3c, width = 10, height = 6)
ggsave("pub_figures/Figure3C_Violin.png", fig3c, width = 10, height = 6, dpi = 300)

cat("  Figure 3 saved!\n")

# =============================================================================
# FIGURE 4: Survival Analysis
# =============================================================================

cat("Creating Figure 4: Survival Analysis...\n")

# Tumor samples only
tumor_idx <- group == "Tumor"
tumor_scores <- scores[tumor_idx, ]
tumor_clin <- example_clinical[tumor_idx, ]

# Figure 4A: Forest plot
surv_res <- data.frame()

for (path in colnames(scores)) {
  tryCatch({
    s <- scale(tumor_scores[[path]])[,1]
    df <- data.frame(time = tumor_clin$OS_time, 
                     status = tumor_clin$OS_status, 
                     score = s)
    df <- na.omit(df)
    
    fit <- coxph(Surv(time, status) ~ score, data = df)
    sm <- summary(fit)
    
    surv_res <- rbind(surv_res, data.frame(
      Pathway = path,
      HR = sm$conf.int[1,1],
      Lower = sm$conf.int[1,3],
      Upper = sm$conf.int[1,4],
      P = sm$coefficients[1,5]
    ))
  }, error = function(e) {})
}

surv_res$P_adj <- p.adjust(surv_res$P, "BH")
surv_res$Sig <- ifelse(surv_res$P_adj < 0.05, "p < 0.05", "NS")
surv_res <- surv_res[order(surv_res$HR), ]
surv_res$Pathway <- factor(surv_res$Pathway, levels = surv_res$Pathway)

fig4a <- ggplot(surv_res, aes(x = HR, y = Pathway)) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray40", size = 0.8) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, color = Sig), 
                 height = 0.25, size = 0.9) +
  geom_point(aes(color = Sig, size = -log10(P_adj)), shape = 18) +
  scale_color_manual(values = c("p < 0.05" = "#E41A1C", "NS" = "gray50")) +
  scale_size_continuous(range = c(3, 8), name = "-log10(p)") +
  scale_x_log10() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Hazard Ratio (95% CI)", y = NULL,
       title = "A. Forest Plot: Pathways and Survival",
       color = "Significance")

ggsave("pub_figures/Figure4A_Forest.pdf", fig4a, width = 10, height = 8)
ggsave("pub_figures/Figure4A_Forest.png", fig4a, width = 10, height = 8, dpi = 300)

# Figure 4B: KM curves
km_list <- list()

for (path in key_paths) {
  s <- tumor_scores[[path]]
  df <- data.frame(
    time = tumor_clin$OS_time,
    status = tumor_clin$OS_status,
    group = factor(ifelse(s > median(s), "High", "Low"), 
                   levels = c("Low", "High"))
  )
  df <- na.omit(df)
  
  fit <- survfit(Surv(time, status) ~ group, data = df)
  
  km_list[[path]] <- ggsurvplot(
    fit, data = df,
    pval = TRUE, pval.size = 4,
    conf.int = TRUE, conf.int.alpha = 0.15,
    risk.table = FALSE,
    palette = c("#3498DB", "#E74C3C"),
    legend.title = "Score",
    legend.labs = c("Low", "High"),
    xlab = "Time", ylab = "Survival",
    title = tools::toTitleCase(path),
    ggtheme = theme_bw(base_size = 11)
  )
}

pdf("pub_figures/Figure4B_KM_Curves.pdf", width = 12, height = 10)
arrange_ggsurvplots(km_list, ncol = 2, nrow = 2)
dev.off()

png("pub_figures/Figure4B_KM_Curves.png", width = 12, height = 10, 
    units = "in", res = 300)
arrange_ggsurvplots(km_list, ncol = 2, nrow = 2)
dev.off()

cat("  Figure 4 saved!\n")

# =============================================================================
# FIGURE 5: Correlation Analysis
# =============================================================================

cat("Creating Figure 5: Pathway Correlation...\n")

cor_mat <- cor(scores, method = "spearman")

# PDF version
pdf("pub_figures/Figure5_Correlation.pdf", width = 10, height = 9)
corrplot(
  cor_mat,
  method = "color",
  type = "lower",
  order = "hclust",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.9,
  addCoef.col = "black",
  number.cex = 0.7,
  col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  title = "Cell Death Pathway Correlation (Spearman)",
  mar = c(0,0,2,0)
)
dev.off()

# PNG version with ggplot
cor_df <- reshape2::melt(cor_mat)
names(cor_df) <- c("Path1", "Path2", "Cor")

fig5 <- ggplot(cor_df, aes(x = Path1, y = Path2, fill = Cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Cor)), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Pathway Correlation (Spearman)") +
  coord_fixed()

ggsave("pub_figures/Figure5_Correlation.png", fig5, width = 10, height = 9, dpi = 300)

cat("  Figure 5 saved!\n")

# =============================================================================
# FIGURE 6: Distribution Analysis
# =============================================================================

cat("Creating Figure 6: Score Distribution...\n")

# Density plots
fig6a <- ggplot(score_long %>% filter(Pathway %in% key_paths),
                aes(x = Score, fill = Group, color = Group)) +
  geom_density(alpha = 0.4, size = 0.8) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "#F0F0F0"),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(x = "Score", y = "Density", title = "A. Score Distribution by Group")

# Ridge plot style
fig6b <- ggplot(score_long, aes(x = Score, y = Pathway, fill = Group)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = group_cols) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top") +
  labs(x = "Pathway Score", y = NULL, title = "B. All Pathways Comparison")

fig6 <- fig6a / fig6b + plot_layout(heights = c(1.2, 1))
ggsave("pub_figures/Figure6_Distribution.pdf", fig6, width = 11, height = 12)
ggsave("pub_figures/Figure6_Distribution.png", fig6, width = 11, height = 12, dpi = 300)

cat("  Figure 6 saved!\n")

# =============================================================================
# FIGURE 7: Enrichment Analysis
# =============================================================================

cat("Creating Figure 7: Enrichment Analysis...\n")

# Create example gene list
set.seed(42)
death_genes <- get_all_death_genes()
test_genes <- c(
  sample(rownames(example_expr)[rownames(example_expr) %in% death_genes], 25),
  sample(rownames(example_expr), 75)
)

ora <- death_enrich_ora(unique(test_genes), min_genes = 2)

if (!is.null(ora) && nrow(ora) > 0) {
  ora <- ora[order(ora$p_value), ]
  ora$pathway <- factor(ora$pathway, levels = rev(ora$pathway))
  
  fig7a <- ggplot(ora, aes(x = -log10(p_adjust), y = pathway)) +
    geom_col(aes(fill = fold_enrichment), color = "black", width = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = 2, color = "red") +
    scale_fill_gradient(low = "#FDDBC7", high = "#B2182B", name = "Fold\nEnrichment") +
    theme_bw(base_size = 12) +
    labs(x = "-log10(Adjusted P-value)", y = NULL,
         title = "A. Over-Representation Analysis")
  
  fig7b <- ggplot(ora, aes(x = fold_enrichment, y = pathway)) +
    geom_point(aes(size = n_overlap, color = -log10(p_adjust))) +
    scale_color_gradient(low = "#4575B4", high = "#D73027", 
                         name = "-log10(p)") +
    scale_size_continuous(range = c(4, 12), name = "Gene\nCount") +
    theme_bw(base_size = 12) +
    labs(x = "Fold Enrichment", y = NULL, title = "B. Enrichment Dot Plot")
  
  fig7 <- fig7a + fig7b
  ggsave("pub_figures/Figure7_Enrichment.pdf", fig7, width = 14, height = 6)
  ggsave("pub_figures/Figure7_Enrichment.png", fig7, width = 14, height = 6, dpi = 300)
}

cat("  Figure 7 saved!\n")

# =============================================================================
# Summary Tables
# =============================================================================

cat("\nCreating summary tables...\n")

# Table 1: Gene set summary
table1 <- pinfo %>% 
  select(pathway, total_genes, year_discovered) %>%
  arrange(desc(total_genes))
write.csv(table1, "pub_figures/Table1_GeneSetSummary.csv", row.names = FALSE)

# Table 2: Survival results
write.csv(surv_res, "pub_figures/Table2_SurvivalResults.csv", row.names = FALSE)

# Table 3: Group comparison statistics
stats_table <- score_long %>%
  group_by(Pathway) %>%
  summarise(
    Mean_Normal = mean(Score[Group == "Normal"]),
    Mean_Tumor = mean(Score[Group == "Tumor"]),
    P_value = wilcox.test(Score ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(P_adj = p.adjust(P_value, "BH"))

write.csv(stats_table, "pub_figures/Table3_GroupComparison.csv", row.names = FALSE)

# =============================================================================
# Final Output
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║           ALL FIGURES GENERATED SUCCESSFULLY!                ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat("║                                                              ║\n")
cat("║  Output folder: ./pub_figures/                               ║\n")
cat("║                                                              ║\n")
cat("║  Main Figures:                                               ║\n")
cat("║    • Figure2_GeneSetStats.pdf/png                            ║\n")
cat("║    • Figure3A_Heatmap.pdf                                    ║\n")
cat("║    • Figure3B_Boxplot.pdf/png                                ║\n")
cat("║    • Figure3C_Violin.pdf/png                                 ║\n")
cat("║    • Figure4A_Forest.pdf/png                                 ║\n")
cat("║    • Figure4B_KM_Curves.pdf/png                              ║\n")
cat("║    • Figure5_Correlation.pdf/png                             ║\n")
cat("║    • Figure6_Distribution.pdf/png                            ║\n")
cat("║    • Figure7_Enrichment.pdf/png                              ║\n")
cat("║                                                              ║\n")
cat("║  Tables:                                                     ║\n")
cat("║    • Table1_GeneSetSummary.csv                               ║\n")
cat("║    • Table2_SurvivalResults.csv                              ║\n")
cat("║    • Table3_GroupComparison.csv                              ║\n")
cat("║                                                              ║\n")
cat("║  Note: Figure 1 (workflow) should be created in              ║\n")
cat("║        BioRender or similar software                         ║\n")
cat("║                                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")

# List files
cat("\nGenerated files:\n")
print(list.files("pub_figures", full.names = FALSE))
