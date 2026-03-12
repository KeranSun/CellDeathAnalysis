# =============================================================================
# 生成示例数据集
# CellDeathAnalysis R Package
# =============================================================================

set.seed(2024)

# -----------------------------------------------------------------------------
# 1. 获取所有细胞死亡相关基因
# -----------------------------------------------------------------------------

# 首先运行 build_data.R 确保基因集已创建
# source("data-raw/build_data.R")

# 收集所有基因
all_death_genes <- unique(c(
  # Ferroptosis
  "ACSL4", "LPCAT3", "ALOX5", "ALOX12", "ALOX15", "ALOX15B",

"NCOA4", "TFRC", "SLC39A14", "SLC39A8", "HMOX1", "HMOX2",
  "STEAP3", "FTH1", "FTL", "IREB2", "ACO1", "NOX1", "NOX4", "CYBB",
  "PTGS2", "POR", "FDFT1", "ATG5", "ATG7", "BECN1", "CHAC1", "CBS", "CTH",
  "SAT1", "GLS2", "TP53", "CDKN1A", "CARS1", "SLC1A5", "GOT1",
  "HSPB1", "CISD1", "CISD2", "VDAC2", "VDAC3",
  "GPX4", "SLC7A11", "SLC3A2", "FSP1", "DHODH", "GCH1",
  "NFE2L2", "KEAP1", "GCLC", "GCLM", "GSS", "GSR",
  "NQO1", "TXNRD1", "TXN", "SLC40A1", "HAMP",
  "PROM2", "MT1G", "MT2A", "HSPA5", "NFS1", "ISCU",
  "SQSTM1", "FANCD2", "CDKN2A", "ARNTL", "SCD", "FADS2", "ACSL3", "HILPDA",
  
 # Cuproptosis
  "FDX1", "LIAS", "LIPT1", "LIPT2", "DLD", "DLAT", 
  "PDHA1", "PDHB", "DLST", "GCSH", "DBT",
  "MTF1", "GLS",
  "SLC31A1", "SLC31A2", "ATP7A", "ATP7B", "ATOX1", "CCS", "COX17", "SCO1", "SCO2",
  
  # Disulfidptosis
  "NDUFS1", "NUBPL", "LRPPRC", "NDUFA11", "OXSM", "GYS1",
  "ACTN4", "MYH9", "MYH10", "FLNA", "FLNB", "TLN1", "IQGAP1",
  "RPN1", "SLC2A1", "SLC2A3", "HK1", "HK2",
  
  # Pyroptosis
  "NLRP1", "NLRP3", "NLRP6", "NLRP7", "NLRP12",
  "NLRC4", "NAIP", "AIM2", "PYCARD", "MEFV",
  "CASP1", "CASP4", "CASP5", "CASP3", "CASP8",
  "GSDMA", "GSDMB", "GSDMC", "GSDMD", "GSDME", "PJVK",
  "IL1A", "IL1B", "IL18", "IL33", "HMGB1", "S100A8", "S100A9",
  "TLR2", "TLR4", "TLR7", "TLR9", "NFKB1", "RELA",
  "NEK7", "P2RX7", "PANX1", "GBP1", "GBP2", "GBP5", "IRGM", "TXNIP", "BRCC3", "SIRT2",
  
  # Necroptosis
  "RIPK1", "RIPK3", "MLKL", "FADD", "TRADD", "CFLAR",
  "TNFRSF1A", "TNFRSF1B", "FAS", "TNFRSF10A", "TNFRSF10B",
  "TNF", "FASLG", "TNFSF10",
  "TLR3", "TICAM1", "ZBP1", "IFNAR1", "IFNAR2",
  "CYLD", "BIRC2", "BIRC3", "XIAP", "SHARPIN",
  "RNF31", "RBCK1", "OTULIN", "TNFAIP3", "SPATA2",
  "MAP3K7", "TAB1", "TAB2", "TAB3", "IKBKB", "CHUK", "IKBKG",
  "MAPK8", "MAPK9", "MAPK14",
  
  # Apoptosis
  "TNFRSF10A", "TNFRSF10B", "CASP10",
  "BAX", "BAK1", "BOK", "BID", "BCL2L11", "BAD", "BIK", "BMF", "HRK", 
  "PMAIP1", "BBC3", "BCL2", "BCL2L1", "BCL2L2", "MCL1", "BCL2A1",
  "CYCS", "DIABLO", "HTRA2", "ENDOG", "AIFM1", "APAF1", "CASP9",
  "CASP6", "CASP7",
  "BIRC5", "BIRC6", "BIRC7",
  "MDM2", "GADD45A", "GADD45B", "GADD45G", "SESN1", "SESN2", "PIDD1", "CASP2",
  "DDIT3", "ERN1", "EIF2AK3", "ATF4", "ATF6", "CASP12",
  
  # Autophagy
  "ULK1", "ULK2", "ATG13", "RB1CC1", "ATG101",
  "PIK3C3", "PIK3R4", "ATG14", "AMBRA1", "UVRAG", "SH3GLB1", "RUBCN",
  "ATG9A", "ATG9B", "WIPI1", "WIPI2", "ATG2A", "ATG2B", "WDR45", "WDR45B",
  "ATG10", "ATG12", "ATG16L1", "ATG16L2",
  "ATG3", "ATG4A", "ATG4B", "ATG4C", "ATG4D",
  "MAP1LC3A", "MAP1LC3B", "MAP1LC3B2", "MAP1LC3C",
  "GABARAP", "GABARAPL1", "GABARAPL2",
  "NBR1", "CALCOCO2", "OPTN", "TAX1BP1", "TOLLIP",
  "BNIP3", "BNIP3L", "FUNDC1", "PHB2", "PINK1", "PRKN", "FAM134B", "SEC62",
  "MTOR", "RPTOR", "MLST8", "AKT1S1", "TSC1", "TSC2", "RHEB", "RPS6KB1", "EIF4EBP1",
  "PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2", "PRKAG1", "PRKAG2", "PRKAG3",
  "STK11", "STRADA", "CAB39",
  "DRAM1", "DAPK1", "TFEB", "TFE3", "FOXO1", "FOXO3", "SIRT1", "HDAC6", "VCP",
  "RAB7A", "LAMP1", "LAMP2", "STX17", "SNAP29", "VAMP8",
  "PLEKHM1", "EPG5", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39", "VPS41",
  
  # PANoptosis
  "IRF1",
  
  # NETosis
  "MPO", "ELANE", "CTSG", "PRTN3", "PADI4", "PADI2",
  "NCF1", "NCF2", "NCF4", "CYBA", "RAC1", "RAC2",
  "H2AFX", "H2AFY", "H3F3A", "H3F3B",
  "PRKCA", "PRKCB", "PRKCD", "RAF1", "MAPK1", "MAPK3", "SYK", "AKT1",
  "CAMP", "S100A12", "DEFA1", "DEFA3", "DEFA4", "LCN2", "MMP9",
  
  # Parthanatos
  "PARP1", "MIF", "PARG", "ADPRHL2", "RNLS",
  "ATM", "ATR", "XRCC1", "PARP2",
  "NAMPT", "NMNAT1", "NMNAT2", "NMNAT3",
  
  # Entosis
  "RHOA", "ROCK1", "ROCK2", "CDC42",
  "ACTB", "ACTG1", "MYL9", "MYL12A", "MYL12B",
  "CDH1", "CTNNB1", "CTNNA1", "JUP",
  "CTSL", "CTSD",
  
  # Oxeiptosis
  "PGAM5",
  "SOD1", "SOD2", "CAT", "GPX1",
  
  # Alkaliptosis
  "CA9", "SLC4A7", "SLC9A1", "ATP6V1A", "ATP6V0D1",
  
  # LDCD
  "CTSB", "CTSH", "LGMN", "ACP2"
))

# 去重
all_death_genes <- unique(all_death_genes)
cat("Total unique genes:", length(all_death_genes), "\n")

# 添加一些背景基因（非细胞死亡相关）
background_genes <- c(
  "GAPDH", "ACTB", "B2M", "HPRT1", "RPL13A", "SDHA", "TBP", "YWHAZ",
  "PPIA", "RPLP0", "GUSB", "HMBS", "HSP90AB1", "PGK1", "TFRC",
  "ALDOA", "ENO1", "PKM", "LDHA", "LDHB", "TPI1", "GPI", "PGM1",
  "PFKL", "PFKM", "PFKP", "HK3", "G6PD", "PGD", "TALDO1", "TKT",
  "EEF1A1", "EEF2", "RPS3", "RPS4X", "RPS6", "RPS18", "RPL3", "RPL4",
  "MYC", "JUN", "FOS", "EGR1", "CREB1", "SP1", "YY1", "E2F1",
  "CDK1", "CDK2", "CDK4", "CDK6", "CCNA2", "CCNB1", "CCND1", "CCNE1",
  "RB1", "CDKN1B", "CDKN2B", "CDC25A", "CDC25B", "CDC25C",
  "EGFR", "ERBB2", "KRAS", "BRAF", "PIK3CA", "PTEN", "AKT2", "AKT3",
  "STAT1", "STAT3", "STAT5A", "STAT5B", "JAK1", "JAK2",
  "VEGFA", "VEGFB", "VEGFC", "FLT1", "KDR", "FLT4",
  "MMP2", "MMP3", "MMP7", "MMP14", "TIMP1", "TIMP2", "TIMP3",
  "CD4", "CD8A", "CD8B", "CD19", "CD20", "CD3D", "CD3E", "CD3G",
  "FOXP3", "CTLA4", "PDCD1", "CD274", "LAG3", "HAVCR2", "TIGIT",
  "IL2", "IL4", "IL6", "IL10", "IL12A", "IL12B", "IL17A", "IFNG", "TGFB1",
  "COL1A1", "COL1A2", "COL3A1", "FN1", "VIM", "ACTA2", "DES", "KRT18", "KRT19"
)

all_genes <- unique(c(all_death_genes, background_genes))
cat("Total genes in dataset:", length(all_genes), "\n")

# -----------------------------------------------------------------------------
# 2. 生成模拟表达数据
# -----------------------------------------------------------------------------

n_samples <- 200  # 100 normal + 100 tumor
n_genes <- length(all_genes)

# 创建基础表达矩阵 (模拟log2 TPM)
base_expr <- matrix(
  rnorm(n_genes * n_samples, mean = 6, sd = 2),
  nrow = n_genes,
  ncol = n_samples
)

rownames(base_expr) <- all_genes
colnames(base_expr) <- paste0("Sample_", sprintf("%03d", 1:n_samples))

# 定义样本组
sample_groups <- c(rep("Normal", 100), rep("Tumor", 100))

# 为肿瘤样本添加差异表达模式
tumor_idx <- 101:200
normal_idx <- 1:100

# 铁死亡相关基因在肿瘤中上调
ferro_drivers <- c("ACSL4", "LPCAT3", "TFRC", "NCOA4", "PTGS2", "CHAC1", "NOX1", "NOX4")
ferro_suppressors <- c("GPX4", "SLC7A11", "FSP1", "DHODH", "GCH1", "NFE2L2")

for (gene in ferro_drivers) {
  if (gene %in% rownames(base_expr)) {
    base_expr[gene, tumor_idx] <- base_expr[gene, tumor_idx] + rnorm(100, mean = 1.5, sd = 0.5)
  }
}

for (gene in ferro_suppressors) {
  if (gene %in% rownames(base_expr)) {
    base_expr[gene, tumor_idx] <- base_expr[gene, tumor_idx] - rnorm(100, mean = 0.8, sd = 0.3)
  }
}

# 焦亡相关基因
pyro_genes <- c("NLRP3", "CASP1", "GSDMD", "IL1B", "IL18", "PYCARD")
for (gene in pyro_genes) {
  if (gene %in% rownames(base_expr)) {
    base_expr[gene, tumor_idx] <- base_expr[gene, tumor_idx] + rnorm(100, mean = 1.2, sd = 0.4)
  }
}

# 铜死亡相关基因
cupro_genes <- c("FDX1", "LIAS", "LIPT1", "DLAT", "DLD")
for (gene in cupro_genes) {
  if (gene %in% rownames(base_expr)) {
    base_expr[gene, tumor_idx] <- base_expr[gene, tumor_idx] + rnorm(100, mean = 0.8, sd = 0.3)
  }
}

# 凋亡相关
apop_pro <- c("BAX", "BAK1", "CASP3", "CASP9", "CYCS")
apop_anti <- c("BCL2", "BCL2L1", "MCL1", "BIRC5")

for (gene in apop_pro) {
  if (gene %in% rownames(base_expr)) {
    base_expr[gene, tumor_idx] <- base_expr[gene, tumor_idx] + rnorm(100, mean = 0.6, sd = 0.3)
  }
}

for (gene in apop_anti) {
  if (gene %in% rownames(base_expr)) {
    base_expr[gene, tumor_idx] <- base_expr[gene, tumor_idx] + rnorm(100, mean = 1.0, sd = 0.4)
  }
}

# 自噬相关
autophagy_genes <- c("BECN1", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1", "ULK1")
for (gene in autophagy_genes) {
  if (gene %in% rownames(base_expr)) {
    base_expr[gene, tumor_idx] <- base_expr[gene, tumor_idx] + rnorm(100, mean = 0.7, sd = 0.3)
  }
}

# 确保表达值为正
base_expr[base_expr < 0] <- 0.01

# 添加一些噪声使数据更真实
base_expr <- base_expr + matrix(rnorm(n_genes * n_samples, mean = 0, sd = 0.3), 
                                 nrow = n_genes, ncol = n_samples)
base_expr[base_expr < 0] <- 0.01

# 创建示例表达矩阵
example_expr <- base_expr

# -----------------------------------------------------------------------------
# 3. 生成临床数据
# -----------------------------------------------------------------------------

# 根据细胞死亡得分模拟生存数据
# 高铁死亡/焦亡得分 -> 较差预后

# 计算简单的死亡评分用于生成生存数据
ferro_score <- colMeans(base_expr[ferro_drivers[ferro_drivers %in% rownames(base_expr)], ])
pyro_score <- colMeans(base_expr[pyro_genes[pyro_genes %in% rownames(base_expr)], ])

# 综合得分
death_score <- scale(ferro_score) + scale(pyro_score)

# 生成生存时间 (与死亡得分负相关)
base_survival <- 60  # 基础生存时间（月）
survival_time <- base_survival - death_score * 10 + rnorm(n_samples, mean = 0, sd = 15)
survival_time[survival_time < 1] <- runif(sum(survival_time < 1), 1, 6)
survival_time <- round(survival_time, 1)

# 生成生存状态 (得分高的更容易死亡)
death_prob <- pnorm(scale(death_score))
survival_status <- rbinom(n_samples, 1, prob = death_prob * 0.7 + 0.1)

# 生成其他临床特征
example_clinical <- data.frame(
  sample_id = colnames(base_expr),
  group = sample_groups,
  
  # 生存数据
  OS_time = survival_time,
  OS_status = survival_status,
  
  # 临床分期 (肿瘤样本)
  stage = c(rep(NA, 100), sample(c("I", "II", "III", "IV"), 100, 
                                  replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))),
  
  # 分级
  grade = c(rep(NA, 100), sample(c("G1", "G2", "G3"), 100, 
                                  replace = TRUE, prob = c(0.25, 0.45, 0.3))),
  
  # 年龄
  age = round(c(rnorm(100, mean = 45, sd = 12), rnorm(100, mean = 58, sd = 10))),
  
  # 性别
  gender = sample(c("Male", "Female"), n_samples, replace = TRUE, prob = c(0.55, 0.45)),
  
  # 治疗响应 (仅肿瘤)
  treatment_response = c(rep(NA, 100), sample(c("CR", "PR", "SD", "PD"), 100,
                                               replace = TRUE, prob = c(0.15, 0.35, 0.30, 0.20))),
  
  stringsAsFactors = FALSE
)

# 调整年龄范围
example_clinical$age[example_clinical$age < 20] <- 20
example_clinical$age[example_clinical$age > 85] <- 85

# -----------------------------------------------------------------------------
# 4. 保存数据
# -----------------------------------------------------------------------------

# 保存表达矩阵
usethis::use_data(example_expr, overwrite = TRUE)

# 保存临床数据
usethis::use_data(example_clinical, overwrite = TRUE)

cat("\n========================================\n")
cat("示例数据已成功创建！\n")
cat("========================================\n\n")
cat("数据集信息:\n")
cat("  - example_expr: ", nrow(example_expr), " genes x ", ncol(example_expr), " samples\n", sep = "")
cat("  - example_clinical: ", nrow(example_clinical), " samples x ", ncol(example_clinical), " variables\n", sep = "")
cat("\n样本分布:\n")
print(table(example_clinical$group))
cat("\n生存状态:\n")
print(table(example_clinical$OS_status))
