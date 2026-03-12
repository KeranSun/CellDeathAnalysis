#!/usr/bin/env Rscript
# =============================================================================
# 生成 CellDeathAnalysis 包的数据文件
# 运行方法: source("data-raw/build_data.R")
# =============================================================================

# 设置工作目录（如果需要）
# setwd("path/to/CellDeathAnalysis")

# -----------------------------------------------------------------------------
# 1. 创建基因集数据
# -----------------------------------------------------------------------------

# Ferroptosis
ferroptosis <- list(
  driver = c(
    "ACSL4", "LPCAT3", "ALOX5", "ALOX12", "ALOX15", "ALOX15B",
    "NCOA4", "TFRC", "SLC39A14", "SLC39A8", "HMOX1", "HMOX2",
    "STEAP3", "FTH1", "FTL", "IREB2", "ACO1",
    "NOX1", "NOX4", "CYBB",
    "PTGS2", "POR", "FDFT1",
    "ATG5", "ATG7", "BECN1",
    "CHAC1", "CBS", "CTH",
    "SAT1", "GLS2",
    "TP53", "CDKN1A",
    "CARS1", "SLC1A5", "GOT1",
    "HSPB1", "CISD1", "CISD2",
    "VDAC2", "VDAC3"
  ),
  suppressor = c(
    "GPX4", "SLC7A11", "SLC3A2",
    "FSP1", "DHODH", "GCH1",
    "NFE2L2", "KEAP1",
    "GCLC", "GCLM", "GSS", "GSR",
    "NQO1", "TXNRD1", "TXN",
    "SLC40A1", "HAMP",
    "PROM2", "MT1G", "MT2A",
    "HSPA5",
    "NFS1", "ISCU",
    "SQSTM1", "FANCD2",
    "CDKN2A", "ARNTL",
    "SCD", "FADS2",
    "ACSL3", "HILPDA"
  ),
  marker = c("PTGS2", "CHAC1", "ACSL4", "TFRC", "SLC7A11", "GPX4", 
             "HMOX1", "FTH1", "FTL", "NCOA4")
)

# Cuproptosis
cuproptosis <- list(
  positive = c("FDX1", "LIAS", "LIPT1", "LIPT2", "DLD", "DLAT", 
               "PDHA1", "PDHB", "DLST", "GCSH", "DBT"),
  negative = c("MTF1", "GLS", "CDKN2A"),
  transport = c("SLC31A1", "SLC31A2", "ATP7A", "ATP7B", "ATOX1", 
                "CCS", "COX17", "SCO1", "SCO2")
)

# Disulfidptosis
disulfidptosis <- list(
  core = c("SLC7A11", "SLC3A2", "NDUFS1", "NUBPL", "LRPPRC", 
           "NDUFA11", "OXSM", "GYS1"),
  cytoskeleton = c("ACTN4", "MYH9", "MYH10", "FLNA", "FLNB", "TLN1", "IQGAP1"),
  metabolism = c("RPN1", "SLC2A1", "SLC2A3", "HK1", "HK2")
)

# Pyroptosis
pyroptosis <- list(
  inflammasome = c("NLRP1", "NLRP3", "NLRP6", "NLRP7", "NLRP12",
                   "NLRC4", "NAIP", "AIM2", "PYCARD", "MEFV"),
  caspases = c("CASP1", "CASP4", "CASP5", "CASP3", "CASP8"),
  gasdermins = c("GSDMA", "GSDMB", "GSDMC", "GSDMD", "GSDME", "PJVK"),
  cytokines = c("IL1A", "IL1B", "IL18", "IL33", "HMGB1", "S100A8", "S100A9"),
  regulators = c("TLR2", "TLR4", "TLR7", "TLR9", "NFKB1", "RELA",
                 "NEK7", "P2RX7", "PANX1", "GBP1", "GBP2", "GBP5",
                 "IRGM", "TXNIP", "BRCC3", "SIRT2")
)

# Necroptosis
necroptosis <- list(
  core = c("RIPK1", "RIPK3", "MLKL", "FADD", "TRADD", "CASP8", "CFLAR"),
  receptors = c("TNFRSF1A", "TNFRSF1B", "FAS", "TNFRSF10A", "TNFRSF10B"),
  ligands = c("TNF", "FASLG", "TNFSF10"),
  prr = c("TLR3", "TLR4", "TICAM1", "ZBP1", "IFNAR1", "IFNAR2"),
  ubiquitination = c("CYLD", "BIRC2", "BIRC3", "XIAP", "SHARPIN",
                     "RNF31", "RBCK1", "OTULIN", "TNFAIP3", "SPATA2"),
  kinases = c("MAP3K7", "TAB1", "TAB2", "TAB3", "IKBKB", "CHUK", "IKBKG",
              "MAPK8", "MAPK9", "MAPK14")
)

# Apoptosis
apoptosis <- list(
  extrinsic = c("FAS", "TNFRSF1A", "TNFRSF10A", "TNFRSF10B",
                "FASLG", "TNF", "TNFSF10",
                "FADD", "TRADD", "CASP8", "CASP10"),
  intrinsic = c("BAX", "BAK1", "BOK",
                "BID", "BCL2L11", "BAD", "BIK", "BMF", "HRK", "PMAIP1", "BBC3",
                "BCL2", "BCL2L1", "BCL2L2", "MCL1", "BCL2A1",
                "CYCS", "DIABLO", "HTRA2", "ENDOG", "AIFM1",
                "APAF1", "CASP9"),
  executioner = c("CASP3", "CASP6", "CASP7"),
  iap = c("XIAP", "BIRC2", "BIRC3", "BIRC5", "BIRC6", "BIRC7", "NAIP"),
  p53_pathway = c("TP53", "MDM2", "CDKN1A", "GADD45A", "GADD45B", "GADD45G",
                  "SESN1", "SESN2", "PIDD1", "CASP2"),
  er_stress = c("DDIT3", "ERN1", "EIF2AK3", "ATF4", "ATF6", "CASP12")
)

# Autophagy
autophagy <- list(
  ulk_complex = c("ULK1", "ULK2", "ATG13", "RB1CC1", "ATG101"),
  pi3k_complex = c("BECN1", "PIK3C3", "PIK3R4", "ATG14", "AMBRA1", 
                   "UVRAG", "SH3GLB1", "RUBCN"),
  atg9_system = c("ATG9A", "ATG9B", "WIPI1", "WIPI2", "ATG2A", "ATG2B",
                  "WDR45", "WDR45B"),
  atg12_system = c("ATG5", "ATG7", "ATG10", "ATG12", "ATG16L1", "ATG16L2"),
  lc3_system = c("ATG3", "ATG4A", "ATG4B", "ATG4C", "ATG4D", "ATG7"),
  lc3_family = c("MAP1LC3A", "MAP1LC3B", "MAP1LC3B2", "MAP1LC3C",
                 "GABARAP", "GABARAPL1", "GABARAPL2"),
  receptors = c("SQSTM1", "NBR1", "CALCOCO2", "OPTN", "TAX1BP1", "TOLLIP",
                "BNIP3", "BNIP3L", "FUNDC1", "PHB2", "PINK1", "PRKN",
                "NCOA4", "FAM134B", "SEC62"),
  mtor_pathway = c("MTOR", "RPTOR", "MLST8", "AKT1S1",
                   "TSC1", "TSC2", "RHEB", "RPS6KB1", "EIF4EBP1"),
  ampk_pathway = c("PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2",
                   "PRKAG1", "PRKAG2", "PRKAG3", "STK11", "STRADA", "CAB39"),
  regulators = c("TP53", "DRAM1", "DAPK1", "TFEB", "TFE3",
                 "FOXO1", "FOXO3", "SIRT1", "HDAC6", "VCP"),
  lysosome_fusion = c("RAB7A", "LAMP1", "LAMP2", "STX17", "SNAP29", "VAMP8",
                      "PLEKHM1", "EPG5", "VPS11", "VPS16", "VPS18", 
                      "VPS33A", "VPS39", "VPS41")
)

# PANoptosis
panoptosis <- list(
  panoptosome = c("ZBP1", "RIPK1", "RIPK3", "CASP8", "FADD",
                  "NLRP3", "PYCARD", "CASP1", "CASP6"),
  effectors = c("MLKL", "GSDMD", "GSDME", "CASP3", "CASP7"),
  crosstalk = c("BID", "BAX", "BAK1", "CYCS", "AIFM1", "CFLAR", "TNFAIP3"),
  upstream = c("TNF", "TNFRSF1A", "IFNG", "IFNAR1", "TLR3", "TLR4", "IRF1")
)

# NETosis
netosis <- list(
  core = c("MPO", "ELANE", "CTSG", "PRTN3", "PADI4", "PADI2"),
  ros_generation = c("CYBB", "NCF1", "NCF2", "NCF4", "CYBA", "RAC1", "RAC2"),
  chromatin = c("H2AFX", "H2AFY", "H3F3A", "H3F3B"),
  signaling = c("PRKCA", "PRKCB", "PRKCD", "RAF1", "MAPK1", "MAPK3", "SYK", "AKT1"),
  released_factors = c("HMGB1", "CAMP", "S100A8", "S100A9", "S100A12",
                       "DEFA1", "DEFA3", "DEFA4", "LCN2", "MMP9")
)

# Parthanatos
parthanatos <- list(
  core = c("PARP1", "AIFM1", "MIF", "PARG", "ADPRHL2", "RNLS"),
  dna_damage = c("ATM", "ATR", "H2AFX", "XRCC1", "PARP2"),
  nad_metabolism = c("NAMPT", "NMNAT1", "NMNAT2", "NMNAT3", "SIRT1", "SIRT2")
)

# Entosis
entosis <- list(
  core = c("RHOA", "ROCK1", "ROCK2", "CDC42", "RAC1",
           "ACTB", "ACTG1", "MYH9", "MYH10", "MYL9", "MYL12A", "MYL12B"),
  adhesion = c("CDH1", "CTNNB1", "CTNNA1", "JUP"),
  degradation = c("LAMP1", "LAMP2", "CTSL", "CTSD", "RAB7A")
)

# Oxeiptosis
oxeiptosis <- list(
  core = c("KEAP1", "PGAM5", "AIFM1"),
  related = c("NFE2L2", "NQO1", "HMOX1", "SOD1", "SOD2", "CAT", "GPX1", "GPX4")
)

# Alkaliptosis
alkaliptosis <- list(
  core = c("CA9", "SLC4A7", "SLC9A1", "ATP6V1A", "ATP6V0D1"),
  related = c("NFKB1", "RELA")
)

# LDCD
ldcd <- list(
  core = c("LAMP1", "LAMP2", "CTSL", "CTSB", "CTSD", "CTSH", "LGMN", "ACP2"),
  regulators = c("TFEB", "TFE3", "HSPA1A", "HSPA1B", "BID", "BAX")
)

# 整合所有基因集
death_genesets <- list(
  ferroptosis = ferroptosis,
  cuproptosis = cuproptosis,
  disulfidptosis = disulfidptosis,
  pyroptosis = pyroptosis,
  necroptosis = necroptosis,
  apoptosis = apoptosis,
  autophagy = autophagy,
  panoptosis = panoptosis,
  netosis = netosis,
  parthanatos = parthanatos,
  entosis = entosis,
  oxeiptosis = oxeiptosis,
  alkaliptosis = alkaliptosis,
  ldcd = ldcd
)

# 添加属性
attr(death_genesets, "version") <- "1.0.0"
attr(death_genesets, "date") <- Sys.Date()
attr(death_genesets, "species") <- "Homo sapiens"

# -----------------------------------------------------------------------------
# 2. 创建简化版基因集
# -----------------------------------------------------------------------------

death_genesets_simple <- list(
  ferroptosis = unique(c(ferroptosis$driver[1:20], ferroptosis$suppressor[1:15], 
                         ferroptosis$marker)),
  cuproptosis = unique(c(cuproptosis$positive, cuproptosis$negative, 
                         cuproptosis$transport)),
  disulfidptosis = unique(c(disulfidptosis$core, disulfidptosis$cytoskeleton)),
  pyroptosis = unique(c(pyroptosis$inflammasome, pyroptosis$caspases, 
                        pyroptosis$gasdermins)),
  necroptosis = unique(c(necroptosis$core, necroptosis$receptors)),
  apoptosis = unique(c(apoptosis$extrinsic, apoptosis$intrinsic[1:15], 
                       apoptosis$executioner)),
  autophagy = unique(c(autophagy$ulk_complex, autophagy$pi3k_complex, 
                       autophagy$lc3_family, autophagy$receptors[1:8])),
  panoptosis = unique(c(panoptosis$panoptosome, panoptosis$effectors)),
  netosis = unique(c(netosis$core, netosis$ros_generation[1:5])),
  parthanatos = parthanatos$core,
  entosis = entosis$core,
  oxeiptosis = oxeiptosis$core
)

attr(death_genesets_simple, "version") <- "1.0.0"
attr(death_genesets_simple, "date") <- Sys.Date()
attr(death_genesets_simple, "species") <- "Homo sapiens"
attr(death_genesets_simple, "description") <- "Simplified gene sets with core genes"

# -----------------------------------------------------------------------------
# 3. 创建通路信息数据框
# -----------------------------------------------------------------------------

death_pathway_info <- data.frame(
  pathway = c("ferroptosis", "cuproptosis", "disulfidptosis", "pyroptosis",
              "necroptosis", "apoptosis", "autophagy", "panoptosis",
              "netosis", "parthanatos", "entosis", "oxeiptosis",
              "alkaliptosis", "ldcd"),
  
  full_name = c("Ferroptosis", "Cuproptosis", "Disulfidptosis", "Pyroptosis",
                "Necroptosis", "Apoptosis", "Autophagy", "PANoptosis",
                "NETosis", "Parthanatos", "Entosis", "Oxeiptosis",
                "Alkaliptosis", "Lysosome-dependent cell death"),
  
  chinese_name = c("铁死亡", "铜死亡", "二硫死亡", "焦亡",
                   "坏死性凋亡", "凋亡", "自噬", "泛凋亡",
                   "中性粒细胞胞外陷阱死亡", "PARP依赖性死亡", 
                   "细胞内吞死亡", "氧化诱导死亡",
                   "碱中毒死亡", "溶酶体依赖性死亡"),
  
  year_discovered = c(2012, 2022, 2023, 2001, 2005, 1972, 1963, 2019,
                      2004, 2009, 2007, 2018, 2018, 2000),
  
  key_paper = c(
    "Dixon et al., Cell 2012",
    "Tsvetkov et al., Science 2022",
    "Liu et al., Nat Cell Biol 2023",
    "Cookson & Brennan, Trends Mol Med 2001",
    "Degterev et al., Nat Chem Biol 2005",
    "Kerr et al., Br J Cancer 1972",
    "de Duve & Wattiaux, Annu Rev Physiol 1963",
    "Malireddi et al., PNAS 2019",
    "Brinkmann et al., Science 2004",
    "Andrabi et al., PNAS 2006",
    "Overholtzer et al., Cell 2007",
    "Holze et al., Nature 2018",
    "Song et al., PNAS 2018",
    "Multiple sources"
  ),
  
  description = c(
    "Iron-dependent lipid peroxidation leading to cell death",
    "Copper-dependent cell death through protein lipoylation",
    "Disulfide stress-induced cell death in glucose-starved cells",
    "Inflammatory cell death with gasdermin pore formation",
    "Programmed necrosis via RIPK1-RIPK3-MLKL pathway",
    "Classical programmed cell death with caspase activation",
    "Self-degradation via autophagosome-lysosome pathway",
    "Combined pyroptosis, apoptosis, and necroptosis",
    "Neutrophil death with extracellular trap release",
    "PARP1-dependent cell death with AIF translocation",
    "Cell-in-cell death through live cell engulfment",
    "ROS-induced death via KEAP1-PGAM5-AIFM1",
    "Cell death from intracellular alkalinization",
    "Cell death from lysosomal membrane permeabilization"
  ),
  
  key_genes = c(
    "GPX4, SLC7A11, ACSL4, TFRC, FSP1",
    "FDX1, LIAS, LIPT1, DLAT, DLD",
    "SLC7A11, NDUFS1, ACTN4, GYS1",
    "NLRP3, CASP1, GSDMD, IL1B, PYCARD",
    "RIPK1, RIPK3, MLKL, FADD, CASP8",
    "CASP3, CASP8, CASP9, BAX, BCL2",
    "BECN1, ATG5, ATG7, LC3B, SQSTM1",
    "ZBP1, RIPK3, CASP8, GSDMD, MLKL",
    "MPO, ELANE, PADI4, CYBB",
    "PARP1, AIFM1, MIF",
    "RHOA, ROCK1, MYH9, CDH1",
    "KEAP1, PGAM5, AIFM1",
    "CA9, SLC4A7, NHE1",
    "CTSL, CTSB, LAMP1, LAMP2"
  ),
  
  morphology = c(
    "Cell shrinkage, mitochondrial damage, lipid droplets",
    "Mitochondrial shrinkage, protein aggregation",
    "Cytoskeletal collapse, F-actin contraction",
    "Cell swelling, membrane rupture, cytokine release",
    "Cell swelling, membrane rupture, organelle swelling",
    "Cell shrinkage, chromatin condensation, apoptotic bodies",
    "Autophagic vacuoles, no inflammation",
    "Mixed features of three death modes",
    "NET release, chromatin decondensation",
    "DNA fragmentation, chromatin condensation",
    "Cell-in-cell structures",
    "Caspase-independent, no inflammation",
    "Cytoplasmic vacuolization",
    "Lysosomal rupture, acidification"
  ),
  
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 4. 保存数据文件
# -----------------------------------------------------------------------------

# 创建data目录
if (!dir.exists("data")) {
  dir.create("data")
}

# 保存数据
save(death_genesets, file = "data/death_genesets.rda", compress = "xz")
save(death_genesets_simple, file = "data/death_genesets_simple.rda", compress = "xz")
save(death_pathway_info, file = "data/death_pathway_info.rda", compress = "xz")

cat("\n========================================\n")
cat("数据文件已成功创建！\n")
cat("========================================\n\n")

cat("创建的文件:\n")
cat("  - data/death_genesets.rda\n")
cat("  - data/death_genesets_simple.rda\n")
cat("  - data/death_pathway_info.rda\n\n")

cat("基因集统计:\n")
for (name in names(death_genesets)) {
  n_genes <- length(unique(unlist(death_genesets[[name]])))
  cat(sprintf("  - %s: %d genes\n", name, n_genes))
}

cat("\n总计: ", length(unique(unlist(death_genesets))), " unique genes\n")

# =============================================================================
# 5. 生成示例数据
# =============================================================================

cat("\n正在生成示例数据...\n")

set.seed(2024)

# 收集所有基因
all_death_genes <- unique(unlist(death_genesets))

# 添加背景基因
background_genes <- c(
  "GAPDH", "ACTB", "B2M", "HPRT1", "RPL13A", "SDHA", "TBP", "YWHAZ",
  "PPIA", "RPLP0", "GUSB", "HMBS", "HSP90AB1", "PGK1",
  "ALDOA", "ENO1", "PKM", "LDHA", "LDHB", "TPI1", "GPI", "PGM1",
  "MYC", "JUN", "FOS", "EGR1", "CREB1", "SP1", "YY1", "E2F1",
  "CDK1", "CDK2", "CDK4", "CDK6", "CCNA2", "CCNB1", "CCND1", "CCNE1",
  "EGFR", "ERBB2", "KRAS", "BRAF", "PIK3CA", "PTEN", "AKT2", "AKT3",
  "VEGFA", "VEGFB", "VEGFC", "FLT1", "KDR", "FLT4",
  "CD4", "CD8A", "CD8B", "CD19", "CD3D", "CD3E", "CD3G",
  "FOXP3", "CTLA4", "PDCD1", "CD274", "LAG3", "HAVCR2", "TIGIT",
  "IL2", "IL4", "IL6", "IL10", "IL12A", "IL12B", "IL17A", "IFNG", "TGFB1"
)

all_genes <- unique(c(all_death_genes, background_genes))

n_samples <- 200
n_genes <- length(all_genes)

# 创建基础表达矩阵
base_expr <- matrix(
  rnorm(n_genes * n_samples, mean = 6, sd = 2),
  nrow = n_genes,
  ncol = n_samples
)

rownames(base_expr) <- all_genes
colnames(base_expr) <- paste0("Sample_", sprintf("%03d", 1:n_samples))

# 定义样本组
sample_groups <- c(rep("Normal", 100), rep("Tumor", 100))
tumor_idx <- 101:200

# 铁死亡相关基因差异表达
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

# 确保表达值为正
base_expr[base_expr < 0] <- 0.01

example_expr <- base_expr

# 生成临床数据
ferro_score <- colMeans(base_expr[ferro_drivers[ferro_drivers %in% rownames(base_expr)], ])
pyro_score <- colMeans(base_expr[pyro_genes[pyro_genes %in% rownames(base_expr)], ])
death_score <- scale(ferro_score) + scale(pyro_score)

base_survival <- 60
survival_time <- base_survival - death_score * 10 + rnorm(n_samples, mean = 0, sd = 15)
survival_time[survival_time < 1] <- runif(sum(survival_time < 1), 1, 6)
survival_time <- round(survival_time, 1)

death_prob <- pnorm(scale(death_score))
survival_status <- rbinom(n_samples, 1, prob = death_prob * 0.7 + 0.1)

example_clinical <- data.frame(
  sample_id = colnames(base_expr),
  group = sample_groups,
  OS_time = survival_time,
  OS_status = survival_status,
  stage = c(rep(NA, 100), sample(c("I", "II", "III", "IV"), 100, 
                                  replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))),
  grade = c(rep(NA, 100), sample(c("G1", "G2", "G3"), 100, 
                                  replace = TRUE, prob = c(0.25, 0.45, 0.3))),
  age = round(c(rnorm(100, mean = 45, sd = 12), rnorm(100, mean = 58, sd = 10))),
  gender = sample(c("Male", "Female"), n_samples, replace = TRUE, prob = c(0.55, 0.45)),
  treatment_response = c(rep(NA, 100), sample(c("CR", "PR", "SD", "PD"), 100,
                                               replace = TRUE, prob = c(0.15, 0.35, 0.30, 0.20))),
  stringsAsFactors = FALSE
)

example_clinical$age[example_clinical$age < 20] <- 20
example_clinical$age[example_clinical$age > 85] <- 85

# 保存示例数据
save(example_expr, file = "data/example_expr.rda", compress = "xz")
save(example_clinical, file = "data/example_clinical.rda", compress = "xz")

cat("\n========================================\n")
cat("所有数据文件已成功创建！\n")
cat("========================================\n\n")

cat("数据文件列表:\n")
cat("  - data/death_genesets.rda      (完整基因集)\n")
cat("  - data/death_genesets_simple.rda (简化基因集)\n")
cat("  - data/death_pathway_info.rda  (通路信息)\n")
cat("  - data/example_expr.rda        (示例表达矩阵)\n")
cat("  - data/example_clinical.rda    (示例临床数据)\n\n")

cat("示例数据信息:\n")
cat("  - example_expr: ", nrow(example_expr), " genes x ", ncol(example_expr), " samples\n", sep = "")
cat("  - example_clinical: ", nrow(example_clinical), " samples\n", sep = "")
cat("\n样本分布:\n")
print(table(example_clinical$group))

cat("\n请运行 devtools::document() 然后 devtools::install() 完成安装\n")
