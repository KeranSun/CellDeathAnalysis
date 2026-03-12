# =============================================================================
# 细胞死亡基因集数据准备脚本
# CellDeathAnalysis R Package
# 
# 运行此脚本生成 data/death_genesets.rda
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Ferroptosis (铁死亡)
# 来源: FerrDb V3 (http://www.zhounan.org/ferrdb/)
# 参考文献: Zhou et al., Nucleic Acids Res. 2023
# -----------------------------------------------------------------------------

ferroptosis <- list(
  # Drivers (促进铁死亡) - 264个基因，这里列出核心基因
  driver = c(
    "ACSL4", "LPCAT3", "ALOX5", "ALOX12", "ALOX15", "ALOX15B",
    "NCOA4", "TFRC", "SLC39A14", "SLC39A8", "HMOX1", "HMOX2",
    "STEAP3", "FTH1", "FTL", "IREB2", "ACO1",
    "NOX1", "NOX2", "NOX4", "CYBB",
    "PTGS2", "POR", "FDFT1",
    "ATG5", "ATG7", "BECN1",
    "CHAC1", "CBS", "CTH",
    "SAT1", "SPERMIDINE", "GLS2",
    "TP53", "CDKN1A", "ALOX12B",
    "CARS1", "SLC1A5", "GOT1",
    "HSPB1", "CISD1", "CISD2",
    "VDAC2", "VDAC3"
  ),
  
  # Suppressors (抑制铁死亡) - 238个基因，这里列出核心基因
  suppressor = c(
    "GPX4", "SLC7A11", "SLC3A2",
    "FSP1", "AIFM2", "DHODH", "GCH1",
    "NFE2L2", "KEAP1",
    "GCLC", "GCLM", "GSS", "GSR",
    "NQO1", "TXNRD1", "TXN",
    "FPN1", "SLC40A1", "HAMP",
    "PROM2", "MT1G", "MT2A",
    "HSPA5", "HSPB1",
    "NFS1", "ISCU",
    "SQSTM1", "FANCD2",
    "CDKN2A", "ARNTL",
    "SCD", "FADS2",
    "ACSL3", "HILPDA",
    "CoQ10A", "CoQ10B"
  ),
  
  # Markers (铁死亡标志物)
  marker = c(
    "PTGS2", "CHAC1", "ACSL4", "TFRC",
    "SLC7A11", "GPX4", "HMOX1",
    "FTH1", "FTL", "NCOA4"
  )
)

# -----------------------------------------------------------------------------
# 2. Cuproptosis (铜死亡)
# 来源: Tsvetkov et al., Science 2022 (DOI: 10.1126/science.abf0529)
# -----------------------------------------------------------------------------

cuproptosis <- list(
  # 正调控基因 - 脂酸途径和脂酰化蛋白
  positive = c(
    "FDX1",    # 铁氧还蛋白1 - 核心调控因子
    "LIAS",    # 硫辛酸合成酶
    "LIPT1",   # 硫辛酸转移酶1
    "LIPT2",   # 硫辛酸转移酶2 (新增)
    "DLD",     # 二氢硫辛酰胺脱氢酶
    "DLAT",    # 丙酮酸脱氢酶E2亚基
    "PDHA1",   # 丙酮酸脱氢酶E1α亚基
    "PDHB",    # 丙酮酸脱氢酶E1β亚基
    "DLST",    # α-酮戊二酸脱氢酶E2亚基
    "GCSH",    # 甘氨酸裂解系统H蛋白
    "DBT"      # 支链α-酮酸脱氢酶E2亚基
  ),
  
  # 负调控基因
  negative = c(
    "MTF1",    # 金属调节转录因子1
    "GLS",     # 谷氨酰胺酶
    "CDKN2A"   # p16/p14ARF
  ),
  
  # 铜转运相关基因
  transport = c(
    "SLC31A1", # CTR1 - 高亲和力铜摄入
    "SLC31A2", # CTR2
    "ATP7A",   # 铜外排 (门克斯病相关)
    "ATP7B",   # 铜外排 (威尔逊病相关)
    "ATOX1",   # 铜伴侣蛋白
    "CCS",     # SOD1铜伴侣
    "COX17",   # 细胞色素c氧化酶铜伴侣
    "SCO1",    # 线粒体铜伴侣
    "SCO2"     # 线粒体铜伴侣
  )
)

# -----------------------------------------------------------------------------
# 3. Disulfidptosis (二硫死亡)
# 来源: Liu et al., Nature Cell Biology 2023 
# DOI: 10.1038/s41556-023-01091-2
# -----------------------------------------------------------------------------

disulfidptosis <- list(
  # 核心基因
  core = c(
    "SLC7A11",  # xCT - 胱氨酸转运体，高表达是前提
    "SLC3A2",   # 4F2hc - SLC7A11的伴侣
    "NDUFS1",   # NADH脱氢酶铁硫蛋白1
    "NUBPL",    # 铁硫簇组装因子
    "LRPPRC",   # 线粒体RNA结合蛋白
    "NDUFA11",  # 线粒体复合物I亚基
    "OXSM",     # 3-酮脂酰ACP合成酶
    "GYS1"      # 糖原合成酶1
  ),
  
  # 细胞骨架相关 (二硫键形成导致骨架塌陷)
  cytoskeleton = c(
    "ACTN4",    # α-辅肌动蛋白4
    "MYH9",     # 肌球蛋白重链9
    "MYH10",    # 肌球蛋白重链10
    "FLNA",     # 丝蛋白A
    "FLNB",     # 丝蛋白B
    "TLN1",     # 踝蛋白1
    "IQGAP1"    # IQ基序GTP酶激活蛋白1
  ),
  
  # 代谢相关
  metabolism = c(
    "RPN1",     # 核糖体蛋白
    "SLC2A1",   # GLUT1 - 葡萄糖转运体
    "SLC2A3",   # GLUT3
    "HK1",      # 己糖激酶1
    "HK2"       # 己糖激酶2
  )
)

# -----------------------------------------------------------------------------
# 4. Pyroptosis (焦亡)
# 来源: MSigDB GOBP_PYROPTOSIS + KEGG + 文献整理
# -----------------------------------------------------------------------------

pyroptosis <- list(
  # 炎症小体成分
  inflammasome = c(
    "NLRP1", "NLRP3", "NLRP6", "NLRP7", "NLRP12",
    "NLRC4", "NAIP",
    "AIM2",     # dsDNA感受器
    "PYCARD",   # ASC
    "MEFV"      # Pyrin
  ),
  
  # Caspases
  caspases = c(
    "CASP1",    # 经典途径
    "CASP4",    # 非经典途径 (人)
    "CASP5",    # 非经典途径 (人)
    "CASP3",    # 可切割GSDME
    "CASP8"     # 交叉调控
  ),
  
  # Gasdermin家族 - 执行者
  gasdermins = c(
    "GSDMA",
    "GSDMB",
    "GSDMC",
    "GSDMD",    # 主要执行者
    "GSDME",    # DFNA5
    "PJVK"      # DFNB59/GSDMF
  ),
  
  # 炎症因子
  cytokines = c(
    "IL1A", "IL1B", "IL18", "IL33",
    "HMGB1",    # DAMPs
    "S100A8", "S100A9"
  ),
  
  # 上游信号和调控因子
  regulators = c(
    "TLR2", "TLR4", "TLR7", "TLR9",
    "NFKB1", "RELA",
    "NEK7",     # NLRP3激活必需
    "P2RX7",    # ATP受体
    "PANX1",    # Pannexin-1
    "GBP1", "GBP2", "GBP5",  # 干扰素诱导
    "IRGM",
    "TXNIP",    # NLRP3激活
    "BRCC3",    # NLRP3去泛素化
    "SIRT2"     # NLRP3乙酰化调控
  )
)

# -----------------------------------------------------------------------------
# 5. Necroptosis (坏死性凋亡)
# 来源: KEGG hsa04217 + MSigDB
# -----------------------------------------------------------------------------

necroptosis <- list(
  # 核心信号复合物
  core = c(
    "RIPK1",    # RIP1
    "RIPK3",    # RIP3
    "MLKL",     # 执行者
    "FADD",
    "TRADD",
    "CASP8",    # 抑制坏死性凋亡
    "CFLAR"     # c-FLIP
  ),
  
  # 死亡受体
  receptors = c(
    "TNFRSF1A", # TNFR1
    "TNFRSF1B", # TNFR2
    "FAS",      # CD95
    "TNFRSF10A", # DR4/TRAIL-R1
    "TNFRSF10B", # DR5/TRAIL-R2
    "TNFRSF6B"   # DcR3
  ),
  
  # 配体
  ligands = c(
    "TNF",
    "FASLG",    # FasL
    "TNFSF10"   # TRAIL
  ),
  
  # 模式识别受体途径
  prr = c(
    "TLR3", "TLR4",
    "TICAM1",   # TRIF
    "ZBP1",     # DAI/DLM-1
    "IFNAR1", "IFNAR2"
  ),
  
  # 泛素化调控
  ubiquitination = c(
    "CYLD",     # 去泛素化酶
    "BIRC2",    # cIAP1
    "BIRC3",    # cIAP2
    "XIAP",
    "LUBAC",
    "SHARPIN",
    "RNF31",    # HOIP
    "RBCK1",    # HOIL-1
    "OTULIN",
    "TNFAIP3",  # A20
    "SPATA2"
  ),
  
  # 激酶信号
  kinases = c(
    "MAP3K7",   # TAK1
    "TAB1", "TAB2", "TAB3",
    "IKBKB",    # IKKβ
    "CHUK",     # IKKα
    "IKBKG",    # NEMO
    "MAPK8",    # JNK1
    "MAPK9",    # JNK2
    "MAPK14"    # p38α
  )
)

# -----------------------------------------------------------------------------
# 6. Apoptosis (凋亡)
# 来源: KEGG hsa04210 + HALLMARK_APOPTOSIS
# -----------------------------------------------------------------------------

apoptosis <- list(
  # 外源性途径 (死亡受体途径)
  extrinsic = c(
    # 死亡受体
    "FAS", "TNFRSF1A", "TNFRSF10A", "TNFRSF10B",
    # 配体
    "FASLG", "TNF", "TNFSF10",
    # 接头蛋白
    "FADD", "TRADD",
    # 启动Caspases
    "CASP8", "CASP10"
  ),
  
  # 内源性途径 (线粒体途径)
  intrinsic = c(
    # 促凋亡 Bcl-2家族
    "BAX", "BAK1", "BOK",
    # BH3-only蛋白
    "BID", "BIM", "BAD", "BIK", "BMF", "HRK",
    "PMAIP1",   # NOXA
    "BBC3",     # PUMA
    # 抗凋亡 Bcl-2家族
    "BCL2", "BCL2L1", "BCL2L2", "MCL1", "BCL2A1",
    # 线粒体释放因子
    "CYCS",     # 细胞色素c
    "DIABLO",   # SMAC
    "HTRA2",    # Omi
    "ENDOG",    # 核酸内切酶G
    "AIFM1",    # AIF
    # Apoptosome
    "APAF1",
    "CASP9"
  ),
  
  # 执行Caspases
  executioner = c(
    "CASP3", "CASP6", "CASP7"
  ),
  
  # IAP家族
  iap = c(
    "XIAP", "BIRC2", "BIRC3", "BIRC5", "BIRC6", "BIRC7", "NAIP"
  ),
  
  # p53途径
  p53_pathway = c(
    "TP53", "MDM2", "CDKN1A",
    "GADD45A", "GADD45B", "GADD45G",
    "SESN1", "SESN2",
    "PIDD1", "RAIDD", "CASP2"
  ),
  
  # ER应激相关
  er_stress = c(
    "DDIT3",    # CHOP
    "ERN1",     # IRE1
    "EIF2AK3",  # PERK
    "ATF4", "ATF6",
    "CASP12"
  )
)

# -----------------------------------------------------------------------------
# 7. Autophagy (自噬)
# 来源: HADb (Human Autophagy Database) + KEGG hsa04140
# -----------------------------------------------------------------------------

autophagy <- list(
  # ULK1复合物 (起始)
  ulk_complex = c(
    "ULK1", "ULK2",
    "ATG13",
    "RB1CC1",   # FIP200
    "ATG101"
  ),
  
  # PI3K复合物
  pi3k_complex = c(
    "BECN1",    # Beclin-1
    "PIK3C3",   # VPS34
    "PIK3R4",   # VPS15
    "ATG14",
    "AMBRA1",
    "UVRAG",
    "SH3GLB1",  # Bif-1
    "RUBCN"     # Rubicon (负调控)
  ),
  
  # ATG9循环系统
  atg9_system = c(
    "ATG9A", "ATG9B",
    "WIPI1", "WIPI2",
    "ATG2A", "ATG2B",
    "WDR45",    # WIPI4
    "WDR45B"    # WIPI3
  ),
  
  # ATG12结合系统
  atg12_system = c(
    "ATG5", "ATG7", "ATG10", "ATG12", "ATG16L1", "ATG16L2"
  ),
  
  # LC3结合系统
  lc3_system = c(
    "ATG3", "ATG4A", "ATG4B", "ATG4C", "ATG4D", "ATG7"
  ),
  
  # LC3/GABARAP家族
  lc3_family = c(
    "MAP1LC3A", "MAP1LC3B", "MAP1LC3B2", "MAP1LC3C",
    "GABARAP", "GABARAPL1", "GABARAPL2"
  ),
  
  # 选择性自噬受体
  receptors = c(
    "SQSTM1",   # p62
    "NBR1",
    "CALCOCO2", # NDP52
    "OPTN",
    "TAX1BP1",
    "TOLLIP",
    # 线粒体自噬
    "BNIP3", "BNIP3L", "FUNDC1", "PHB2",
    "PINK1", "PRKN",  # Parkin
    # 其他
    "NCOA4",    # 铁蛋白自噬
    "FAM134B",  # ER自噬
    "SEC62"
  ),
  
  # mTOR信号 (负调控)
  mtor_pathway = c(
    "MTOR", "RPTOR", "MLST8", "AKT1S1",
    "TSC1", "TSC2", "RHEB",
    "RPS6KB1",  # S6K1
    "EIF4EBP1"
  ),
  
  # AMPK信号 (正调控)
  ampk_pathway = c(
    "PRKAA1", "PRKAA2",  # AMPKα
    "PRKAB1", "PRKAB2",  # AMPKβ
    "PRKAG1", "PRKAG2", "PRKAG3",  # AMPKγ
    "STK11",    # LKB1
    "STRADA", "STRADB", "CAB39"
  ),
  
  # 其他调控因子
  regulators = c(
    "TP53", "DRAM1", "DAPK1",
    "TFEB", "TFE3",  # 溶酶体生成
    "FOXO1", "FOXO3",
    "SIRT1",
    "HDAC6",
    "VCP",      # p97
    "SQSTM1", "KEAP1", "NFE2L2"
  ),
  
  # 溶酶体融合
  lysosome_fusion = c(
    "RAB7A", "LAMP1", "LAMP2",
    "STX17", "SNAP29", "VAMP8",
    "PLEKHM1", "EPG5",
    "HOPS",     # 复合物
    "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39", "VPS41"
  )
)

# -----------------------------------------------------------------------------
# 8. PANoptosis (泛凋亡)
# 整合pyroptosis + apoptosis + necroptosis
# -----------------------------------------------------------------------------

panoptosis <- list(
  # 核心复合物成分 (PANoptosome)
  panoptosome = c(
    "ZBP1",     # Z-DNA结合蛋白 - 关键感受器
    "RIPK1",
    "RIPK3",
    "CASP8",    # 分子开关
    "FADD",
    "NLRP3",
    "PYCARD",   # ASC
    "CASP1",
    "CASP6"
  ),
  
  # 效应分子
  effectors = c(
    "MLKL",     # 坏死性凋亡执行
    "GSDMD",    # 焦亡执行
    "GSDME",    # 焦亡执行
    "CASP3",    # 凋亡执行
    "CASP7"
  ),
  
  # 交叉调控因子
  crosstalk = c(
    "BID",      # 连接外源和内源凋亡
    "BAX", "BAK1",
    "CYCS",
    "AIFM1",
    "CFLAR",    # c-FLIP
    "TNFAIP3"   # A20
  ),
  
  # 上游信号
  upstream = c(
    "TNF", "TNFRSF1A",
    "IFNG", "IFNAR1",
    "TLR3", "TLR4",
    "IRF1"
  )
)

# -----------------------------------------------------------------------------
# 9. NETosis (中性粒细胞胞外陷阱)
# -----------------------------------------------------------------------------

netosis <- list(
  core = c(
    "MPO",      # 髓过氧化物酶
    "ELANE",    # 中性粒细胞弹性蛋白酶
    "CTSG",     # 组织蛋白酶G
    "PRTN3",    # 蛋白酶3
    "PADI4",    # PAD4 - 组蛋白瓜氨酸化
    "PADI2"
  ),
  
  ros_generation = c(
    "CYBB",     # NOX2/gp91phox
    "NCF1", "NCF2", "NCF4",  # p47phox, p67phox, p40phox
    "CYBA",     # p22phox
    "RAC1", "RAC2"
  ),
  
  chromatin = c(
    "H2AFX", "H2AFY",
    "H3F3A", "H3F3B",
    "HIST1H2BC", "HIST1H4A"
  ),
  
  signaling = c(
    "PKC",      # 蛋白激酶C
    "PRKCA", "PRKCB", "PRKCD",
    "RAF1", "MAPK1", "MAPK3",
    "SYK",
    "PI3K", "AKT1"
  ),
  
  released_factors = c(
    "HMGB1",
    "CAMP",     # LL-37/cathelicidin
    "S100A8", "S100A9", "S100A12",
    "DEFA1", "DEFA3", "DEFA4",
    "LCN2",     # NGAL
    "MMP9"
  )
)

# -----------------------------------------------------------------------------
# 10. Parthanatos (PARP依赖性细胞死亡)
# -----------------------------------------------------------------------------

parthanatos <- list(
  core = c(
    "PARP1",    # 核心因子
    "AIFM1",    # AIF - 执行者
    "MIF",      # 巨噬细胞迁移抑制因子 - AIF核酸酶活性激活
    "PARG",     # PAR糖水解酶
    "ADPRHL2",  # ARH3
    "RNLS"      # Renalase
  ),
  
  dna_damage = c(
    "ATM", "ATR",
    "H2AFX",    # γH2AX
    "XRCC1",
    "PARP2"
  ),
  
  nad_metabolism = c(
    "NAMPT",
    "NMNAT1", "NMNAT2", "NMNAT3",
    "SIRT1", "SIRT2"
  )
)

# -----------------------------------------------------------------------------
# 11. Entosis (细胞内吞)
# -----------------------------------------------------------------------------

entosis <- list(
  core = c(
    "RHOA",
    "ROCK1", "ROCK2",
    "CDC42",
    "RAC1",
    "ACTB", "ACTG1",
    "MYH9", "MYH10",
    "MYL9", "MYL12A", "MYL12B"
  ),
  
  adhesion = c(
    "CDH1",     # E-cadherin
    "CTNNB1",   # β-catenin
    "CTNNA1",   # α-catenin
    "JUP"       # γ-catenin
  ),
  
  degradation = c(
    "LAMP1", "LAMP2",
    "CTSL", "CTSD",
    "RAB7A"
  )
)

# -----------------------------------------------------------------------------
# 12. Oxeiptosis (氧化诱导细胞死亡)
# 来源: Holze et al., Nature 2018
# -----------------------------------------------------------------------------

oxeiptosis <- list(
  core = c(
    "KEAP1",    # 感受器
    "PGAM5",    # 磷酸酶
    "AIFM1"     # AIF - 执行者
  ),
  
  related = c(
    "NFE2L2",   # NRF2
    "NQO1",
    "HMOX1",
    "SOD1", "SOD2",
    "CAT",
    "GPX1", "GPX4"
  )
)

# -----------------------------------------------------------------------------
# 13. Alkaliptosis (碱中毒细胞死亡)
# 来源: Song et al., PNAS 2018
# -----------------------------------------------------------------------------

alkaliptosis <- list(
  core = c(
    "CA9",      # 碳酸酐酶9
    "SLC4A7",   # NBC3 - 碳酸氢盐转运
    "SLC9A1",   # NHE1 - Na+/H+交换
    "ATP6V1A",  # V-ATPase
    "ATP6V0D1"
  ),
  
  related = c(
    "JTC801",   # 诱导剂靶点
    "NFKB1", "RELA"
  )
)

# -----------------------------------------------------------------------------
# 14. Lysosome-dependent cell death (溶酶体依赖性细胞死亡)
# -----------------------------------------------------------------------------

ldcd <- list(
  core = c(
    "LAMP1", "LAMP2",
    "CTSL", "CTSB", "CTSD", "CTSH",
    "LGMN",     # Legumain
    "ACP2"      # 酸性磷酸酶
  ),
  
  regulators = c(
    "TFEB", "TFE3",
    "HSP70",    # HSPA1A
    "HSPA1A", "HSPA1B",
    "BID",
    "BAX"
  )
)

# =============================================================================
# 整合所有基因集
# =============================================================================

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

# 添加元数据
attr(death_genesets, "version") <- "1.0.0"
attr(death_genesets, "date") <- Sys.Date()
attr(death_genesets, "species") <- "Homo sapiens"

# =============================================================================
# 创建简化版基因集 (每种死亡类型只包含核心基因)
# =============================================================================

death_genesets_simple <- list(
  ferroptosis = unique(c(ferroptosis$driver[1:20], ferroptosis$suppressor[1:15], ferroptosis$marker)),
  cuproptosis = unique(c(cuproptosis$positive, cuproptosis$negative, cuproptosis$transport)),
  disulfidptosis = unique(c(disulfidptosis$core, disulfidptosis$cytoskeleton)),
  pyroptosis = unique(c(pyroptosis$inflammasome, pyroptosis$caspases, pyroptosis$gasdermins)),
  necroptosis = unique(c(necroptosis$core, necroptosis$receptors[1:4])),
  apoptosis = unique(c(apoptosis$extrinsic, apoptosis$intrinsic[1:15], apoptosis$executioner)),
  autophagy = unique(c(autophagy$ulk_complex, autophagy$pi3k_complex, autophagy$lc3_family, autophagy$receptors[1:8])),
  panoptosis = unique(c(panoptosis$panoptosome, panoptosis$effectors)),
  netosis = unique(c(netosis$core, netosis$ros_generation[1:4])),
  parthanatos = parthanatos$core,
  entosis = entosis$core,
  oxeiptosis = oxeiptosis$core
)

attr(death_genesets_simple, "version") <- "1.0.0"
attr(death_genesets_simple, "date") <- Sys.Date()
attr(death_genesets_simple, "species") <- "Homo sapiens"
attr(death_genesets_simple, "description") <- "Simplified gene sets containing core genes only"

# =============================================================================
# 保存数据
# =============================================================================

usethis::use_data(death_genesets, overwrite = TRUE)
usethis::use_data(death_genesets_simple, overwrite = TRUE)

cat("基因集数据已保存!\n")
cat("包含", length(death_genesets), "种细胞死亡类型\n")
cat("完整版基因集: death_genesets\n")
cat("简化版基因集: death_genesets_simple\n")
