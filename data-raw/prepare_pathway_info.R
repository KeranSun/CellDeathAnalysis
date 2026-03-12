# =============================================================================
# 细胞死亡通路元数据
# =============================================================================

death_pathway_info <- data.frame(
  pathway = c(
    "ferroptosis",
    "cuproptosis", 
    "disulfidptosis",
    "pyroptosis",
    "necroptosis",
    "apoptosis",
    "autophagy",
    "panoptosis",
    "netosis",
    "parthanatos",
    "entosis",
    "oxeiptosis",
    "alkaliptosis",
    "ldcd"
  ),
  
  full_name = c(
    "Ferroptosis",
    "Cuproptosis",
    "Disulfidptosis", 
    "Pyroptosis",
    "Necroptosis",
    "Apoptosis",
    "Autophagy",
    "PANoptosis",
    "NETosis",
    "Parthanatos",
    "Entosis",
    "Oxeiptosis",
    "Alkaliptosis",
    "Lysosome-dependent cell death"
  ),
  
  chinese_name = c(
    "铁死亡",
    "铜死亡",
    "二硫死亡",
    "焦亡",
    "坏死性凋亡",
    "凋亡",
    "自噬",
    "泛凋亡",
    "中性粒细胞胞外陷阱相关死亡",
    "PARP依赖性细胞死亡",
    "细胞内吞死亡",
    "氧化诱导细胞死亡",
    "碱中毒细胞死亡",
    "溶酶体依赖性细胞死亡"
  ),
  
  year_discovered = c(
    2012,  # ferroptosis
    2022,  # cuproptosis
    2023,  # disulfidptosis
    2001,  # pyroptosis (named)
    2005,  # necroptosis
    1972,  # apoptosis (named)
    1963,  # autophagy
    2019,  # panoptosis
    2004,  # netosis
    2009,  # parthanatos (named)
    2007,  # entosis
    2018,  # oxeiptosis
    2018,  # alkaliptosis
    2000   # ldcd
  ),
  
  key_paper = c(
    "Dixon et al., Cell 2012",
    "Tsvetkov et al., Science 2022",
    "Liu et al., Nat Cell Biol 2023",
    "Cookson & Bhrennan, Trends Mol Med 2001",
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
    "Iron-dependent lipid peroxidation leading to cell death. Characterized by accumulation of lipid ROS and failure of GPX4 antioxidant system.",
    "Copper-dependent cell death through lipoylated protein aggregation and Fe-S cluster protein destabilization in mitochondria.",
    "Cell death caused by disulfide stress due to abnormal accumulation of cystine under glucose starvation in SLC7A11-high cells.",
    "Inflammatory cell death mediated by gasdermin pore formation following inflammasome activation and caspase-1 cleavage.",
    "Programmed necrosis mediated by RIPK1-RIPK3-MLKL signaling pathway, independent of caspases.",
    "Classical programmed cell death characterized by caspase activation, DNA fragmentation, and apoptotic body formation.",
    "Self-degradation process involving formation of autophagosomes and fusion with lysosomes for cellular component recycling.",
    "Coordinated cell death involving simultaneous activation of pyroptosis, apoptosis, and necroptosis pathways.",
    "Cell death in neutrophils accompanied by release of neutrophil extracellular traps (NETs) containing DNA, histones, and antimicrobial proteins.",
    "Cell death dependent on PARP1 hyperactivation and AIF translocation to nucleus.",
    "Cell-in-cell death where one cell engulfs and kills another live cell.",
    "ROS-induced cell death mediated by KEAP1-PGAM5-AIFM1 pathway.",
    "Cell death induced by intracellular alkalinization.",
    "Cell death caused by lysosomal membrane permeabilization and release of cathepsins."
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
    "Cell shrinkage, mitochondrial damage, lipid droplet accumulation",
    "Mitochondrial shrinkage, protein aggregation",
    "Cytoskeletal collapse, F-actin contraction",
    "Cell swelling, plasma membrane rupture, IL-1β/IL-18 release",
    "Cell swelling, plasma membrane rupture, organelle swelling",
    "Cell shrinkage, chromatin condensation, apoptotic bodies",
    "Autophagic vacuoles, no inflammation",
    "Mixed features of pyroptosis, apoptosis, and necroptosis",
    "NET release, chromatin decondensation",
    "Large-scale DNA fragmentation, chromatin condensation",
    "Cell-in-cell structures",
    "Caspase-independent, no inflammation",
    "Cytoplasmic vacuolization",
    "Lysosomal rupture, cytoplasmic acidification"
  ),
  
  stringsAsFactors = FALSE
)

usethis::use_data(death_pathway_info, overwrite = TRUE)
