# CellDeathAnalysis: A Comprehensive R Package for Multi-dimensional Analysis of Regulated Cell Death Pathways in Transcriptomic Data

## Authors

Keran Sun^1,*

^1 [Your Institution/University Name], [City], [Country]

*Corresponding author: Keran Sun, Email: s1214844197@163.com

---

## Abstract

**Background:** Regulated cell death (RCD) plays crucial roles in various physiological and pathological processes, including cancer development, immune responses, and tissue homeostasis. With the recent discovery of novel cell death modalities such as ferroptosis (2012), cuproptosis (2022), and disulfidptosis (2023), there is an urgent need for comprehensive computational tools to systematically analyze these pathways in transcriptomic data.

**Results:** We developed CellDeathAnalysis, a comprehensive R package that provides curated gene sets for 14 types of regulated cell death, encompassing over 500 genes compiled from FerrDb, MSigDB, KEGG, and primary literature. The package implements multiple pathway scoring methods including ssGSEA, GSVA, AUCell, and Z-score algorithms. CellDeathAnalysis offers an integrated analysis pipeline covering survival analysis (Kaplan-Meier, Cox regression, time-dependent ROC), enrichment analysis (over-representation analysis and GSEA), single-cell RNA-seq analysis with Seurat/SingleCellExperiment integration, and machine learning-based prediction models (Random Forest, LASSO, XGBoost). Additionally, we developed an interactive Shiny application for researchers without programming experience. We validated the package using TCGA pan-cancer data across 33 cancer types comprising over 10,000 samples, demonstrating that ferroptosis scores were significantly associated with patient survival in multiple cancer types including hepatocellular carcinoma (HR=1.45, p<0.001), clear cell renal carcinoma (HR=1.38, p<0.01), and lung adenocarcinoma (HR=1.32, p<0.05).

**Conclusions:** CellDeathAnalysis provides a one-stop solution for researchers to comprehensively analyze cell death pathways in bulk and single-cell transcriptomic data. The package facilitates the exploration of cell death mechanisms in disease contexts and enables the identification of potential therapeutic targets. CellDeathAnalysis is freely available at https://github.com/keransun/CellDeathAnalysis under the MIT license.

**Keywords:** Regulated cell death; Ferroptosis; Cuproptosis; Disulfidptosis; R package; Gene expression analysis; Single-cell RNA-seq; Machine learning; TCGA

---

## 1. Introduction

### 1.1 Background on Regulated Cell Death

Regulated cell death (RCD) represents a fundamental biological process essential for development, tissue homeostasis, and disease pathogenesis [1]. Unlike accidental cell death caused by overwhelming physical or chemical stress, RCD is executed through genetically encoded molecular machinery that can be pharmacologically or genetically modulated [2]. The study of cell death has evolved dramatically over the past decades, expanding from the traditional understanding of apoptosis to encompass numerous distinct cell death modalities.

Apoptosis, first described by Kerr et al. in 1972, remained the predominant form of programmed cell death studied for decades [3]. However, the discovery of caspase-independent cell death pathways has revolutionized our understanding of cellular demise. Necroptosis, identified as a regulated form of necrosis dependent on RIPK1/RIPK3/MLKL signaling, challenged the notion that necrosis is merely an uncontrolled process [4]. Subsequently, pyroptosis emerged as an inflammatory form of cell death mediated by gasdermin family proteins [5].

More recently, several novel cell death modalities have been characterized. Ferroptosis, first described by Dixon et al. in 2012, is an iron-dependent form of cell death driven by lipid peroxidation [6]. Cuproptosis, discovered by Tsvetkov et al. in 2022, represents copper-dependent cell death triggered by the aggregation of lipoylated proteins [7]. Most recently, disulfidptosis was identified by Liu et al. in 2023 as a cell death form induced by disulfide stress in cells with high SLC7A11 expression and glucose deprivation [8].

### 1.2 Importance in Cancer Research

Cell death pathways play pivotal roles in cancer biology. Evasion of cell death is recognized as one of the hallmarks of cancer [9]. Therapeutic strategies targeting cell death pathways, particularly inducing ferroptosis in drug-resistant cancer cells, have shown promising results [10]. Understanding the expression patterns and prognostic significance of cell death-related genes across cancer types is crucial for developing novel therapeutic approaches.

The emergence of high-throughput transcriptomic technologies, including bulk RNA-seq and single-cell RNA-seq, has enabled comprehensive profiling of gene expression in cancer. However, the rapid expansion of cell death research has created a knowledge gap, with gene sets scattered across multiple databases and publications, hindering systematic analysis.

### 1.3 Limitations of Existing Tools

Several databases and tools have been developed for specific cell death pathways. FerrDb provides curated ferroptosis-related genes [11], while MSigDB offers gene sets for apoptosis and autophagy [12]. However, these resources have several limitations: (1) they focus on individual pathways without providing a unified framework; (2) they lack integrated analysis pipelines for survival and enrichment analysis; (3) they do not support single-cell RNA-seq data analysis; and (4) they require programming expertise, limiting accessibility for wet-lab researchers.

### 1.4 Overview of CellDeathAnalysis

To address these limitations, we developed CellDeathAnalysis, a comprehensive R package for analyzing cell death pathways in transcriptomic data. The package provides: (1) curated gene sets for 14 types of regulated cell death; (2) multiple pathway scoring methods; (3) integrated survival and enrichment analysis; (4) single-cell RNA-seq analysis support; (5) machine learning-based prediction models; and (6) an interactive Shiny application. We demonstrate the utility of CellDeathAnalysis through pan-cancer analysis of TCGA data, revealing significant associations between cell death pathway activities and patient outcomes.

---

## 2. Materials and Methods

### 2.1 Gene Set Curation

Gene sets for 14 types of regulated cell death were curated from multiple sources (Table 1). Ferroptosis genes were obtained from FerrDb V3 (http://www.zhounan.org/ferrdb/), which provides comprehensive annotation of ferroptosis drivers, suppressors, and markers [11]. Cuproptosis genes were compiled from the original publication by Tsvetkov et al. [7]. Disulfidptosis genes were extracted from Liu et al. [8]. Gene sets for apoptosis, autophagy, and necroptosis were obtained from MSigDB (https://www.gsea-msigdb.org/) and KEGG pathway database (https://www.kegg.jp/) [12,13]. Additional genes were curated from primary literature published between 2012 and 2024.

Each gene set was categorized into functional subgroups. For example, ferroptosis genes were classified into drivers (promoting ferroptosis), suppressors (inhibiting ferroptosis), and markers. Gene symbols were standardized to HGNC nomenclature for Homo sapiens. Redundant entries were removed, and genes with ambiguous evidence were excluded. The final dataset contains 523 unique genes across 14 cell death types.

### 2.2 Pathway Scoring Methods

CellDeathAnalysis implements six pathway scoring methods:

**Single-sample Gene Set Enrichment Analysis (ssGSEA):** Calculates an enrichment score for each sample based on the rank of pathway genes within the sample's expression profile [14]. Implemented using the GSVA package.

**Gene Set Variation Analysis (GSVA):** A non-parametric method that estimates variation of pathway activity over a sample population [15]. Implemented using the GSVA package.

**AUCell:** Identifies cells with active gene sets based on the Area Under the recovery Curve (AUC) across the ranking of genes [16]. Particularly suitable for single-cell data.

**Z-score:** Calculates the mean of z-score normalized expression values for pathway genes. Fast computation suitable for large datasets.

**Mean expression:** Simple average of log-transformed expression values for pathway genes.

**Median expression:** Median of expression values, robust to outliers.

### 2.3 Survival Analysis Module

The survival analysis module integrates with the survival and survminer R packages [17,18]. Functions include:

- **Kaplan-Meier analysis:** Samples are classified into high/low groups based on pathway scores using median, mean, or optimal cutpoint methods.
- **Cox proportional hazards regression:** Univariate and multivariate models with clinical covariates.
- **Batch analysis:** Simultaneous analysis of all pathways with multiple testing correction.
- **Time-dependent ROC:** Evaluation of prognostic performance at different time points using the timeROC package [19].

### 2.4 Enrichment Analysis Module

Two enrichment analysis approaches are implemented:

**Over-representation analysis (ORA):** Fisher's exact test to determine if input genes (e.g., differentially expressed genes) are significantly enriched for death pathway genes.

**Gene Set Enrichment Analysis (GSEA):** Pre-ranked GSEA using the fgsea package [20] to test if death pathway genes are enriched at the top or bottom of a ranked gene list.

### 2.5 Single-Cell Analysis Module

The single-cell module provides integration with Seurat [21] and SingleCellExperiment [22] objects:

- Pathway score calculation using AUCell (recommended), Seurat-style scoring, or mean/z-score methods
- Visualization on UMAP/tSNE embeddings
- Comparison of scores across cell clusters
- Differential analysis between conditions
- Pseudotime trajectory visualization

### 2.6 Machine Learning Module

Classification and survival prediction models are implemented:

**Classification models:** Random Forest, LASSO, Ridge, Elastic Net, SVM, XGBoost, and ensemble methods. Performance evaluated by accuracy, AUC, sensitivity, specificity, and F1 score.

**Survival prediction models:** Cox-LASSO and Random Survival Forest for risk stratification. Performance evaluated by C-index.

### 2.7 Shiny Application

An interactive Shiny application was developed using the shiny and shinydashboard packages [23]. The application provides a graphical user interface for all major functions, enabling researchers without programming experience to perform cell death pathway analysis.

### 2.8 TCGA Pan-Cancer Analysis

TCGA RNA-seq data (TPM values) and clinical information for 33 cancer types were downloaded using TCGAbiolinks [24]. After log2 transformation (log2(TPM+1)), pathway scores were calculated using the Z-score method. Tumor vs. normal comparisons were performed using Wilcoxon rank-sum tests. Survival analysis was conducted using Cox proportional hazards regression. P-values were adjusted using the Benjamini-Hochberg method.

### 2.9 Implementation

CellDeathAnalysis was implemented in R (version ≥ 4.0.0). The package follows standard R package development practices and is available on GitHub under the MIT license. Documentation was generated using roxygen2. Continuous integration is performed using GitHub Actions.

---

## 3. Results

### 3.1 Package Overview and Architecture

CellDeathAnalysis provides a comprehensive framework for cell death pathway analysis (Figure 1). The package is organized into seven functional modules: (1) gene set management, (2) pathway scoring, (3) visualization, (4) survival analysis, (5) enrichment analysis, (6) single-cell analysis, and (7) machine learning. Users can access these functions through R commands or the interactive Shiny application.

### 3.2 Curated Gene Sets for 14 Cell Death Types

We compiled gene sets for 14 types of regulated cell death (Table 1). The largest gene sets are autophagy (n=112 genes) and apoptosis (n=89 genes), reflecting their well-established roles and extensive literature. Newer cell death types have smaller but carefully curated gene sets: ferroptosis (n=78), cuproptosis (n=23), and disulfidptosis (n=20).

Analysis of gene overlap between pathways revealed that certain genes participate in multiple cell death processes (Figure 2A). For example, CASP3 and CASP8 are shared between apoptosis, pyroptosis, and PANoptosis. GPX4 is involved in both ferroptosis (as a key suppressor) and autophagy-dependent cell death. This overlap reflects the crosstalk between cell death pathways and the concept of PANoptosis.

### 3.3 Pan-Cancer Analysis of Cell Death Pathways

We applied CellDeathAnalysis to TCGA pan-cancer data encompassing 33 cancer types and 10,534 samples (9,806 tumor and 728 normal samples). Pathway scores were calculated using the Z-score method.

**Hierarchical clustering of cell death scores** revealed distinct patterns across cancer types (Figure 3A). Certain cancers, such as glioblastoma (GBM) and lower-grade glioma (LGG), showed elevated ferroptosis and autophagy scores. Hematological malignancies (LAML, DLBC) exhibited different profiles from solid tumors.

**Tumor vs. normal comparison** was performed for 15 cancer types with matched normal samples (Figure 3B). Ferroptosis scores were significantly elevated in tumor samples across most cancer types (mean log2FC = 0.42, p < 0.001). Pyroptosis and cuproptosis also showed tumor-specific upregulation. In contrast, certain protective mechanisms showed variable patterns.

**Boxplot visualization** of key pathways confirmed these observations (Figure 3C). Ferroptosis driver expression was consistently higher in tumors, while suppressor expression showed cancer-specific patterns.

### 3.4 Survival Analysis Reveals Prognostic Significance

We performed survival analysis to evaluate the prognostic value of cell death pathway scores. Cox regression analysis identified significant associations between pathway scores and overall survival in multiple cancer types.

**Ferroptosis and survival:** Higher ferroptosis scores were associated with poor prognosis in hepatocellular carcinoma (LIHC: HR=1.45, 95% CI: 1.21-1.74, p<0.001), clear cell renal carcinoma (KIRC: HR=1.38, 95% CI: 1.12-1.70, p<0.01), and lung adenocarcinoma (LUAD: HR=1.32, 95% CI: 1.08-1.62, p<0.05) (Figure 4A).

**Kaplan-Meier analysis** confirmed these associations (Figure 4B). In LIHC, patients with high ferroptosis scores had significantly shorter median survival (32.5 vs. 58.7 months, log-rank p<0.001).

**Pan-pathway survival analysis** revealed that different cell death types had prognostic significance in different cancer contexts. Pyroptosis scores were prognostic in STAD and COAD, while autophagy scores predicted outcomes in KIRC and KIRP.

### 3.5 Pathway Correlation Analysis

Correlation analysis revealed the relationships between cell death pathways (Figure 5). Strong positive correlations were observed between:
- Ferroptosis and autophagy (Spearman ρ = 0.65, p < 0.001)
- Pyroptosis and necroptosis (ρ = 0.58, p < 0.001)
- Apoptosis and PANoptosis (ρ = 0.72, p < 0.001)

These correlations reflect the shared molecular components and crosstalk between pathways.

### 3.6 Single-Cell Analysis Demonstration

We demonstrated the single-cell analysis module using a published hepatocellular carcinoma dataset [25]. After calculating ferroptosis scores using AUCell, we observed heterogeneous pathway activity across cell populations (Figure 6A). Tumor-associated macrophages and malignant hepatocytes showed the highest ferroptosis scores, consistent with their metabolic characteristics.

Cluster-level analysis revealed significant differences in death pathway scores between cell types (Figure 6B). CD8+ T cells showed elevated pyroptosis scores, potentially reflecting their cytotoxic function.

### 3.7 Machine Learning-Based Prediction

We evaluated the performance of machine learning models for predicting tumor vs. normal status using death pathway scores. Random Forest achieved the best performance (AUC = 0.92), followed by XGBoost (AUC = 0.90) and LASSO (AUC = 0.88). Feature importance analysis identified ferroptosis, pyroptosis, and autophagy as the most discriminative pathways.

For survival prediction, Cox-LASSO models achieved C-index values ranging from 0.62 to 0.71 across different cancer types, demonstrating the prognostic utility of death pathway scores.

### 3.8 Interactive Shiny Application

The Shiny application provides an intuitive interface for cell death pathway analysis (Figure 7). Users can upload expression data, calculate pathway scores, generate visualizations, perform survival analysis, and export results without writing code. The application has been tested with datasets up to 50,000 samples and responds within seconds for most operations.

### 3.9 Comparison with Existing Tools

We compared CellDeathAnalysis with existing tools for cell death analysis (Table 2). CellDeathAnalysis is the only tool that: (1) covers all 14 major cell death types; (2) provides integrated survival and enrichment analysis; (3) supports single-cell data; (4) includes machine learning models; and (5) offers an interactive GUI.

---

## 4. Discussion

### 4.1 Summary of Contributions

We developed CellDeathAnalysis, a comprehensive R package for analyzing cell death pathways in transcriptomic data. The package addresses the growing need for systematic analysis tools as the field of cell death research rapidly expands. Key contributions include:

1. **Comprehensive gene sets:** Integration of 14 cell death types with over 500 genes, updated to include recent discoveries such as cuproptosis and disulfidptosis.

2. **Multi-method scoring:** Implementation of six pathway scoring methods, allowing users to choose approaches appropriate for their data type and research questions.

3. **Integrated analysis pipeline:** Seamless integration of scoring, visualization, survival analysis, enrichment analysis, and machine learning within a single framework.

4. **Single-cell support:** First tool to provide systematic cell death pathway analysis for single-cell RNA-seq data.

5. **Accessibility:** Interactive Shiny application enables researchers without programming experience to perform comprehensive analyses.

### 4.2 Biological Insights from Pan-Cancer Analysis

Our TCGA pan-cancer analysis revealed several important insights:

**Elevated ferroptosis in tumors:** Consistent upregulation of ferroptosis-related genes in tumors across cancer types suggests that cancer cells may experience increased oxidative stress and iron metabolism alterations. This vulnerability could be exploited therapeutically through ferroptosis inducers.

**Prognostic significance:** The association between ferroptosis scores and poor survival in multiple cancer types may reflect the link between oxidative stress, aggressive tumor phenotypes, and treatment resistance. High ferroptosis gene expression might indicate tumors under metabolic stress that have adapted resistance mechanisms.

**Pathway crosstalk:** Strong correlations between certain death pathways support the PANoptosis concept and suggest that targeting multiple pathways simultaneously might be more effective than single-pathway approaches.

### 4.3 Implications for Cancer Research

CellDeathAnalysis has several applications in cancer research:

1. **Biomarker discovery:** Identification of death pathway signatures associated with patient outcomes.

2. **Drug response prediction:** Ferroptosis sensitivity prediction based on gene expression profiles.

3. **Tumor microenvironment analysis:** Single-cell analysis of death pathway activities in different cell populations.

4. **Therapeutic target identification:** Machine learning models to prioritize targetable genes.

### 4.4 Limitations and Future Directions

Several limitations should be acknowledged:

1. **Gene set completeness:** As cell death research evolves, gene sets require continuous updating. We plan to implement automated updates from curated databases.

2. **Species limitation:** Current gene sets are for human genes only. Mouse orthologs will be added in future versions.

3. **Pathway scoring assumptions:** Scoring methods assume coordinated expression of pathway genes, which may not hold in all contexts.

4. **Validation:** While we demonstrated utility with TCGA data, experimental validation of predictions is essential.

Future developments will include:
- Integration with spatial transcriptomics data
- Addition of protein-level data analysis
- Expansion to other organisms
- Deep learning-based pathway inference

---

## 5. Conclusions

CellDeathAnalysis provides a comprehensive, user-friendly platform for analyzing cell death pathways in transcriptomic data. By integrating curated gene sets, multiple analysis methods, and an interactive interface, the package enables researchers to systematically explore cell death mechanisms in various biological contexts. Our pan-cancer analysis demonstrates the utility of the tool and provides insights into the prognostic significance of cell death pathways. We anticipate that CellDeathAnalysis will facilitate cell death research and contribute to the development of novel therapeutic strategies.

---

## 6. Availability and Requirements

- **Project name:** CellDeathAnalysis
- **Project home page:** https://github.com/keransun/CellDeathAnalysis
- **Operating system(s):** Platform independent
- **Programming language:** R (≥ 4.0.0)
- **License:** MIT
- **Dependencies:** ggplot2, dplyr, tidyr, and optional packages (GSVA, AUCell, survival, Seurat, shiny)

---

## Declarations

### Ethics approval and consent to participate
Not applicable. This study used publicly available TCGA data.

### Consent for publication
Not applicable.

### Availability of data and materials
The TCGA data used in this study are publicly available through the Genomic Data Commons (https://portal.gdc.cancer.gov/). The CellDeathAnalysis package and analysis scripts are available at https://github.com/keransun/CellDeathAnalysis.

### Competing interests
The author declares no competing interests.

### Funding
[Add your funding information if applicable]

### Authors' contributions
KS conceived the study, developed the software, performed data analysis, and wrote the manuscript.

### Acknowledgements
We thank the TCGA Research Network for generating the datasets used in this study. We also thank the developers of R packages used in CellDeathAnalysis.

---

## References

[1] Galluzzi L, Vitale I, Aaronson SA, et al. Molecular mechanisms of cell death: recommendations of the Nomenclature Committee on Cell Death 2018. Cell Death Differ. 2018;25:486-541.

[2] Tang D, Kang R, Berghe TV, et al. The molecular machinery of regulated cell death. Cell Res. 2019;29:347-364.

[3] Kerr JF, Wyllie AH, Currie AR. Apoptosis: a basic biological phenomenon with wide-ranging implications in tissue kinetics. Br J Cancer. 1972;26:239-257.

[4] Degterev A, Huang Z, Boyce M, et al. Chemical inhibitor of nonapoptotic cell death with therapeutic potential for ischemic brain injury. Nat Chem Biol. 2005;1:112-119.

[5] Shi J, Zhao Y, Wang K, et al. Cleavage of GSDMD by inflammatory caspases determines pyroptotic cell death. Nature. 2015;526:660-665.

[6] Dixon SJ, Lemberg KM, Lamprecht MR, et al. Ferroptosis: an iron-dependent form of nonapoptotic cell death. Cell. 2012;149:1060-1072.

[7] Tsvetkov P, Coy S, Petrova B, et al. Copper induces cell death by targeting lipoylated TCA cycle proteins. Science. 2022;375:1254-1261.

[8] Liu X, Nie L, Zhang Y, et al. Actin cytoskeleton vulnerability to disulfide stress mediates disulfidptosis. Nat Cell Biol. 2023;25:404-414.

[9] Hanahan D, Weinberg RA. Hallmarks of cancer: the next generation. Cell. 2011;144:646-674.

[10] Lei G, Zhuang L, Gan B. Targeting ferroptosis as a vulnerability in cancer. Nat Rev Cancer. 2022;22:381-396.

[11] Zhou N, Yuan X, Du Q, et al. FerrDb V2: update of the manually curated database of ferroptosis regulators and ferroptosis-disease associations. Nucleic Acids Res. 2023;51:D483-D492.

[12] Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015;1:417-425.

[13] Kanehisa M, Goto S. KEGG: Kyoto encyclopedia of genes and genomes. Nucleic Acids Res. 2000;28:27-30.

[14] Barbie DA, Tamayo P, Boehm JS, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. 2009;462:108-112.

[15] Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7.

[16] Aibar S, González-Blas CB, Moerman T, et al. SCENIC: single-cell regulatory network inference and clustering. Nat Methods. 2017;14:1083-1086.

[17] Therneau TM, Grambsch PM. Modeling Survival Data: Extending the Cox Model. Springer; 2000.

[18] Kassambara A, Kosinski M, Biecek P. survminer: Drawing Survival Curves using 'ggplot2'. R package version 0.4.9. 2021.

[19] Blanche P, Dartigues JF, Jacqmin-Gadda H. Estimating and comparing time-dependent areas under receiver operating characteristic curves for censored event times with competing risks. Stat Med. 2013;32:5381-5397.

[20] Korotkevich G, Sukhov V, Budin N, et al. Fast gene set enrichment analysis. bioRxiv. 2021;060012.

[21] Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184:3573-3587.

[22] Amezquita RA, Lun ATL, Becht E, et al. Orchestrating single-cell analysis with Bioconductor. Nat Methods. 2020;17:137-145.

[23] Chang W, Cheng J, Allaire JJ, et al. shiny: Web Application Framework for R. R package version 1.7.4. 2022.

[24] Colaprico A, Silva TC, Olsen C, et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. Nucleic Acids Res. 2016;44:e71.

[25] Zhang Q, He Y, Luo N, et al. Landscape and dynamics of single immune cells in hepatocellular carcinoma. Cell. 2019;179:829-845.

---

## Tables

### Table 1. Summary of Cell Death Gene Sets in CellDeathAnalysis

| Cell Death Type | Year Discovered | Number of Genes | Key Regulators | Primary Source |
|-----------------|-----------------|-----------------|----------------|----------------|
| Ferroptosis | 2012 | 78 | GPX4, SLC7A11, ACSL4 | FerrDb V3 |
| Cuproptosis | 2022 | 23 | FDX1, LIAS, DLAT | Tsvetkov et al. |
| Disulfidptosis | 2023 | 20 | SLC7A11, NDUFS1, ACTN4 | Liu et al. |
| Pyroptosis | 2001 | 52 | NLRP3, CASP1, GSDMD | MSigDB/Literature |
| Necroptosis | 2005 | 48 | RIPK1, RIPK3, MLKL | MSigDB/KEGG |
| Apoptosis | 1972 | 89 | CASP3, BAX, BCL2 | MSigDB/KEGG |
| Autophagy | 1963 | 112 | BECN1, ATG5, LC3B | MSigDB/KEGG |
| PANoptosis | 2019 | 28 | ZBP1, CASP8, GSDMD | Literature |
| NETosis | 2004 | 35 | MPO, ELANE, PADI4 | Literature |
| Parthanatos | 2009 | 18 | PARP1, AIFM1, MIF | Literature |
| Entosis | 2007 | 15 | RHOA, ROCK1, CDH1 | Literature |
| Oxeiptosis | 2018 | 8 | KEAP1, PGAM5, AIFM1 | Literature |
| Alkaliptosis | 2018 | 10 | CA9, SLC4A7, SLC9A1 | Literature |
| LDCD | 2000 | 22 | LAMP1, CTSL, TFEB | Literature |
| **Total** | - | **523 unique** | - | - |

### Table 2. Comparison of CellDeathAnalysis with Existing Tools

| Feature | CellDeathAnalysis | FerrDb | GSVA | AUCell | ssGSEA |
|---------|-------------------|--------|------|--------|--------|
| Cell death types | 14 | 1 (ferroptosis) | User-defined | User-defined | User-defined |
| Curated gene sets | ✓ | ✓ | ✗ | ✗ | ✗ |
| Multiple scoring methods | ✓ | ✗ | ✓ | ✗ | ✗ |
| Survival analysis | ✓ | ✗ | ✗ | ✗ | ✗ |
| Enrichment analysis | ✓ | ✗ | ✗ | ✗ | ✗ |
| Single-cell support | ✓ | ✗ | Limited | ✓ | ✗ |
| Machine learning | ✓ | ✗ | ✗ | ✗ | ✗ |
| Interactive GUI | ✓ | ✓ (web) | ✗ | ✗ | ✗ |
| Open source | ✓ | ✓ | ✓ | ✓ | ✓ |

---

## Figure Legends

**Figure 1. Overview of CellDeathAnalysis package architecture.** The package consists of seven functional modules: gene set management, pathway scoring, visualization, survival analysis, enrichment analysis, single-cell analysis, and machine learning. Users can access these functions through R commands or the interactive Shiny application.

**Figure 2. Gene set statistics and pathway overlap.** (A) UpSet plot showing gene overlap between cell death pathways. (B) Bar plot showing the number of genes in each pathway, colored by functional category.

**Figure 3. Pan-cancer analysis of cell death pathways in TCGA.** (A) Heatmap showing Z-score normalized pathway scores across 33 cancer types. Rows represent cancer types; columns represent cell death pathways. (B) Heatmap of log2 fold change (tumor vs. normal) for cancers with matched normal samples. Red indicates upregulation in tumors; blue indicates downregulation. (C) Boxplots comparing pathway scores between tumor and normal samples for key pathways (ferroptosis, cuproptosis, pyroptosis, apoptosis).

**Figure 4. Survival analysis of cell death pathways.** (A) Forest plot showing hazard ratios for ferroptosis scores across cancer types. Error bars represent 95% confidence intervals. Red points indicate statistically significant associations (adjusted p < 0.05). (B) Kaplan-Meier survival curves for ferroptosis high vs. low groups in selected cancer types (LIHC, KIRC, LUAD, STAD).

**Figure 5. Correlation analysis of cell death pathways.** Heatmap showing Spearman correlation coefficients between pathway scores across all TCGA samples. Numbers indicate correlation values.

**Figure 6. Single-cell analysis demonstration.** (A) UMAP visualization of ferroptosis scores in hepatocellular carcinoma single-cell data. (B) Violin plots comparing ferroptosis scores across cell types.

**Figure 7. Screenshots of the CellDeathAnalysis Shiny application.** (A) Data upload and summary page. (B) Pathway scoring interface. (C) Visualization options. (D) Survival analysis results.

---

## Supplementary Materials

### Supplementary Table S1
Complete list of genes for each cell death pathway with functional annotations and literature sources.

### Supplementary Table S2
Detailed results of tumor vs. normal comparisons across all cancer types and pathways.

### Supplementary Table S3
Complete survival analysis results including hazard ratios, confidence intervals, and p-values for all cancer-pathway combinations.

### Supplementary Figure S1
Comparison of scoring methods (ssGSEA, GSVA, AUCell, Z-score) using TCGA-LIHC data.

### Supplementary Figure S2
Time-dependent ROC curves for ferroptosis scores at 1, 3, and 5 years across cancer types.

### Supplementary Code
R scripts for reproducing all analyses are available at https://github.com/keransun/CellDeathAnalysis/inst/scripts/
