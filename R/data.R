#' Cell Death Pathway Gene Sets
#'
#' A comprehensive collection of gene sets for various cell death pathways.
#' This dataset contains genes associated with 14 different types of regulated
#' cell death, curated from multiple databases and literature sources.
#'
#' @format A named list containing 14 cell death types, each with sub-lists
#' of genes organized by function:
#' \describe{
#'   \item{ferroptosis}{Iron-dependent cell death. Contains drivers, suppressors,
#'   and markers. Source: FerrDb V3}
#'   \item{cuproptosis}{Copper-dependent cell death. Contains positive regulators,
#'   negative regulators, and transport genes. Source: Tsvetkov et al., Science 2022}
#'   \item{disulfidptosis}{Disulfide stress-induced cell death. Contains core genes,
#'   cytoskeleton-related genes, and metabolism genes. Source: Liu et al., Nat Cell Biol 2023}
#'   \item{pyroptosis}{Inflammatory cell death. Contains inflammasome components,
#'   caspases, gasdermins, cytokines, and regulators}
#'   \item{necroptosis}{Programmed necrosis. Contains core signaling components,
#'   receptors, ligands, PRR pathway genes, ubiquitination regulators, and kinases}
#'   \item{apoptosis}{Programmed cell death. Contains extrinsic pathway, intrinsic
#'   pathway, executioner caspases, IAPs, p53 pathway, and ER stress genes}
#'   \item{autophagy}{Self-degradation process. Contains ULK complex, PI3K complex,
#'   ATG9 system, conjugation systems, LC3 family, receptors, mTOR/AMPK pathways}
#'   \item{panoptosis}{Combined cell death involving pyroptosis, apoptosis, and
#'   necroptosis. Contains PANoptosome components, effectors, and crosstalk factors}
#'   \item{netosis}{Neutrophil extracellular trap-associated cell death}
#'   \item{parthanatos}{PARP-dependent cell death}
#'   \item{entosis}{Cell-in-cell death}
#'   \item{oxeiptosis}{Oxidation-induced cell death}
#'   \item{alkaliptosis}{Alkaline pH-induced cell death}
#'   \item{ldcd}{Lysosome-dependent cell death}
#' }
#'
#' @details
#' Gene symbols follow HGNC nomenclature for Homo sapiens.
#'
#' @source
#' \itemize{
#'   \item FerrDb V3: \url{http://www.zhounan.org/ferrdb/}
#'   \item MSigDB: \url{https://www.gsea-msigdb.org/gsea/msigdb/}
#'   \item KEGG: \url{https://www.kegg.jp/}
#'   \item Primary literature (see package vignette for full references)
#' }
#'
#' @examples
#' # Load the gene sets
#' data(death_genesets)
#'
#' # View available death types
#' names(death_genesets)
#'
#' # Get ferroptosis drivers
#' death_genesets$ferroptosis$driver
#'
#' # Get all pyroptosis gasdermins
#' death_genesets$pyroptosis$gasdermins
#'
#' @seealso \code{\link{death_genesets_simple}} for simplified gene sets
#'
"death_genesets"


#' Simplified Cell Death Pathway Gene Sets
#'
#' A simplified version of cell death gene sets containing only core genes
#' for each pathway. This is useful for quick analyses and when computational
#' resources are limited.
#'
#' @format A named list containing 12 cell death types, each as a character
#' vector of core gene symbols:
#' \describe{
#'   \item{ferroptosis}{~45 core genes}
#'   \item{cuproptosis}{~20 core genes}
#'   \item{disulfidptosis}{~15 core genes}
#'   \item{pyroptosis}{~20 core genes}
#'   \item{necroptosis}{~15 core genes}
#'   \item{apoptosis}{~30 core genes}
#'   \item{autophagy}{~30 core genes}
#'   \item{panoptosis}{~15 core genes}
#'   \item{netosis}{~10 core genes}
#'   \item{parthanatos}{~6 core genes}
#'   \item{entosis}{~10 core genes}
#'   \item{oxeiptosis}{~3 core genes}
#' }
#'
#' @examples
#' # Load simplified gene sets
#' data(death_genesets_simple)
#'
#' # Get all ferroptosis core genes
#' death_genesets_simple$ferroptosis
#'
#' # Count genes in each pathway
#' sapply(death_genesets_simple, length)
#'
#' @seealso \code{\link{death_genesets}} for complete gene sets
#'
"death_genesets_simple"


#' Cell Death Pathway Information
#'
#' Metadata and descriptions for each cell death pathway.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{pathway}{Name of the cell death pathway}
#'   \item{full_name}{Full name in English}
#'   \item{chinese_name}{Chinese name}
#'   \item{year_discovered}{Year the pathway was first described}
#'   \item{key_paper}{Key reference paper}
#'   \item{description}{Brief description of the mechanism}
#'   \item{key_genes}{Comma-separated list of key genes}
#' }
#'
"death_pathway_info"


#' Example Expression Data
#'
#' A simulated gene expression matrix for demonstrating package functionality.
#' Contains 200 samples (100 normal, 100 tumor) and genes covering all cell
#' death pathways.
#'
#' @format A numeric matrix with genes as rows and samples as columns:
#' \describe{
#'   \item{Rows}{Gene symbols (HGNC format)}
#'   \item{Columns}{Sample IDs (Sample_001 to Sample_200)}
#'   \item{Values}{Log2-transformed expression values (simulated TPM)}
#' }
#'
#' @details
#' This dataset is designed to demonstrate all package functions. The expression
#' patterns are simulated to show differential expression of cell death genes
#' between normal and tumor samples:
#' \itemize{
#'   \item Samples 1-100: Normal tissue
#'   \item Samples 101-200: Tumor tissue
#'   \item Ferroptosis drivers are upregulated in tumors
#'   \item Ferroptosis suppressors are downregulated in tumors
#'   \item Pyroptosis genes are upregulated in tumors
#' }
#'
#' @examples
#' # Load example data
#' data(example_expr)
#'
#' # Check dimensions
#' dim(example_expr)
#'
#' # View first few genes and samples
#' example_expr[1:5, 1:5]
#'
#' # Calculate death scores
#' scores <- calculate_death_score(example_expr, method = "zscore")
#'
#' @seealso \code{\link{example_clinical}} for corresponding clinical data
#'
"example_expr"


#' Example Clinical Data
#'
#' Clinical information corresponding to the example expression data.
#' Contains sample annotations, survival data, and clinical features.
#'
#' @format A data frame with 200 rows (samples) and the following columns:
#' \describe{
#'   \item{sample_id}{Sample identifier matching column names of example_expr}
#'   \item{group}{Sample group: "Normal" or "Tumor"}
#'   \item{OS_time}{Overall survival time in months}
#'   \item{OS_status}{Overall survival status: 0 = alive/censored, 1 = deceased}
#'   \item{stage}{Tumor stage: I, II, III, IV (NA for normal samples)}
#'   \item{grade}{Tumor grade: G1, G2, G3 (NA for normal samples)}
#'   \item{age}{Patient age in years}
#'   \item{gender}{Patient gender: "Male" or "Female"}
#'   \item{treatment_response}{Treatment response: CR, PR, SD, PD (NA for normal)}
#' }
#'
#' @details
#' This dataset accompanies \code{\link{example_expr}} and provides clinical
#' annotations for survival analysis and group comparisons.
#'
#' Treatment response codes:
#' \itemize{
#'   \item CR: Complete Response
#'   \item PR: Partial Response
#'   \item SD: Stable Disease
#'   \item PD: Progressive Disease
#' }
#'
#' @examples
#' # Load clinical data
#' data(example_clinical)
#'
#' # View structure
#' str(example_clinical)
#'
#' # Sample distribution
#' table(example_clinical$group)
#'
#' # Survival status
#' table(example_clinical$OS_status)
#'
#' # Use with expression data for analysis
#' data(example_expr)
#' scores <- calculate_death_score(example_expr, method = "zscore")
#' plot_death_boxplot(scores, group = example_clinical$group)
#'
#' @seealso \code{\link{example_expr}} for corresponding expression data
#'
"example_clinical"
