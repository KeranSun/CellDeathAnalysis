#' Get Cell Death Gene Set
#'
#' Retrieve gene sets for specified cell death pathway(s).
#'
#' @param pathway Character vector of pathway names. Use "all" to get all pathways.
#'   Available options: "ferroptosis", "cuproptosis", "disulfidptosis", "pyroptosis",
#'   "necroptosis", "apoptosis", "autophagy", "panoptosis", "netosis", "parthanatos",
#'   "entosis", "oxeiptosis", "alkaliptosis", "ldcd".
#' @param type Character. Type of gene set to return:
#'   \itemize{
#'     \item "full": Return complete gene set with all sub-categories (default)
#'     \item "simple": Return simplified core genes only
#'     \item "all": Return all genes as a single vector (flattened)
#'   }
#' @param category Character. For "full" type, specify sub-category to return.
#'   E.g., "driver", "suppressor" for ferroptosis. Use NULL for all categories.
#'
#' @return A list of gene sets (if multiple pathways) or character vector (if single pathway
#'   with type="all" or specific category)
#'
#' @export
#'
#' @examples
#' # Get all ferroptosis genes
#' get_death_geneset("ferroptosis")
#'
#' # Get simplified gene sets for multiple pathways
#' get_death_geneset(c("ferroptosis", "pyroptosis"), type = "simple")
#'
#' # Get only ferroptosis drivers
#' get_death_geneset("ferroptosis", category = "driver")
#'
#' # Get all available pathways
#' get_death_geneset("all", type = "simple")
#'
get_death_geneset <- function(pathway = "all", 
                               type = c("full", "simple", "all"),
                               category = NULL) {
  

  type <- match.arg(type)
  
 # 加载数据
  if (type == "simple") {
    genesets <- death_genesets_simple
  } else {
    genesets <- death_genesets
  }
  
  # 获取可用通路
 available_pathways <- names(genesets)
  
  # 处理 "all" 参数
  if (length(pathway) == 1 && pathway == "all") {
    pathway <- available_pathways
  }
  
  # 检查通路名称
  invalid <- setdiff(pathway, available_pathways)
  if (length(invalid) > 0) {
    stop("Invalid pathway name(s): ", paste(invalid, collapse = ", "),
         "\nAvailable pathways: ", paste(available_pathways, collapse = ", "))
  }
  
  # 提取基因集
  result <- genesets[pathway]
  
  # 如果指定了category
  if (!is.null(category) && type == "full") {
    result <- lapply(result, function(x) {
      if (category %in% names(x)) {
        x[[category]]
      } else {
        warning("Category '", category, "' not found in some pathways")
        NULL
      }
    })
  }
  
  # 如果type是"all"，展平列表
  if (type == "all") {
    result <- lapply(result, function(x) {
      if (is.list(x)) {
        unique(unlist(x))
      } else {
        unique(x)
      }
    })
  }
  
  # 如果只有一个通路，返回向量而不是列表
  if (length(result) == 1 && (type == "all" || !is.null(category) || type == "simple")) {
    return(result[[1]])
  }
  
  return(result)
}


#' List Available Cell Death Pathways
#'
#' Display all available cell death pathways in the package.
#'
#' @param detailed Logical. If TRUE, return detailed information including
#'   Chinese names and gene counts. Default is FALSE.
#'
#' @return If detailed=FALSE, returns a character vector of pathway names.
#'   If detailed=TRUE, returns a data frame with pathway information.
#'
#' @export
#'
#' @examples
#' # Simple list
#' list_death_pathways()
#'
#' # Detailed information
#' list_death_pathways(detailed = TRUE)
#'
list_death_pathways <- function(detailed = FALSE) {
  
  if (!detailed) {
    return(names(death_genesets))
  }
  
  # 计算每个通路的基因数
  gene_counts <- sapply(death_genesets, function(x) {
    length(unique(unlist(x)))
  })
  
  simple_counts <- sapply(death_genesets_simple, function(x) {
    length(unique(unlist(x)))
  })
  
  # 合并信息
  info <- death_pathway_info[, c("pathway", "full_name", "chinese_name", 
                                  "year_discovered", "key_genes")]
  info$total_genes <- gene_counts[info$pathway]
  info$core_genes <- simple_counts[match(info$pathway, names(simple_counts))]
  
  return(info)
}


#' Get Pathway Information
#'
#' Retrieve detailed information about a specific cell death pathway.
#'
#' @param pathway Character. Name of the pathway.
#'
#' @return A list containing pathway information including description,
#'   key genes, morphology, and references.
#'
#' @export
#'
#' @examples
#' get_pathway_info("ferroptosis")
#' get_pathway_info("cuproptosis")
#'
get_pathway_info <- function(pathway) {
  
  if (!pathway %in% death_pathway_info$pathway) {
    stop("Unknown pathway: ", pathway,
         "\nAvailable pathways: ", paste(death_pathway_info$pathway, collapse = ", "))
  }
  
  info <- death_pathway_info[death_pathway_info$pathway == pathway, ]
  
  result <- list(
    pathway = info$pathway,
    full_name = info$full_name,
    chinese_name = info$chinese_name,
    year_discovered = info$year_discovered,
    key_paper = info$key_paper,
    description = info$description,
    key_genes = strsplit(info$key_genes, ", ")[[1]],
    morphology = info$morphology,
    gene_count = length(unique(unlist(death_genesets[[pathway]])))
  )
  
  class(result) <- c("death_pathway_info", "list")
  return(result)
}


#' Print Method for death_pathway_info
#' @param x A death_pathway_info object
#' @param ... Additional arguments (ignored)
#' @export
print.death_pathway_info <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(" ", x$full_name, " (", x$chinese_name, ")\n", sep = "")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("\n")
  cat("Year discovered:", x$year_discovered, "\n")
  cat("Key paper:", x$key_paper, "\n")
  cat("Total genes:", x$gene_count, "\n")
  cat("\n")
  cat("Description:\n")
  cat(strwrap(x$description, width = 60), sep = "\n")
  cat("\n")
  cat("Key genes:", paste(x$key_genes, collapse = ", "), "\n")
  cat("\n")
  cat("Morphology:\n")
  cat(strwrap(x$morphology, width = 60), sep = "\n")
  cat("\n")
  invisible(x)
}


#' Add Custom Gene Set
#'
#' Add a user-defined gene set to the analysis.
#'
#' @param name Character. Name for the custom gene set.
#' @param genes Character vector. Gene symbols to include.
#' @param description Character. Optional description of the gene set.
#'
#' @return Invisibly returns the updated gene set list.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Add a custom gene set
#' add_custom_geneset("my_pathway", c("GENE1", "GENE2", "GENE3"),
#'                    description = "My custom pathway")
#' }
#'
add_custom_geneset <- function(name, genes, description = NULL) {
  
  if (!is.character(genes)) {
    stop("genes must be a character vector")
  }
  
  # 创建自定义基因集环境（如果不存在）
  if (!exists("custom_genesets", envir = .GlobalEnv)) {
    assign("custom_genesets", list(), envir = .GlobalEnv)
  }
  
  custom <- get("custom_genesets", envir = .GlobalEnv)
  custom[[name]] <- list(
    genes = unique(genes),
    description = description,
    added_date = Sys.Date()
  )
  
  assign("custom_genesets", custom, envir = .GlobalEnv)
  
  message("Added custom gene set '", name, "' with ", length(genes), " genes")
  invisible(custom)
}


#' Get All Genes
#'
#' Get a union of all genes from specified pathways.
#'
#' @param pathways Character vector of pathway names. Default is "all".
#'
#' @return Character vector of unique gene symbols.
#'
#' @export
#'
#' @examples
#' # Get all cell death genes
#' all_genes <- get_all_death_genes()
#' length(all_genes)
#'
#' # Get genes from specific pathways
#' genes <- get_all_death_genes(c("ferroptosis", "pyroptosis"))
#'
get_all_death_genes <- function(pathways = "all") {
  
  genesets <- get_death_geneset(pathways, type = "all")
  
  if (is.list(genesets)) {
    return(unique(unlist(genesets)))
  } else {
    return(unique(genesets))
  }
}


#' Check Gene Overlap Between Pathways
#'
#' Calculate overlap between different cell death pathways.
#'
#' @param pathways Character vector of pathway names to compare.
#'   Default compares all pathways.
#'
#' @return A matrix showing the number of overlapping genes between pathways.
#'
#' @export
#'
#' @examples
#' # Check overlap between ferroptosis and autophagy
#' check_pathway_overlap(c("ferroptosis", "autophagy"))
#'
#' # Check all pathway overlaps
#' overlap_matrix <- check_pathway_overlap()
#'
check_pathway_overlap <- function(pathways = "all") {
  
  genesets <- get_death_geneset(pathways, type = "all")
  
  n <- length(genesets)
  pathway_names <- names(genesets)
  
  overlap_matrix <- matrix(0, n, n, 
                           dimnames = list(pathway_names, pathway_names))
  
  for (i in 1:n) {
    for (j in 1:n) {
      overlap_matrix[i, j] <- length(intersect(genesets[[i]], genesets[[j]]))
    }
  }
  
  return(overlap_matrix)
}
