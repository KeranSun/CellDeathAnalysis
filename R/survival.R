#' =============================================================================
#' Survival Analysis Module for CellDeathAnalysis
#' =============================================================================

#' Perform Survival Analysis for Cell Death Pathways
#'
#' Conduct Kaplan-Meier survival analysis and Cox regression for cell death
#' pathway scores.
#'
#' @param scores A death_scores object or data frame of pathway scores.
#' @param time Numeric vector of survival times.
#' @param status Numeric vector of event status (0 = censored, 1 = event).
#' @param pathway Character. Name of the pathway to analyze. If NULL, analyzes
#'   all pathways.
#' @param method Character. Classification method for grouping: "median", "mean",
#'   "quantile", or "optimal".
#' @param quantile Numeric. Quantile value when method = "quantile".
#'
#' @return A list containing:
#'   \itemize{
#'     \item km_fit: survfit object for Kaplan-Meier analysis
#'     \item cox_fit: coxph object for Cox regression
#'     \item log_rank_p: Log-rank test p-value
#'     \item cox_summary: Summary of Cox regression
#'     \item groups: Factor of high/low groups
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_expr)
#' data(example_clinical)
#'
#' scores <- calculate_death_score(example_expr, method = "zscore")
#'
#' # Analyze ferroptosis in tumor samples
#' tumor_idx <- example_clinical$group == "Tumor"
#' result <- death_survival(
#'   scores[tumor_idx, ],
#'   time = example_clinical$OS_time[tumor_idx],
#'   status = example_clinical$OS_status[tumor_idx],
#'   pathway = "ferroptosis"
#' )
#'
#' # View results
#' print(result)
#' }
#'
death_survival <- function(scores,
                            time,
                            status,
                            pathway = NULL,
                            method = c("median", "mean", "quantile", "optimal"),
                            quantile = 0.5) {
  
  # Check for survival package

if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }
  
  method <- match.arg(method)
  
  # If pathway is specified, extract that pathway
  if (!is.null(pathway)) {
    if (!pathway %in% names(scores)) {
      stop("Pathway '", pathway, "' not found in scores")
    }
    score_vec <- scores[[pathway]]
    pathway_name <- pathway
  } else {
    # Use first pathway if not specified
    score_vec <- scores[[1]]
    pathway_name <- names(scores)[1]
  }
  
  # Remove NA values
  complete_idx <- !is.na(score_vec) & !is.na(time) & !is.na(status)
  score_vec <- score_vec[complete_idx]
  time <- time[complete_idx]
  status <- status[complete_idx]
  
  # Classify samples
  if (method == "optimal" && requireNamespace("survminer", quietly = TRUE)) {
    # Use optimal cutpoint
    surv_data <- data.frame(time = time, status = status, score = score_vec)
    cutpoint <- survminer::surv_cutpoint(surv_data, time = "time", 
                                          event = "status", variables = "score")
    groups <- survminer::surv_categorize(cutpoint)$score
    cutoff <- cutpoint$cutpoint$cutpoint
  } else {
    groups <- classify_by_score(score_vec, method = method, quantile = quantile)
    cutoff <- attr(groups, "cutoff")
  }
  
  # Create survival data
  surv_obj <- survival::Surv(time, status)
  
  # Kaplan-Meier analysis
  km_fit <- survival::survfit(surv_obj ~ groups)
  
  # Log-rank test
  log_rank <- survival::survdiff(surv_obj ~ groups)
  log_rank_p <- 1 - pchisq(log_rank$chisq, df = length(unique(groups)) - 1)
  
  # Cox regression
  cox_fit <- survival::coxph(surv_obj ~ score_vec)
  cox_summary <- summary(cox_fit)
  
  # Compile results
  result <- list(
    pathway = pathway_name,
    km_fit = km_fit,
    cox_fit = cox_fit,
    log_rank_p = log_rank_p,
    cox_hr = cox_summary$conf.int[1, 1],
    cox_hr_lower = cox_summary$conf.int[1, 3],
    cox_hr_upper = cox_summary$conf.int[1, 4],
    cox_p = cox_summary$coefficients[1, 5],
    cutoff = cutoff,
    groups = groups,
    n_high = sum(groups == "High"),
    n_low = sum(groups == "Low"),
    method = method
  )
  
  class(result) <- c("death_survival", "list")
  
  return(result)
}


#' Print Method for death_survival
#' @param x A death_survival object
#' @param ... Additional arguments (ignored)
#' @export
print.death_survival <- function(x, ...) {
  cat("\n")
  cat("========================================\n")
  cat("  Survival Analysis:", x$pathway, "\n")
  cat("========================================\n\n")
  
  cat("Sample Size:\n")
  cat("  High group:", x$n_high, "\n")
  cat("  Low group:", x$n_low, "\n")
  cat("  Cutoff (", x$method, "):", round(x$cutoff, 3), "\n\n")
  
  cat("Log-rank Test:\n")
  cat("  p-value:", format.pval(x$log_rank_p, digits = 3), "\n\n")
  
  cat("Cox Regression:\n")
  cat("  Hazard Ratio:", round(x$cox_hr, 3), "\n")
  cat("  95% CI: [", round(x$cox_hr_lower, 3), ", ", round(x$cox_hr_upper, 3), "]\n")
  cat("  p-value:", format.pval(x$cox_p, digits = 3), "\n")
  cat("\n")
  
  invisible(x)
}


#' Plot Survival Curves for Cell Death Pathways
#'
#' Create Kaplan-Meier survival curves with risk table.
#'
#' @param surv_result A death_survival object from death_survival().
#' @param palette Character vector of colors for groups.
#' @param risk_table Logical. Whether to show risk table. Default is TRUE.
#' @param pval Logical. Whether to show p-value. Default is TRUE.
#' @param conf_int Logical. Whether to show confidence intervals. Default is FALSE.
#' @param title Character. Plot title. If NULL, uses pathway name.
#' @param xlab Character. X-axis label. Default is "Time".
#' @param ylab Character. Y-axis label. Default is "Survival Probability".
#' @param legend_title Character. Legend title. Default is "Group".
#' @param ... Additional arguments passed to ggsurvplot.
#'
#' @return A ggsurvplot object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- death_survival(scores, time, status, pathway = "ferroptosis")
#' plot_survival(result)
#' }
#'
plot_survival <- function(surv_result,
                           palette = c("#E41A1C", "#377EB8"),
                           risk_table = TRUE,
                           pval = TRUE,
                           conf_int = FALSE,
                           title = NULL,
                           xlab = "Time",
                           ylab = "Survival Probability",
                           legend_title = "Group",
                           ...) {
  
  if (!inherits(surv_result, "death_survival")) {
    stop("surv_result must be a death_survival object")
  }
  
  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop("Package 'survminer' is required. Install with: install.packages('survminer')")
  }
  
  if (is.null(title)) {
    title <- paste("Survival Analysis:", surv_result$pathway)
  }
  
  # Create plot
  p <- survminer::ggsurvplot(
    surv_result$km_fit,
    data = data.frame(groups = surv_result$groups),
    palette = palette,
    risk.table = risk_table,
    pval = pval,
    conf.int = conf_int,
    title = title,
    xlab = xlab,
    ylab = ylab,
    legend.title = legend_title,
    legend.labs = c("Low", "High"),
    ggtheme = ggplot2::theme_bw(),
    risk.table.height = 0.25,
    ...
  )
  
  return(p)
}


#' Batch Survival Analysis for Multiple Pathways
#'
#' Perform survival analysis for all cell death pathways.
#'
#' @param scores A death_scores data frame.
#' @param time Numeric vector of survival times.
#' @param status Numeric vector of event status.
#' @param method Classification method.
#' @param p_adjust Method for p-value adjustment. Default is "BH".
#'
#' @return A data frame with survival analysis results for each pathway.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- batch_survival(scores, time, status)
#' print(results)
#' }
#'
batch_survival <- function(scores,
                            time,
                            status,
                            method = "median",
                            p_adjust = "BH") {
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required")
  }
  
  pathways <- names(scores)
  results <- list()
  
  for (pathway in pathways) {
    tryCatch({
      res <- death_survival(scores, time, status, 
                            pathway = pathway, method = method)
      results[[pathway]] <- data.frame(
        pathway = pathway,
        n_high = res$n_high,
        n_low = res$n_low,
        log_rank_p = res$log_rank_p,
        cox_hr = res$cox_hr,
        cox_hr_lower = res$cox_hr_lower,
        cox_hr_upper = res$cox_hr_upper,
        cox_p = res$cox_p,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning("Error analyzing ", pathway, ": ", e$message)
    })
  }
  
  # Combine results
  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  
  # Adjust p-values
  result_df$log_rank_p_adj <- p.adjust(result_df$log_rank_p, method = p_adjust)
  result_df$cox_p_adj <- p.adjust(result_df$cox_p, method = p_adjust)
  
  # Sort by p-value
  result_df <- result_df[order(result_df$log_rank_p), ]
  
  # Add significance
  result_df$significance <- ifelse(result_df$log_rank_p_adj < 0.001, "***",
                                    ifelse(result_df$log_rank_p_adj < 0.01, "**",
                                           ifelse(result_df$log_rank_p_adj < 0.05, "*", "")))
  
  class(result_df) <- c("batch_survival", "data.frame")
  
  return(result_df)
}


#' Print Method for batch_survival
#' @param x A batch_survival object
#' @param ... Additional arguments
#' @export
print.batch_survival <- function(x, ...) {
  cat("\n")
  cat("========================================\n")
  cat("  Batch Survival Analysis Results\n")
  cat("========================================\n\n")
  
  cat("Pathways analyzed:", nrow(x), "\n")
  cat("Significant (adj.p < 0.05):", sum(x$log_rank_p_adj < 0.05), "\n\n")
  
  # Print top results
  print_df <- x[, c("pathway", "cox_hr", "log_rank_p", "log_rank_p_adj", "significance")]
  print_df$cox_hr <- round(print_df$cox_hr, 3)
  print_df$log_rank_p <- format.pval(print_df$log_rank_p, digits = 2)
  print_df$log_rank_p_adj <- format.pval(print_df$log_rank_p_adj, digits = 2)
  
  print(print_df, row.names = FALSE)
  cat("\n")
  
  invisible(x)
}


#' Plot Forest Plot for Batch Survival Results
#'
#' Create a forest plot showing hazard ratios for all pathways.
#'
#' @param batch_result A batch_survival object.
#' @param show_p Logical. Whether to show p-values. Default is TRUE.
#' @param order_by Character. Order pathways by "hr", "p", or "name".
#' @param colors Named vector of colors for significant/non-significant.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
plot_forest <- function(batch_result,
                         show_p = TRUE,
                         order_by = c("p", "hr", "name"),
                         colors = c("significant" = "#E41A1C", "ns" = "grey50")) {
  
  order_by <- match.arg(order_by)
  
  # Prepare data
  plot_data <- batch_result
  plot_data$significant <- ifelse(plot_data$log_rank_p_adj < 0.05, "significant", "ns")
  
  # Order
  if (order_by == "p") {
    plot_data <- plot_data[order(plot_data$log_rank_p), ]
  } else if (order_by == "hr") {
    plot_data <- plot_data[order(plot_data$cox_hr, decreasing = TRUE), ]
  }
  
  plot_data$pathway <- factor(plot_data$pathway, levels = rev(plot_data$pathway))
  
  # Create plot
  p <- ggplot(plot_data, aes(x = cox_hr, y = pathway, color = significant)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = cox_hr_lower, xmax = cox_hr_upper), 
                   height = 0.2, linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = colors, guide = "none") +
    scale_x_log10() +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 10)
    ) +
    labs(
      x = "Hazard Ratio (log scale)",
      y = NULL,
      title = "Forest Plot: Cell Death Pathways and Survival"
    )
  
  # Add p-values
  if (show_p) {
    plot_data$p_label <- paste0("p=", format.pval(plot_data$log_rank_p_adj, digits = 2))
    p <- p + 
      geom_text(aes(x = max(cox_hr_upper) * 1.2, label = p_label),
                hjust = 0, size = 3, color = "black")
  }
  
  return(p)
}


#' Multivariate Cox Regression with Clinical Variables
#'
#' Perform multivariate Cox regression including death scores and clinical variables.
#'
#' @param scores A death_scores object or numeric vector.
#' @param time Survival time.
#' @param status Event status.
#' @param clinical_vars Data frame of clinical variables to include.
#' @param pathway Character. Pathway name if scores is a data frame.
#'
#' @return A list containing cox model and summary.
#'
#' @export
#'
death_cox_multivariate <- function(scores,
                                    time,
                                    status,
                                    clinical_vars,
                                    pathway = NULL) {
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required")
  }
  
  # Extract score
  if (is.data.frame(scores) && !is.null(pathway)) {
    score_vec <- scores[[pathway]]
    pathway_name <- pathway
  } else if (is.numeric(scores)) {
    score_vec <- scores
    pathway_name <- "death_score"
  } else {
    score_vec <- scores[[1]]
    pathway_name <- names(scores)[1]
  }
  
  # Combine data
  cox_data <- data.frame(
    time = time,
    status = status,
    death_score = score_vec
  )
  cox_data <- cbind(cox_data, clinical_vars)
  
  # Remove NA
  cox_data <- na.omit(cox_data)
  
  # Build formula
  var_names <- c("death_score", names(clinical_vars))
  formula_str <- paste("survival::Surv(time, status) ~", paste(var_names, collapse = " + "))
  
  # Fit model
  cox_fit <- survival::coxph(as.formula(formula_str), data = cox_data)
  
  result <- list(
    pathway = pathway_name,
    model = cox_fit,
    summary = summary(cox_fit),
    n = nrow(cox_data),
    variables = var_names
  )
  
  class(result) <- c("death_cox_multi", "list")
  
  return(result)
}


#' Print Method for death_cox_multi
#' @param x A death_cox_multi object
#' @param ... Additional arguments
#' @export
print.death_cox_multi <- function(x, ...) {
  cat("\n")
  cat("========================================\n")
  cat("  Multivariate Cox Regression\n")
  cat("========================================\n\n")
  
  cat("Pathway:", x$pathway, "\n")
  cat("Sample size:", x$n, "\n")
  cat("Variables:", paste(x$variables, collapse = ", "), "\n\n")
  
  # Print coefficients
  coef_df <- as.data.frame(x$summary$coefficients)
  coef_df$HR <- round(exp(coef_df$coef), 3)
  coef_df$`p-value` <- format.pval(coef_df$`Pr(>|z|)`, digits = 3)
  
  print(coef_df[, c("HR", "p-value")])
  cat("\n")
  
  cat("Concordance:", round(x$summary$concordance[1], 3), "\n")
  cat("Likelihood ratio test p:", format.pval(x$summary$logtest[3], digits = 3), "\n")
  cat("\n")
  
  invisible(x)
}


#' Time-dependent ROC Analysis
#'
#' Perform time-dependent ROC analysis for death scores.
#'
#' @param scores Numeric vector of scores.
#' @param time Survival time.
#' @param status Event status.
#' @param times Numeric vector of time points for ROC evaluation.
#' @param pathway Character. Pathway name for labeling.
#'
#' @return A list containing ROC curves and AUC values.
#'
#' @export
#'
death_time_roc <- function(scores,
                            time,
                            status,
                            times = c(12, 36, 60),
                            pathway = "Death Score") {
  
  if (!requireNamespace("timeROC", quietly = TRUE)) {
    stop("Package 'timeROC' is required. Install with: install.packages('timeROC')")
  }
  
  # Remove NA
  complete_idx <- !is.na(scores) & !is.na(time) & !is.na(status)
  scores <- scores[complete_idx]
  time <- time[complete_idx]
  status <- status[complete_idx]
  
  # Compute time-dependent ROC
  roc_result <- timeROC::timeROC(
    T = time,
    delta = status,
    marker = scores,
    cause = 1,
    times = times,
    iid = TRUE
  )
  
  result <- list(
    pathway = pathway,
    roc = roc_result,
    times = times,
    auc = roc_result$AUC,
    auc_ci = confint(roc_result)
  )
  
  class(result) <- c("death_time_roc", "list")
  
  return(result)
}


#' Plot Time-dependent ROC Curves
#'
#' @param roc_result A death_time_roc object.
#' @param colors Colors for different time points.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
plot_time_roc <- function(roc_result, colors = NULL) {
  
  if (!inherits(roc_result, "death_time_roc")) {
    stop("roc_result must be a death_time_roc object")
  }
  
  roc <- roc_result$roc
  times <- roc_result$times
  
  if (is.null(colors)) {
    colors <- scales::hue_pal()(length(times))
  }
  
  # Prepare data for plotting
  plot_data <- data.frame()
  for (i in seq_along(times)) {
    t <- times[i]
    fpr <- roc$FP[, i]
    tpr <- roc$TP[, i]
    auc <- round(roc$AUC[i], 3)
    
    temp <- data.frame(
      FPR = fpr,
      TPR = tpr,
      Time = paste0(t, " months (AUC=", auc, ")")
    )
    plot_data <- rbind(plot_data, temp)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = FPR, y = TPR, color = Time)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = colors) +
    coord_equal() +
    theme_bw() +
    theme(
      legend.position = c(0.7, 0.2),
      legend.background = element_rect(fill = "white", color = "grey80")
    ) +
    labs(
      x = "1 - Specificity (FPR)",
      y = "Sensitivity (TPR)",
      title = paste("Time-dependent ROC:", roc_result$pathway),
      color = "Time Point"
    )
  
  return(p)
}
