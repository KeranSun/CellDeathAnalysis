#' =============================================================================
#' Machine Learning Module for CellDeathAnalysis
#' =============================================================================

#' Build Predictive Model Using Death Pathway Scores
#'
#' Train machine learning models to predict clinical outcomes using cell death
#' pathway scores.
#'
#' @param scores A death_scores data frame or matrix.
#' @param outcome Outcome variable (factor for classification, numeric for regression).
#' @param method Model type: "rf" (Random Forest), "lasso", "ridge", "elasticnet",
#'   "svm", "xgboost", or "ensemble".
#' @param train_ratio Proportion of data for training. Default is 0.7.
#' @param nfolds Number of folds for cross-validation. Default is 5.
#' @param seed Random seed for reproducibility.
#' @param ... Additional arguments passed to the model function.
#'
#' @return A death_model object containing the trained model and performance metrics.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_expr)
#' data(example_clinical)
#'
#' scores <- calculate_death_score(example_expr, method = "zscore")
#' outcome <- factor(example_clinical$group)
#'
#' model <- death_build_model(scores, outcome, method = "rf")
#' print(model)
#' }
#'
death_build_model <- function(scores,
                               outcome,
                               method = c("rf", "lasso", "ridge", "elasticnet",
                                          "svm", "xgboost", "ensemble"),
                               train_ratio = 0.7,
                               nfolds = 5,
                               seed = 123,
                               ...) {
 
  method <- match.arg(method)
  set.seed(seed)
 
  # Prepare data
  scores <- as.data.frame(scores)
 
  # Remove NA
  complete_idx <- complete.cases(scores) & !is.na(outcome)
  scores <- scores[complete_idx, ]
  outcome <- outcome[complete_idx]
 
  # Determine task type
  if (is.factor(outcome) || is.character(outcome)) {
    task <- "classification"
    outcome <- as.factor(outcome)
  } else {
    task <- "regression"
  }
 
  # Split data
  n <- nrow(scores)
  train_idx <- sample(1:n, size = round(n * train_ratio))
  test_idx <- setdiff(1:n, train_idx)
 
  X_train <- scores[train_idx, ]
  X_test <- scores[test_idx, ]
  y_train <- outcome[train_idx]
  y_test <- outcome[test_idx]
 
  # Train model
  if (method == "rf") {
    model_result <- .train_rf(X_train, y_train, X_test, y_test, task, ...)
   
  } else if (method %in% c("lasso", "ridge", "elasticnet")) {
    model_result <- .train_glmnet(X_train, y_train, X_test, y_test, task, method, nfolds, ...)
   
  } else if (method == "svm") {
    model_result <- .train_svm(X_train, y_train, X_test, y_test, task, ...)
   
  } else if (method == "xgboost") {
    model_result <- .train_xgboost(X_train, y_train, X_test, y_test, task, nfolds, ...)
   
  } else if (method == "ensemble") {
    model_result <- .train_ensemble(X_train, y_train, X_test, y_test, task, nfolds, ...)
  }
 
  # Create result object
  result <- list(
    model = model_result$model,
    method = method,
    task = task,
    predictions = model_result$predictions,
    performance = model_result$performance,
    importance = model_result$importance,
    train_idx = train_idx,
    test_idx = test_idx,
    feature_names = colnames(scores),
    n_train = length(train_idx),
    n_test = length(test_idx),
    seed = seed
  )
 
  class(result) <- c("death_model", "list")
 
  return(result)
}


#' Internal: Train Random Forest
#' @keywords internal
.train_rf <- function(X_train, y_train, X_test, y_test, task, ...) {
 
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' is required. Install with: install.packages('randomForest')")
  }
 
  # Train model
  model <- randomForest::randomForest(
    x = X_train,
    y = y_train,
    importance = TRUE,
    ...
  )
 
  # Predict
  predictions <- predict(model, X_test)
 
  # Calculate performance
  if (task == "classification") {
    pred_prob <- predict(model, X_test, type = "prob")
    performance <- .calc_class_performance(y_test, predictions, pred_prob)
  } else {
    performance <- .calc_reg_performance(y_test, predictions)
  }
 
  # Get importance
  importance <- randomForest::importance(model)
  importance <- data.frame(
    feature = rownames(importance),
    importance = importance[, 1],
    stringsAsFactors = FALSE
  )
  importance <- importance[order(-importance$importance), ]
 
  return(list(
    model = model,
    predictions = predictions,
    performance = performance,
    importance = importance
  ))
}


#' Internal: Train glmnet (Lasso/Ridge/Elastic Net)
#' @keywords internal
.train_glmnet <- function(X_train, y_train, X_test, y_test, task, method, nfolds, ...) {
 
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required. Install with: install.packages('glmnet')")
  }
 
  # Set alpha based on method
  alpha <- switch(method,
                   "lasso" = 1,
                   "ridge" = 0,
                   "elasticnet" = 0.5)
 
  # Convert to matrix
  X_train_mat <- as.matrix(X_train)
  X_test_mat <- as.matrix(X_test)
 
  # Set family
  if (task == "classification") {
    family <- "binomial"
    y_train_num <- as.numeric(y_train) - 1
    y_test_num <- as.numeric(y_test) - 1
  } else {
    family <- "gaussian"
    y_train_num <- y_train
    y_test_num <- y_test
  }
 
  # Cross-validation
  cv_model <- glmnet::cv.glmnet(
    x = X_train_mat,
    y = y_train_num,
    alpha = alpha,
    family = family,
    nfolds = nfolds,
    ...
  )
 
  # Best model
  model <- cv_model
 
  # Predict
  if (task == "classification") {
    pred_prob <- predict(model, X_test_mat, s = "lambda.min", type = "response")
    predictions <- factor(ifelse(pred_prob > 0.5, levels(y_test)[2], levels(y_test)[1]),
                          levels = levels(y_test))
    pred_prob_mat <- cbind(1 - pred_prob, pred_prob)
    colnames(pred_prob_mat) <- levels(y_test)
    performance <- .calc_class_performance(y_test, predictions, pred_prob_mat)
  } else {
    predictions <- as.vector(predict(model, X_test_mat, s = "lambda.min"))
    pred_prob_mat <- NULL
    performance <- .calc_reg_performance(y_test, predictions)
  }
 
  # Get importance (coefficients)
  coef_mat <- as.matrix(coef(model, s = "lambda.min"))
  importance <- data.frame(
    feature = rownames(coef_mat)[-1],
    importance = abs(coef_mat[-1, 1]),
    coefficient = coef_mat[-1, 1],
    stringsAsFactors = FALSE
  )
  importance <- importance[order(-importance$importance), ]
 
  return(list(
    model = model,
    predictions = predictions,
    performance = performance,
    importance = importance
  ))
}


#' Internal: Train SVM
#' @keywords internal
.train_svm <- function(X_train, y_train, X_test, y_test, task, ...) {
 
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("Package 'e1071' is required. Install with: install.packages('e1071')")
  }
 
  # Train model
  model <- e1071::svm(
    x = X_train,
    y = y_train,
    probability = TRUE,
    ...
  )
 
  # Predict
  predictions <- predict(model, X_test, probability = TRUE)
 
  # Calculate performance
  if (task == "classification") {
    pred_prob <- attr(predictions, "probabilities")
    performance <- .calc_class_performance(y_test, predictions, pred_prob)
  } else {
    performance <- .calc_reg_performance(y_test, as.numeric(predictions))
  }
 
  # SVM doesn't have built-in importance
  importance <- data.frame(
    feature = colnames(X_train),
    importance = NA,
    stringsAsFactors = FALSE
  )
 
  return(list(
    model = model,
    predictions = predictions,
    performance = performance,
    importance = importance
  ))
}


#' Internal: Train XGBoost
#' @keywords internal
.train_xgboost <- function(X_train, y_train, X_test, y_test, task, nfolds, ...) {
 
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required. Install with: install.packages('xgboost')")
  }
 
  X_train_mat <- as.matrix(X_train)
  X_test_mat <- as.matrix(X_test)
 
  if (task == "classification") {
    y_train_num <- as.numeric(y_train) - 1
    y_test_num <- as.numeric(y_test) - 1
    objective <- "binary:logistic"
    eval_metric <- "auc"
  } else {
    y_train_num <- y_train
    y_test_num <- y_test
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
  }
 
  # Create DMatrix
  dtrain <- xgboost::xgb.DMatrix(data = X_train_mat, label = y_train_num)
  dtest <- xgboost::xgb.DMatrix(data = X_test_mat, label = y_test_num)
 
  # Cross-validation for best nrounds
  cv_result <- xgboost::xgb.cv(
    data = dtrain,
    nrounds = 500,
    nfold = nfolds,
    objective = objective,
    eval_metric = eval_metric,
    early_stopping_rounds = 20,
    verbose = 0
  )
 
  best_nrounds <- cv_result$best_iteration
 
  # Train final model
  model <- xgboost::xgb.train(
    data = dtrain,
    nrounds = best_nrounds,
    objective = objective,
    eval_metric = eval_metric,
    verbose = 0,
    ...
  )
 
  # Predict
  pred_raw <- predict(model, dtest)
 
  if (task == "classification") {
    predictions <- factor(ifelse(pred_raw > 0.5, levels(y_test)[2], levels(y_test)[1]),
                          levels = levels(y_test))
    pred_prob <- cbind(1 - pred_raw, pred_raw)
    colnames(pred_prob) <- levels(y_test)
    performance <- .calc_class_performance(y_test, predictions, pred_prob)
  } else {
    predictions <- pred_raw
    pred_prob <- NULL
    performance <- .calc_reg_performance(y_test, predictions)
  }
 
  # Get importance
  imp <- xgboost::xgb.importance(model = model)
  importance <- data.frame(
    feature = imp$Feature,
    importance = imp$Gain,
    stringsAsFactors = FALSE
  )
 
  return(list(
    model = model,
    predictions = predictions,
    performance = performance,
    importance = importance
  ))
}


#' Internal: Train Ensemble
#' @keywords internal
.train_ensemble <- function(X_train, y_train, X_test, y_test, task, nfolds, ...) {
 
  # Train multiple models
  rf_result <- tryCatch(
    .train_rf(X_train, y_train, X_test, y_test, task),
    error = function(e) NULL
  )
 
  lasso_result <- tryCatch(
    .train_glmnet(X_train, y_train, X_test, y_test, task, "lasso", nfolds),
    error = function(e) NULL
  )
 
  xgb_result <- tryCatch(
    .train_xgboost(X_train, y_train, X_test, y_test, task, nfolds),
    error = function(e) NULL
  )
 
  # Combine predictions (majority voting for classification, average for regression)
  models <- list(rf = rf_result, lasso = lasso_result, xgb = xgb_result)
  models <- models[!sapply(models, is.null)]
 
  if (length(models) == 0) {
    stop("No models could be trained")
  }
 
  if (task == "classification") {
    # Majority voting
    pred_matrix <- sapply(models, function(m) as.character(m$predictions))
    predictions <- apply(pred_matrix, 1, function(x) {
      tab <- table(x)
      names(tab)[which.max(tab)]
    })
    predictions <- factor(predictions, levels = levels(y_test))
   
    # Average probabilities (simplified)
    performance <- .calc_class_performance(y_test, predictions, NULL)
   
  } else {
    # Average predictions
    pred_matrix <- sapply(models, function(m) m$predictions)
    predictions <- rowMeans(pred_matrix)
    performance <- .calc_reg_performance(y_test, predictions)
  }
 
  # Combine importance
  importance <- do.call(rbind, lapply(names(models), function(name) {
    imp <- models[[name]]$importance
    imp$model <- name
    imp
  }))
 
  # Aggregate importance
  agg_imp <- aggregate(importance ~ feature, data = importance, FUN = mean, na.rm = TRUE)
  agg_imp <- agg_imp[order(-agg_imp$importance), ]
 
  return(list(
    model = models,
    predictions = predictions,
    performance = performance,
    importance = agg_imp
  ))
}


#' Internal: Calculate classification performance
#' @keywords internal
.calc_class_performance <- function(actual, predicted, prob_matrix = NULL) {
 
  # Confusion matrix
  conf_mat <- table(Actual = actual, Predicted = predicted)
 
  # Accuracy
  accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
 
  # For binary classification
  if (length(levels(actual)) == 2) {
    tp <- conf_mat[2, 2]
    tn <- conf_mat[1, 1]
    fp <- conf_mat[1, 2]
    fn <- conf_mat[2, 1]
   
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    precision <- tp / (tp + fp)
    f1 <- 2 * precision * sensitivity / (precision + sensitivity)
   
    # AUC
    if (!is.null(prob_matrix) && requireNamespace("pROC", quietly = TRUE)) {
      roc_obj <- pROC::roc(actual, prob_matrix[, 2], quiet = TRUE)
      auc <- as.numeric(pROC::auc(roc_obj))
    } else {
      auc <- NA
    }
   
    return(list(
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      precision = precision,
      f1 = f1,
      auc = auc,
      confusion_matrix = conf_mat
    ))
  }
 
  return(list(
    accuracy = accuracy,
    confusion_matrix = conf_mat
  ))
}


#' Internal: Calculate regression performance
#' @keywords internal
.calc_reg_performance <- function(actual, predicted) {
 
  # RMSE
  rmse <- sqrt(mean((actual - predicted)^2))
 
  # MAE
  mae <- mean(abs(actual - predicted))
 
  # R-squared
  ss_res <- sum((actual - predicted)^2)
  ss_tot <- sum((actual - mean(actual))^2)
  r2 <- 1 - ss_res / ss_tot
 
  # Correlation
  cor_val <- cor(actual, predicted)
 
  return(list(
    rmse = rmse,
    mae = mae,
    r2 = r2,
    correlation = cor_val
  ))
}


#' Print Method for death_model
#' @param x A death_model object
#' @param ... Additional arguments
#' @export
print.death_model <- function(x, ...) {
  cat("\n")
  cat("========================================\n")
  cat("  Cell Death Prediction Model\n")
  cat("========================================\n\n")
 
  cat("Method:", toupper(x$method), "\n")
  cat("Task:", x$task, "\n")
  cat("Features:", length(x$feature_names), "\n")
  cat("Training samples:", x$n_train, "\n")
  cat("Test samples:", x$n_test, "\n\n")
 
  cat("Performance:\n")
  if (x$task == "classification") {
    cat("  Accuracy:", round(x$performance$accuracy, 4), "\n")
    if (!is.null(x$performance$auc)) {
      cat("  AUC:", round(x$performance$auc, 4), "\n")
    }
    cat("  Sensitivity:", round(x$performance$sensitivity, 4), "\n")
    cat("  Specificity:", round(x$performance$specificity, 4), "\n")
    cat("  F1 Score:", round(x$performance$f1, 4), "\n")
  } else {
    cat("  R-squared:", round(x$performance$r2, 4), "\n")
    cat("  RMSE:", round(x$performance$rmse, 4), "\n")
    cat("  Correlation:", round(x$performance$correlation, 4), "\n")
  }
 
  cat("\nTop 5 Important Features:\n")
  print(head(x$importance, 5), row.names = FALSE)
  cat("\n")
 
  invisible(x)
}


#' Plot Model Performance
#'
#' Visualize model performance metrics and feature importance.
#'
#' @param model A death_model object.
#' @param type Plot type: "importance", "roc", "confusion", or "all".
#'
#' @return A ggplot object or list of plots.
#'
#' @import ggplot2
#' @export
#'
plot_model <- function(model, type = c("importance", "roc", "confusion", "all")) {
 
  type <- match.arg(type)
 
  plots <- list()
 
  # Feature importance plot
  if (type %in% c("importance", "all")) {
    imp_data <- head(model$importance, 10)
    imp_data$feature <- factor(imp_data$feature, levels = rev(imp_data$feature))
   
    p_imp <- ggplot(imp_data, aes(x = importance, y = feature)) +
      geom_col(fill = "steelblue") +
      theme_bw() +
      labs(
        x = "Importance",
        y = NULL,
        title = "Top 10 Important Features"
      )
    plots$importance <- p_imp
  }
 
  # ROC curve (classification only)
  if (type %in% c("roc", "all") && model$task == "classification") {
    if (!is.na(model$performance$auc)) {
      # Would need actual probabilities stored
      # Placeholder message
      message("ROC curve requires probability predictions to be stored")
    }
  }
 
  # Confusion matrix
  if (type %in% c("confusion", "all") && model$task == "classification") {
    conf_mat <- model$performance$confusion_matrix
    conf_df <- as.data.frame(conf_mat)
   
    p_conf <- ggplot(conf_df, aes(x = Predicted, y = Actual, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Freq), size = 6) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      theme_minimal() +
      labs(
        title = "Confusion Matrix",
        fill = "Count"
      )
    plots$confusion <- p_conf
  }
 
  if (length(plots) == 1) {
    return(plots[[1]])
  }
 
  return(plots)
}


#' Predict Using Trained Model
#'
#' Make predictions on new data using a trained death_model.
#'
#' @param model A death_model object.
#' @param newdata New data (expression matrix or death scores).
#' @param type For classification: "class" or "prob".
#'
#' @return Predictions.
#'
#' @export
#'
death_predict <- function(model, newdata, type = c("class", "prob")) {
 
  type <- match.arg(type)
 
  # Ensure newdata has same features
  newdata <- as.data.frame(newdata)
  missing_features <- setdiff(model$feature_names, colnames(newdata))
 
  if (length(missing_features) > 0) {
    warning("Missing features: ", paste(missing_features, collapse = ", "))
    # Add missing features as 0
    for (f in missing_features) {
      newdata[[f]] <- 0
    }
  }
 
  newdata <- newdata[, model$feature_names, drop = FALSE]
 
  # Predict based on model type
  if (model$method == "rf") {
    if (type == "prob" && model$task == "classification") {
      pred <- predict(model$model, newdata, type = "prob")
    } else {
      pred <- predict(model$model, newdata)
    }
   
  } else if (model$method %in% c("lasso", "ridge", "elasticnet")) {
    newdata_mat <- as.matrix(newdata)
    if (type == "prob" && model$task == "classification") {
      pred <- predict(model$model, newdata_mat, s = "lambda.min", type = "response")
    } else {
      pred_raw <- predict(model$model, newdata_mat, s = "lambda.min")
      if (model$task == "classification") {
        pred <- factor(ifelse(pred_raw > 0.5, "High", "Low"))
      } else {
        pred <- as.vector(pred_raw)
      }
    }
   
  } else if (model$method == "svm") {
    pred <- predict(model$model, newdata, probability = (type == "prob"))
   
  } else if (model$method == "xgboost") {
    newdata_mat <- as.matrix(newdata)
    dmat <- xgboost::xgb.DMatrix(data = newdata_mat)
    pred_raw <- predict(model$model, dmat)
   
    if (model$task == "classification") {
      if (type == "prob") {
        pred <- pred_raw
      } else {
        pred <- factor(ifelse(pred_raw > 0.5, "High", "Low"))
      }
    } else {
      pred <- pred_raw
    }
   
  } else if (model$method == "ensemble") {
    # Average predictions from all models
    preds <- lapply(model$model, function(m) {
      death_predict(list(model = m$model, method = names(model$model)[1],
                         task = model$task, feature_names = model$feature_names),
                    newdata, type = "class")
    })
    # Majority voting
    pred_mat <- do.call(cbind, preds)
    pred <- apply(pred_mat, 1, function(x) {
      tab <- table(x)
      names(tab)[which.max(tab)]
    })
    pred <- factor(pred)
  }
 
  return(pred)
}


#' Cross-Validation for Model Evaluation
#'
#' Perform k-fold cross-validation for model evaluation.
#'
#' @param scores A death_scores data frame.
#' @param outcome Outcome variable.
#' @param method Model method.
#' @param nfolds Number of folds. Default is 5.
#' @param seed Random seed.
#'
#' @return A list with cross-validation results.
#'
#' @export
#'
death_cv <- function(scores,
                      outcome,
                      method = "rf",
                      nfolds = 5,
                      seed = 123) {
 
  set.seed(seed)
 
  scores <- as.data.frame(scores)
  n <- nrow(scores)
 
  # Create folds
  folds <- sample(rep(1:nfolds, length.out = n))
 
  # Determine task
  if (is.factor(outcome) || is.character(outcome)) {
    task <- "classification"
    outcome <- as.factor(outcome)
  } else {
    task <- "regression"
  }
 
  # Store results
  cv_results <- list()
 
  for (k in 1:nfolds) {
    # Split
    test_idx <- which(folds == k)
    train_idx <- which(folds != k)
   
    X_train <- scores[train_idx, ]
    X_test <- scores[test_idx, ]
    y_train <- outcome[train_idx]
    y_test <- outcome[test_idx]
   
    # Train and evaluate
    tryCatch({
      if (method == "rf") {
        result <- .train_rf(X_train, y_train, X_test, y_test, task)
      } else if (method %in% c("lasso", "ridge", "elasticnet")) {
        result <- .train_glmnet(X_train, y_train, X_test, y_test, task, method, 5)
      } else if (method == "xgboost") {
        result <- .train_xgboost(X_train, y_train, X_test, y_test, task, 5)
      }
     
      cv_results[[k]] <- result$performance
    }, error = function(e) {
      warning("Fold ", k, " failed: ", e$message)
    })
  }
 
  # Summarize
  if (task == "classification") {
    metrics <- c("accuracy", "auc", "sensitivity", "specificity", "f1")
  } else {
    metrics <- c("rmse", "r2", "correlation")
  }
 
  summary_df <- data.frame(
    metric = metrics,
    mean = sapply(metrics, function(m) {
      vals <- sapply(cv_results, function(r) r[[m]])
      mean(vals, na.rm = TRUE)
    }),
    sd = sapply(metrics, function(m) {
      vals <- sapply(cv_results, function(r) r[[m]])
      sd(vals, na.rm = TRUE)
    }),
    stringsAsFactors = FALSE
  )
 
  return(list(
    method = method,
    nfolds = nfolds,
    summary = summary_df,
    fold_results = cv_results
  ))
}


#' Build Survival Prediction Model
#'
#' Build a Cox-LASSO or Random Survival Forest model for survival prediction.
#'
#' @param scores A death_scores data frame.
#' @param time Survival time.
#' @param status Event status.
#' @param method "cox_lasso" or "rsf".
#' @param train_ratio Training data proportion.
#' @param seed Random seed.
#'
#' @return A death_surv_model object.
#'
#' @export
#'
death_surv_model <- function(scores,
                              time,
                              status,
                              method = c("cox_lasso", "rsf"),
                              train_ratio = 0.7,
                              seed = 123) {
 
  method <- match.arg(method)
  set.seed(seed)
 
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required")
  }
 
  scores <- as.data.frame(scores)
 
  # Remove NA
  complete_idx <- complete.cases(scores) & !is.na(time) & !is.na(status)
  scores <- scores[complete_idx, ]
  time <- time[complete_idx]
  status <- status[complete_idx]
 
  # Split
  n <- nrow(scores)
  train_idx <- sample(1:n, size = round(n * train_ratio))
  test_idx <- setdiff(1:n, train_idx)
 
  if (method == "cox_lasso") {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package 'glmnet' is required")
    }
   
    X_train <- as.matrix(scores[train_idx, ])
    X_test <- as.matrix(scores[test_idx, ])
   
    y_train <- survival::Surv(time[train_idx], status[train_idx])
    y_test <- survival::Surv(time[test_idx], status[test_idx])
   
    # CV for lambda
    cv_fit <- glmnet::cv.glmnet(X_train, y_train, family = "cox", alpha = 1)
    model <- cv_fit
   
    # Predict risk scores
    risk_train <- predict(model, X_train, s = "lambda.min", type = "link")
    risk_test <- predict(model, X_test, s = "lambda.min", type = "link")
   
    # C-index
    c_index <- survival::concordance(y_test ~ risk_test)$concordance
   
    # Coefficients
    coef_mat <- as.matrix(coef(model, s = "lambda.min"))
    importance <- data.frame(
      feature = rownames(coef_mat),
      coefficient = coef_mat[, 1],
      importance = abs(coef_mat[, 1]),
      stringsAsFactors = FALSE
    )
    importance <- importance[importance$importance > 0, ]
    importance <- importance[order(-importance$importance), ]
   
  } else if (method == "rsf") {
    if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
      stop("Package 'randomForestSRC' is required. Install with: install.packages('randomForestSRC')")
    }
   
    train_data <- cbind(time = time[train_idx], status = status[train_idx],
                        scores[train_idx, ])
    test_data <- cbind(time = time[test_idx], status = status[test_idx],
                       scores[test_idx, ])
   
    model <- randomForestSRC::rfsrc(
      Surv(time, status) ~ .,
      data = train_data,
      importance = TRUE
    )
   
    # Predict
    pred <- predict(model, test_data)
    risk_test <- pred$predicted
   
    # C-index
    c_index <- 1 - model$err.rate[length(model$err.rate)]
   
    # Importance
    vimp <- model$importance
    importance <- data.frame(
      feature = names(vimp),
      importance = vimp,
      stringsAsFactors = FALSE
    )
    importance <- importance[order(-importance$importance), ]
  }
 
  result <- list(
    model = model,
    method = method,
    c_index = c_index,
    importance = importance,
    n_train = length(train_idx),
    n_test = length(test_idx),
    feature_names = colnames(scores)
  )
 
  class(result) <- c("death_surv_model", "list")
 
  return(result)
}


#' Print Method for death_surv_model
#' @param x A death_surv_model object
#' @param ... Additional arguments
#' @export
print.death_surv_model <- function(x, ...) {
  cat("\n")
  cat("========================================\n")
  cat("  Survival Prediction Model\n")
  cat("========================================\n\n")
 
  cat("Method:", toupper(x$method), "\n")
  cat("Features:", length(x$feature_names), "\n")
  cat("Training samples:", x$n_train, "\n")
  cat("Test samples:", x$n_test, "\n\n")
 
  cat("Performance:\n")
  cat("  C-index:", round(x$c_index, 4), "\n\n")
 
  cat("Top Features:\n")
  print(head(x$importance, 5), row.names = FALSE)
  cat("\n")
 
  invisible(x)
}
