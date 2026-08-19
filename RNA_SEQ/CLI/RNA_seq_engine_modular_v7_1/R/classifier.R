# Extracted from RNA_SEQ_CLI_V7_1_STAT_CLEARANCE_CANDIDATE_2026_08_10
# Modularization only: statistical behavior intentionally unchanged.

classifier_logcpm_sample_local <- function(counts_matrix, library_sizes = NULL, prior_count = 0.5) {
  # Leakage-safe classifier preprocessing: each sample is transformed using only
  # its own counts and full raw library size. No cross-sample TMM factor, centering,
  # scaling, imputation, or other statistic is learned from held-out samples.
  counts_matrix <- as.matrix(counts_matrix)
  storage.mode(counts_matrix) <- "double"
  if (is.null(library_sizes)) library_sizes <- colSums(counts_matrix)
  library_sizes <- as.numeric(library_sizes)
  if (length(library_sizes) != ncol(counts_matrix)) stop("Classifier library-size vector does not match sample count.")
  if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) stop("Classifier sample has non-positive or invalid library size.")
  log2(sweep(counts_matrix + prior_count, 2, library_sizes, "/") * 1e6)
}

select_train_features_by_de <- function(train_counts, train_pheno, phenotype_column, reference_group, contrast_arg, covariates, top_n) {
  resolved <- resolve_groups(train_pheno, phenotype_column, reference_group, contrast_arg)
  group <- resolved$group
  dge <- DGEList(counts = train_counts, group = group)
  design <- make_design(train_pheno, group, covariates)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  if (nrow(dge) < 2) stop("Too few train-set genes after filterByExpr for RF feature selection.")
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot = FALSE)
  contrast_string <- paste0(make.names(resolved$numerator), "-", make.names(resolved$denominator))
  contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design)
  fit <- eBayes(contrasts.fit(lmFit(v, design), contrast_matrix))
  res <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
  res$Ensembl_IDs <- rownames(res)
  sig <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < fdr_cutoff & abs(res$logFC) >= logfc_cutoff, , drop = FALSE]
  if (nrow(sig) < 2) {
    log_msg("RF train-only feature selection found fewer than 2 significant genes; using top ", top_n, " train-set DE-ranked genes as exploratory features.")
    sig <- head(res[order(res$adj.P.Val), , drop = FALSE], top_n)
  }
  strip_ens_version(head(sig$Ensembl_IDs, top_n))
}

run_rf_classifier <- function(counts_raw, pheno, phenotype_column, reference_group, contrast_arg, covariates, top_n, model_choice = "auto") {
  pheno_vec <- droplevels(as.factor(pheno[[phenotype_column]]))
  positive_class <- resolve_groups(pheno, phenotype_column, reference_group, contrast_arg)$numerator
  if (nrow(pheno) < 6) stop("Too few samples for train/test classifier.")
  if (any(table(pheno_vec) < 2)) stop("Each phenotype class needs at least 2 samples for train/test classifier.")

  set.seed(seed)
  train_idx <- caret::createDataPartition(pheno_vec, p = 0.7, list = FALSE)
  train_samples <- rownames(pheno)[train_idx]
  test_samples <- setdiff(rownames(pheno), train_samples)
  if (length(test_samples) < 1) stop("Train/test split produced empty test set.")

  train_counts <- counts_raw[, train_samples, drop = FALSE]
  test_counts <- counts_raw[, test_samples, drop = FALSE]
  train_pheno <- pheno[train_samples, , drop = FALSE]

  feature_ids <- select_train_features_by_de(train_counts, train_pheno, phenotype_column, reference_group, contrast_arg, covariates, top_n)
  train_library_sizes <- colSums(train_counts)
  test_library_sizes <- colSums(test_counts)
  rownames(train_counts) <- strip_ens_version(rownames(train_counts))
  rownames(test_counts) <- strip_ens_version(rownames(test_counts))
  feature_ids <- intersect(feature_ids, rownames(train_counts))
  if (length(feature_ids) < 2) stop("Fewer than two train-selected classifier features remained.")

  x_train <- data.frame(t(classifier_logcpm_sample_local(train_counts[feature_ids, , drop = FALSE], train_library_sizes)), check.names = FALSE)
  x_test <- data.frame(t(classifier_logcpm_sample_local(test_counts[feature_ids, , drop = FALSE], test_library_sizes)), check.names = FALSE)
  y_train <- droplevels(as.factor(pheno[train_samples, phenotype_column]))

  # Model choice/fallback and fitting use training data only. Predict before reading
  # the held-out labels; labels are touched only for final scoring below.
  set.seed(seed)
  pred_obj <- train_classifier_model(x_train, x_test, y_train, model_choice = model_choice, positive_class = positive_class)
  actual <- factor(pheno[test_samples, phenotype_column], levels = levels(y_train))

  preds <- data.frame(
    Fold = 1L,
    Sample = test_samples,
    Actual = as.character(actual),
    Predicted = as.character(pred_obj$pred),
    Prob_Positive = pred_obj$prob_pos,
    Positive_Class = positive_class,
    Model = pred_obj$model,
    Validation = "train_test",
    stringsAsFactors = FALSE
  )
  imp <- pred_obj$importance
  importance_rows <- list()
  if (nrow(imp) > 0) {
    imp$Fold <- 1L
    imp$Model <- pred_obj$model
    importance_rows[[1]] <- imp
  }
  compute_classifier_outputs(preds, importance_rows, levels(y_train), validation_label = "train_test", positive_class = positive_class)
}

train_classifier_model <- function(train_data, test_data, y_train, model_choice = "auto", positive_class = NULL) {
  y_train <- droplevels(as.factor(y_train))
  if (length(levels(y_train)) == 2) {
    if (is.null(positive_class)) positive_class <- levels(y_train)[2]
    if (!positive_class %in% levels(y_train)) stop("Requested classifier positive class is absent from training labels: ", positive_class)
    negative_class <- setdiff(levels(y_train), positive_class)[[1]]
  } else {
    positive_class <- NULL
    negative_class <- NULL
  }
  model_choice <- tolower(model_choice)
  if (model_choice == "auto") model_choice <- "rf"

  if (model_choice == "logistic" && length(levels(y_train)) == 2 && nrow(train_data) > 3) {
    dat <- data.frame(y = as.numeric(y_train == positive_class), train_data, check.names = FALSE)
    fit <- tryCatch(stats::glm(y ~ ., data = dat, family = stats::binomial()), error = function(e) NULL)
    if (!is.null(fit)) {
      prob <- tryCatch(as.numeric(stats::predict(fit, newdata = test_data, type = "response")), error = function(e) NULL)
      if (!is.null(prob) && all(is.finite(prob))) {
        pred <- factor(ifelse(prob >= 0.5, positive_class, negative_class), levels = levels(y_train))
        return(list(pred = pred, prob_pos = prob, positive_class = positive_class, model = "logistic", importance = data.frame(Feature = colnames(train_data), MeanDecreaseAccuracy = NA_real_)))
      }
    }
    log_msg("Logistic classifier failed or was unstable; falling back to random forest.")
  }

  rf_model <- randomForest(x = train_data, y = y_train, ntree = 500, importance = TRUE)
  pred <- predict(rf_model, newdata = test_data, type = "response")
  probs <- tryCatch(predict(rf_model, newdata = test_data, type = "prob"), error = function(e) NULL)
  prob_pos <- if (!is.null(probs) && !is.null(positive_class) && positive_class %in% colnames(probs)) as.numeric(probs[, positive_class]) else rep(NA_real_, nrow(test_data))
  imp <- randomForest::importance(rf_model, type = 1)
  imp_df <- data.frame(Feature = rownames(imp), MeanDecreaseAccuracy = as.numeric(imp[, 1]), stringsAsFactors = FALSE)
  list(pred = pred, prob_pos = prob_pos, positive_class = positive_class, model = "rf", importance = imp_df)
}

compute_classifier_outputs <- function(preds, importance_rows, pheno_levels, validation_label, positive_class = NULL) {
  pred_factor <- factor(preds$Predicted, levels = pheno_levels)
  actual_factor <- factor(preds$Actual, levels = pheno_levels)
  if (length(pheno_levels) == 2) {
    if (is.null(positive_class)) positive_class <- pheno_levels[[2]]
    if (!positive_class %in% pheno_levels) stop("Classifier positive class is not a phenotype level: ", positive_class)
    negative_class <- setdiff(pheno_levels, positive_class)[[1]]
    cm <- caret::confusionMatrix(pred_factor, actual_factor, positive = positive_class)
    sens <- unname(cm$byClass["Sensitivity"])
    spec <- unname(cm$byClass["Specificity"])
    auc_val <- NA_real_
    roc_plot <- plotly::plot_ly() %>% layout(title = "ROC unavailable")
    if (!all(is.na(preds$Prob_Positive))) {
      roc_obj <- pROC::roc(actual_factor, preds$Prob_Positive, levels = c(negative_class, positive_class), direction = "<", quiet = TRUE)
      auc_val <- round(as.numeric(roc_obj$auc), 3)
      roc_df <- data.frame(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities)
      roc_plot <- plot_ly(data = roc_df, x = ~FPR, y = ~TPR, type = "scatter", mode = "lines", line = list(color = "#1f77b4", width = 2)) %>%
        layout(title = paste("Classifier ROC (AUC =", auc_val, ")"), xaxis = list(title = "False Positive Rate"), yaxis = list(title = "True Positive Rate"), showlegend = FALSE)
    }
    metrics <- data.frame(Accuracy = unname(cm$overall["Accuracy"]), Sensitivity = round(sens, 3), Specificity = round(spec, 3), `AUC (ROC)` = ifelse(is.na(auc_val), "N/A", auc_val), Positive_Class = positive_class, Validation = validation_label, Note = "Exploratory classifier; feature selection is inside the training split/fold and classifier preprocessing is sample-local, with metrics computed from held-out predictions. Sensitivity, specificity, probability, and ROC all use Positive_Class consistently.", check.names = FALSE)
  } else {
    cm <- caret::confusionMatrix(pred_factor, actual_factor)
    by_class <- as.data.frame(cm$byClass)
    metrics <- data.frame(Accuracy = unname(cm$overall["Accuracy"]), MacroSensitivity = round(mean(by_class$Sensitivity, na.rm = TRUE), 3), MacroSpecificity = round(mean(by_class$Specificity, na.rm = TRUE), 3), `AUC (ROC)` = "N/A", Validation = validation_label, Note = "Exploratory multi-class classifier.", check.names = FALSE)
    roc_plot <- plotly::plot_ly() %>% layout(title = "ROC only supported for 2-class problems")
  }

  if (length(importance_rows) > 0) {
    imp_all <- do.call(rbind, importance_rows)
    imp_all$Rank <- ave(-imp_all$MeanDecreaseAccuracy, imp_all$Fold, FUN = rank, ties.method = "first")
    stability <- aggregate(list(selection_count = imp_all$Feature), by = list(Feature = imp_all$Feature), FUN = length)
    stability$selection_frequency <- stability$selection_count / length(unique(imp_all$Fold))
    mean_imp <- aggregate(MeanDecreaseAccuracy ~ Feature, imp_all, mean, na.rm = TRUE)
    mean_rank <- aggregate(Rank ~ Feature, imp_all, mean, na.rm = TRUE)
    stability <- merge(stability, mean_imp, by = "Feature", all.x = TRUE)
    stability <- merge(stability, mean_rank, by = "Feature", all.x = TRUE)
    stability <- stability[order(-stability$selection_frequency, stability$Rank), , drop = FALSE]
  } else {
    imp_all <- data.frame()
    stability <- data.frame()
  }

  list(predictions = preds, metrics = metrics, roc_plot = roc_plot, importance = imp_all, stability = stability)
}

run_cv_classifier <- function(counts_raw, pheno, phenotype_column, reference_group, contrast_arg, covariates, top_n, model_choice = "auto") {
  pheno_vec <- droplevels(as.factor(pheno[[phenotype_column]]))
  positive_class <- resolve_groups(pheno, phenotype_column, reference_group, contrast_arg)$numerator
  if (nrow(pheno) < 4) stop("Too few samples for CV classifier.")
  if (any(table(pheno_vec) < 2)) stop("Each phenotype class needs at least 2 samples for CV classifier.")
  pheno_levels <- levels(pheno_vec)
  k <- min(5, as.integer(min(table(pheno_vec))))
  k <- max(2, k)
  set.seed(seed)
  folds <- caret::createFolds(pheno_vec, k = k, returnTrain = FALSE)

  prediction_rows <- list()
  importance_rows <- list()
  for (fold_i in seq_along(folds)) {
    test_samples <- rownames(pheno)[folds[[fold_i]]]
    train_samples <- setdiff(rownames(pheno), test_samples)
    train_counts <- counts_raw[, train_samples, drop = FALSE]
    test_counts <- counts_raw[, test_samples, drop = FALSE]
    train_pheno <- pheno[train_samples, , drop = FALSE]

    feature_ids <- tryCatch(select_train_features_by_de(train_counts, train_pheno, phenotype_column, reference_group, contrast_arg, covariates, top_n), error = function(e) {
      log_msg("Fold ", fold_i, " DE feature selection failed: ", conditionMessage(e))
      character(0)
    })
    train_library_sizes <- colSums(train_counts)
    test_library_sizes <- colSums(test_counts)
    rownames(train_counts) <- strip_ens_version(rownames(train_counts))
    rownames(test_counts) <- strip_ens_version(rownames(test_counts))
    feature_ids <- intersect(feature_ids, rownames(train_counts))
    if (length(feature_ids) < 2) {
      log_msg("Fold ", fold_i, " skipped: fewer than 2 train-selected features.")
      next
    }

    x_train <- data.frame(t(classifier_logcpm_sample_local(train_counts[feature_ids, , drop = FALSE], train_library_sizes)), check.names = FALSE)
    x_test <- data.frame(t(classifier_logcpm_sample_local(test_counts[feature_ids, , drop = FALSE], test_library_sizes)), check.names = FALSE)
    y_train <- droplevels(as.factor(pheno[train_samples, phenotype_column]))

    # Predict before reading held-out labels for scoring. Fold assignment may be
    # stratified by class, but held-out labels do not enter feature selection,
    # preprocessing, fitting, model fallback, or prediction.
    pred_obj <- train_classifier_model(x_train, x_test, y_train, model_choice = model_choice, positive_class = positive_class)
    actual <- factor(pheno[test_samples, phenotype_column], levels = pheno_levels)
    prediction_rows[[length(prediction_rows) + 1]] <- data.frame(Fold = fold_i, Sample = test_samples, Actual = as.character(actual), Predicted = as.character(pred_obj$pred), Prob_Positive = pred_obj$prob_pos, Positive_Class = positive_class, Model = pred_obj$model, Validation = "cv", stringsAsFactors = FALSE)
    imp <- pred_obj$importance
    if (nrow(imp) > 0) {
      imp$Fold <- fold_i
      imp$Model <- pred_obj$model
      importance_rows[[length(importance_rows) + 1]] <- imp
    }
  }

  if (length(prediction_rows) == 0) stop("No CV folds produced classifier predictions.")
  preds <- do.call(rbind, prediction_rows)
  compute_classifier_outputs(preds, importance_rows, pheno_levels, validation_label = paste0("cv_", length(unique(preds$Fold)), "fold"), positive_class = positive_class)
}

