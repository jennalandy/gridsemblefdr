#' Compute Classification Metrics
#'
#' @param fdr vector, local fdr values
#' @param test_statistics vector, test statistics
#' @param truth vector, truth values
#' @param cutoff double, fdr <  cutoff predicted as not-null
#'
#' @return named list of metrics based on specified cutoff
classification_metrics <- function(fdr, test_statistics, truth, cutoff) {
  pred = rep(1, length(test_statistics))
  pred[fdr > cutoff] = 0

  TP = sum(pred == 1 & truth == 1)
  FP = sum(pred == 1 & truth == 0)
  TN = sum(pred == 0 & truth == 0)
  FN = sum(pred == 0 & truth == 1)

  accuracy = (TP + TN)/(TP + TN + FP + FN)
  precision = TP/(TP+FP)
  recall = TP/(TP+FN) # fnr = 1-recall, same as sensitivity
  specificity = TN/(TN + FP) # fpr = 1-specificity
  f1 = (2 * precision * recall) / (precision+recall)

  return(list(
    'pred_pos' = sum(pred == 1),
    'cutoff' = cutoff,
    'accuracy' = accuracy,
    'precision' = precision,
    'recall' = recall,
    'specificity' = specificity,
    'f1' = f1
  ))
}

#' Compute all metrics
#'
#' @param fdr vector, estimated fdr values
#' @param Fdr vector, estimated (tail end) Fdr values
#' @param truth vector, truth values, vector of 0 and 1s
#' @param true_Fdr vector, true Fdr values
#' @param true_fdr vector, true fdr values
#' @param topq vector, booleans for whether each test statistic
#' is in the largest absolute quantile
#'
#' @return labeled list of metrics
metrics <- function(
  fdr, Fdr,
  test_statistics, truth,
  true_Fdr, true_fdr,
  topq,
  cutoff = 0.2
) {
  classification_metrics <- classification_metrics(
    fdr = fdr,
    test_statistics = test_statistics,
    truth = truth,
    cutoff = cutoff
  )
  list(
    'pr' = get_prauc(fdr, truth),
    'roc' = get_roc(fdr, truth),
    'brier' = get_brier(fdr, truth),
    'Fdrerror' = get_MSE(estimate = Fdr, true = true_Fdr),
    'fdrerror' = get_MSE(estimate = fdr, true = true_fdr),
    'pr_topq' = get_prauc(fdr[topq], truth[topq]),
    'roc_topq' = get_roc(fdr[topq], truth[topq]),
    'brier_toq' = get_brier(fdr[topq], truth[topq]),
    'Fdrerror_topq' = get_MSE(estimate = Fdr[topq], true = true_Fdr[topq]),
    'fdrerror_topq' = get_MSE(estimate = fdr[topq], true = true_fdr[topq]),
    'accuracy' = classification_metrics$accuracy,
    'precision' = classification_metrics$precision,
    'recall' = classification_metrics$recall,
    'specificity' = classification_metrics$specificity,
    'f1' = classification_metrics$f1
  )
}

#' Calculate PR AUC
#'
#' @param fdr vector, estimated fdr values
#' @param truth vector, truth values, vector of 0 and 1s
#'
#' @return double, PR AUC
#'
#' @importFrom PRROC pr.curve
get_prauc <- function(fdr, truth) {
  tryCatch({
    pr <- PRROC::pr.curve(
        scores.class0 = fdr[!as.logical(truth)],
        scores.class1 = fdr[as.logical(truth)]
    )
    return(as.numeric(pr$auc.integral))
  }, error = function(e) {
    return(NA)
  })
}

#' Calculate ROC AUC
#'
#' @param fdr vector, estimated fdr values
#' @param truth vector, truth values
#'
#' @return double, ROC AUC
#'
#' @importFrom pROC roc
get_roc <- function(fdr, truth) {

  if (sum(truth) == 0 | sum(!truth) == 0) {
    return(0)
  }

  # direction = ">" accounts for inverse
  # relationship between fdr and y
  tryCatch({
    r <- pROC::roc(
      truth ~ fdr,
      direction = ">",
      levels = levels(as.factor(truth))
    )
    return(as.numeric(r$auc))
  }, error = function(e) {
    return(NA)
  })
}

#' Calculate Brier Score
#'
#' @param fdr vector, estimated fdr values
#' @param truth vector, truth values
#'
#' @return double, Brier Score
get_brier <- function(fdr, truth) {
  prob_1 = 1-fdr
  mean((prob_1 - as.numeric(truth))**2)
}

#' Calculate MSE
#'
#' @param estimate vector, estimated values
#' @param true vector, true values
#'
#' @return double, MSE
get_MSE <- function(estimate, true) {
  mean((true - estimate)^2)
}
