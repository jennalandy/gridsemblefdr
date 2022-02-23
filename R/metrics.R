#' Compute all metrics
#'
#' @param fdr estimated fdr values
#' @param Fdr estimated (tail end) Fdr values
#' @param truth truth values, vector of 0 and 1s
#' @param true_Fdr true Fdr values
#' @param topq vector of booleans for whether each test statistic
#' is in the largest absolute quantile
#'
#' @return labeled list of metrics
#' @export
metrics <- function(fdr, Fdr, truth, true_Fdr, topq) {
  list(
    'pr' = get_prauc(fdr, truth),
    'roc' = get_roc(fdr, truth),
    'brier' = get_brier(fdr, truth),
    'Fdrerror' = get_Fdr_error(Fdr, true_Fdr),
    'pr_topq' = get_prauc(fdr[topq], truth[topq]),
    'roc_topq' = get_roc(fdr[topq], truth[topq]),
    'brier_toq' = get_brier(fdr[topq], truth[topq]),
    'Fdrerror_topq' = get_Fdr_error(Fdr[topq], true_Fdr[topq])
  )
}

#' Calculate PR AUC
#'
#' @param fdr estimated fdr values
#' @param truth truth values, vector of 0 and 1s
#'
#' @return PR AUC
#'
#' @importFrom PRROC pr.curve
#' @export
get_prauc <- function(fdr, truth) {
  pr <- pr.curve(
      scores.class0 = fdr[!as.logical(truth)],
      scores.class1 = fdr[as.logical(truth)]
  )
  return(as.numeric(pr$auc.integral))
}

#' Calculate ROC AUC
#'
#' @param fdr estimated fdr values
#' @param truth truth values, vector of 0 and 1s
#'
#' @return ROC AUC
#'
#' @importFrom pROC roc
#' @export
get_roc <- function(fdr, truth) {

  # direction = ">" accounts for inverse
  # relationship between fdr and y
  r <- roc(
    truth ~ fdr,
    direction = ">",
    levels = levels(as.factor(truth))
  )
  return(as.numeric(r$auc))
}

#' Calculate Brier Score
#'
#' @param fdr estimated fdr values
#' @param truth truth values, vector of 0 and 1s
#'
#' @return Brier Score
#' @export
get_brier <- function(fdr, truth) {
  prob_1 = 1-fdr
  mean((prob_1 - as.numeric(truth))**2)
}

#' Calculate Fdr error
#'
#' @param Fdr estimated (tail end) Fdr values
#' @param true_Fdr true Fdr values
#'
#' @return Fdr error
#' @export
get_Fdr_error <- function(Fdr, true_Fdr) {
  mean((true_Fdr - Fdr)^2)
}
