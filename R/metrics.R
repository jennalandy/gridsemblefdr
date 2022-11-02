#' Compute all metrics
#'
#' @param fdr vector, estimated fdr values
#' @param true_fdr vector, true fdr values
#' @param topq vector, booleans for whether each test statistic
#' is in the largest absolute quantile
#'
#' @return labeled list of metrics
metrics <- function(
  fdr, true_fdr, topq
) {
  if (length(fdr) > 0) {
    list(
      'fdrerror' = get_MSE(estimate = fdr, true = true_fdr),
      'fdrerror_topq' = get_MSE(estimate = fdr[topq], true = true_fdr[topq])
    )
  } else {
    list(
      'fdrerror' = NA,
      'fdrerror_topq' = NA
    )
  }
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
