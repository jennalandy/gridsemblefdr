#' Compute all metrics
#'
#' @param fdr vector, estimated fdr values
#' @param true_fdr vector, true fdr values
#'
#' @return labeled list of metrics
#' @noRd
metrics <- function(
  fdr, true_fdr,
  Fdr, true_Fdr
) {
  if (length(fdr) > 0) {
    list(
      'fdrerror' = get_MSE(estimate = fdr, true = true_fdr),
      'Fdrerror' = get_MSE(estimate = Fdr, true = true_Fdr)
    )
  } else {
    list(
      'fdrerror' = NA,
      'Fdrerror' = NA
    )
  }
}

#' Calculate MSE
#'
#' @param estimate vector, estimated values
#' @param true vector, true values
#'
#' @return double, MSE
#' @noRd
get_MSE <- function(estimate, true) {
  mean((true - estimate)^2)
}
