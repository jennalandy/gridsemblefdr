#' Ensemble
#'
#' @param test_statistics vector of test statistics
#' @param top_grid data frame of method, row number combinations to ensemble
#' @param locfdr_grid data frame where each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data frame where each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data frame where each row is a set of hyperparameters for qvalue
#' @param df degrees of freedom of test statistics, if known
#' @param verbose
#'
#' @return
#' \itemize{
#'   \item fdr - estimated local false discovery rates
#'   \item pi0 - estimated proportion of tests that are null
#' }
#' @export
ensemble <- function(
  test_statistics,
  top_grid,
  locfdr_grid,
  fdrtool_grid,
  qvalue_grid,
  df = NULL,
  verbose = TRUE
) {

  # initialize fdr vector, pi0 estimates, and count (to normalize after)
  topn_fdr <- rep(0, length(test_statistics))
  pi0_est <- 0
  total <- 0

  # iterate through models to ensemble over
  for (i in 1:nrow(top_grid)) {
    method = top_grid[i,'method']
    if (method == 'locfdr') {
      i_fdr <- run_locfdr_row(
        test_statistics = test_statistics,
        locfdr_grid = locfdr_grid,
        row = as.numeric(top_grid[i, 'row'])
      )
    } else if (method == 'fdrtool') {
      i_fdr <- run_fdrtool_row(
        test_statistics = test_statistics,
        fdrtool_grid = fdrtool_grid,
        row = as.numeric(top_grid[i, 'row'])
      )
    } else if (method == 'qvalue') {
      i_fdr <- run_qvalue_row(
        test_statistics = test_statistics,
        qvalue_grid = qvalue_grid,
        row = as.numeric(top_grid[i, 'row']),
        df = df
      )
    }

    if (!is.null(i_fdr)) {
      if ((sum(i_fdr$fdr != 1) > 0)&(sum(i_fdr$fdr != 0) > 0)) {
        topn_fdr <- topn_fdr + i_fdr$fdr
        pi0_est <- pi0_est + i_fdr$pi0
        total <- total + 1
      }
    }

  }

  topn_fdr <- topn_fdr/total
  pi0_est <- pi0_est/total

  return(list(
    fdr = topn_fdr,
    pi0 = pi0_est
  ))
}
