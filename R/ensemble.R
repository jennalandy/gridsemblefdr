#' Ensemble
#'
#' @param test_statistics vector, test statistics
#' @param top_grid data.frame, specifies packages and row numbers within that grid of
#' the best methods used in the ensemble
#' @param locfdr_grid data.frame, each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data.frame, each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data.frame, each row is a set of hyperparameters for qvalue
#' @param df integer, degrees of freedom of test statistics, if known
#'
#' @param focus_metric string, one of one of c('fdrerror'),
#' which metric to optimize in the grid search
#' @param large_abs_metric boolean, if TRUE, only consider focus_metric looking at the
#' large absolute value test statistics (top quartile of abs(t))
#'
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#' @param verbose boolean
#'
#' @return
#' \itemize{
#'   \item fdr - estimated local false discovery rates
#'   \item pi0 - estimated proportion of tests that are null
#' }
#'
#' @importFrom stats var
ensemble <- function(
  test_statistics,
  top_grid,
  locfdr_grid,
  fdrtool_grid,
  qvalue_grid,
  df = NULL,
  focus_metric = 'fdrerror',
  large_abs_metric = FALSE,
  parallel_param = NULL,
  verbose = TRUE
) {

  if (large_abs_metric) {
    focus_metric = paste0(focus_metric, '_topq')
  }
  if (ncol(top_grid) > 2) {
    # means there was simulation and grid search -> there are metrics columns
    # small fdrerror is good -> large weight
    weights = (1 - top_grid[,focus_metric])/sum((1 - top_grid[,focus_metric]))
  } else {
    # no metrics column, equal weighting
    weights = rep(1/nrow(top_grid), nrow(top_grid))
  }

  # nrow(top_grid) rows x len(test_statistics) + 1 columns
  pi0_and_fdr_table <- do.call(rbind, parlapply(
    X = 1:nrow(top_grid),
    parallel_param = NULL,
    FUN = function(i) {
      i_fdr <- run_row(
        test_statistics = test_statistics,
        df = df,
        grids = list(
          'locfdr' = locfdr_grid,
          'fdrtool' = fdrtool_grid,
          'qvalue' = qvalue_grid
        ),
        method = as.character(top_grid[i,'method']),
        row = as.numeric(top_grid[i, 'row']),
        returnFdr = FALSE,
        verbose = verbose
      )

      return(c(unname(i_fdr$pi0), i_fdr$fdr))
    }
  ))

  pi0_and_fdr_table = data.frame(pi0_and_fdr_table)

  # take the weighted mean of each column
  # len(test_statistics) + 1
  pi0_and_fdr_means = sapply(
    pi0_and_fdr_table,
    function(col) {sum(col*weights)}
  )
  # variance of values included in weighted means
  # should this be weighted var??
  pi0_and_fdr_vars = sapply(
    pi0_and_fdr_table,
    stats::var
  )

  fdr <- pi0_and_fdr_means[2:length(pi0_and_fdr_means)]
  fdr_var <- pi0_and_fdr_vars[2:length(pi0_and_fdr_vars)]
  pi0 <- pi0_and_fdr_means[1]
  pi0_var <- pi0_and_fdr_vars[1]

  return(list(
    fdr = fdr,
    fdr_var = fdr_var,
    pi0 = pi0,
    pi0_var = pi0_var
  ))
}
