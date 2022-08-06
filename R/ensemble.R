#' Ensemble
#'
#' @param test_statistics vector, test statistics
#' @param top_grid data.frame, specifies packages and row numbers within that grid of
#' the best methods used in the ensemble
#' @param locfdr_grid data.frame, each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data.frame, each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data.frame, each row is a set of hyperparameters for qvalue
#' @param df integer, degrees of freedom of test statistics, if known
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#' @param verbose boolean
#'
#' @return
#' \itemize{
#'   \item fdr - estimated local false discovery rates
#'   \item pi0 - estimated proportion of tests that are null
#' }
ensemble <- function(
  test_statistics,
  top_grid,
  focus_metric,
  large_abs_metric,
  locfdr_grid,
  fdrtool_grid,
  qvalue_grid,
  df = NULL,
  parallel_param = NULL,
  verbose = TRUE
) {

  if (large_abs_metric & focus_metric %in% c('brier','Fdrerror','fdrerror','pr','roc')) {
    focus_metric = paste(focus_metric, '_topq', sep = '')
  }
  if (ncol(top_grid) > 2) {
    # means there was simulation and grid search
    # there are metrics columns
    if (any(startsWith(
      # smaller is better for these three
      focus_metric, c('Fdrerror','brier','fdrerror')
    ))) {
      # small metric -> large weight
      weights = (1 - top_grid[,focus_metric])/sum((1 - top_grid[,focus_metric]))
    } else {
      # large metric -> large weight
      weights = top_grid[,focus_metric]/sum(top_grid[,focus_metric])
    }
  } else {
    # no metrics column, equal weighting
    weights = rep(1/nrow(top_grid), nrow(top_grid))
  }

  # print(top_grid)

  # nrow(top_grid) rows x len(test_statistics) + 1 columns
  pi0_and_fdr_table <- do.call(rbind, parlapply(
    X = 1:nrow(top_grid),
    parallel_param = NULL,
    FUN = function(i) {
      method = top_grid[i,'method']
      if (is.na(method)) {
        print('null method??')
        print(top_grid[i,])
      }
      if (method == 'locfdr') {
        i_fdr <- run_locfdr_row(
          test_statistics = test_statistics,
          locfdr_grid = locfdr_grid,
          row = as.numeric(top_grid[i, 'row']),
          returnFdr = FALSE
        )
      } else if (method == 'fdrtool') {
        i_fdr <- run_fdrtool_row(
          test_statistics = test_statistics,
          fdrtool_grid = fdrtool_grid,
          row = as.numeric(top_grid[i, 'row']),
          returnFdr = FALSE
        )
      } else if (method == 'qvalue') {
        i_fdr <- run_qvalue_row(
          test_statistics = test_statistics,
          qvalue_grid = qvalue_grid,
          row = as.numeric(top_grid[i, 'row']),
          df = df,
          returnFdr = FALSE
        )
      }

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
    var
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
