#' Ensemble
#'
#' @param test_statistics vector, test statistics
#' @param top_grid data.frame, specifies packages and row numbers within that
#' grid of the best methods used in the ensemble
#' @param locfdr_grid data.frame, rows are possible hyperparameters for locfdr
#' @param fdrtool_grid data.frame, rows are possible hyperparameters for fdrtool
#' @param qvalue_grid data.frame, rows are possible hyperparameters for qvalue
#' @param df integer, degrees of freedom of test statistics, if known
#'
#' @param focus_metric string, one of one of c('fdrerror'),
#' which metric to optimize in the grid search
#' @param large_abs_metric boolean, if TRUE, only consider focus_metric looking
#' at the large absolute value test statistics (top quartile of abs(t))
#'
#' @param parallel_param BiocParallel object
#' @param verbose boolean
#'
#' @return
#' \itemize{
#'   \item fdr - estimated local false discovery rates
#'   \item pi0 - estimated proportion of tests that are null
#' }
#'
#' @importFrom stats var
#' @noRd
ensemble <- function(
  test_statistics, top_grid, locfdr_grid, fdrtool_grid, qvalue_grid, df = NULL,
  focus_metric = 'fdrerror', large_abs_metric = FALSE, parallel_param = NULL,
  verbose = TRUE
) {
  if (verbose) { message('Ensembling') }
  if (large_abs_metric) { focus_metric = paste0(focus_metric, '_topq') }

  if (ncol(top_grid) > 2) { # small fdrerror is good -> large weight
    weights = (1 - top_grid[,focus_metric])/sum((1 - top_grid[,focus_metric]))
  } else { # no metrics column, equal weighting
    weights = rep(1/nrow(top_grid), nrow(top_grid))
  }

  pi0_and_fdr_table <- do.call(rbind, parlapply(
    X = seq_len(nrow(top_grid)), parallel_param = NULL,
    FUN = function(i) {
      i_fdr <- run_row(
        test_statistics = test_statistics, df = df,
        grids = list(
          'locfdr' = locfdr_grid,
          'fdrtool' = fdrtool_grid,
          'qvalue' = qvalue_grid
        ), method = as.character(top_grid[i,'method']),
        row = as.numeric(top_grid[i, 'row']), returnFdr = FALSE,
        verbose = verbose
      )
      return(c(unname(i_fdr$pi0), i_fdr$fdr))
    }
  ))

  pi0_and_fdr_table = data.frame(pi0_and_fdr_table)
  pi0_and_fdr_means = vapply(
    pi0_and_fdr_table, function(col) {sum(col*weights)}, numeric(1)
  )
  pi0_and_fdr_vars = vapply(
    pi0_and_fdr_table, stats::var, numeric(1)
  )

  fdr <- pi0_and_fdr_means[2:length(pi0_and_fdr_means)]
  fdr_var <- pi0_and_fdr_vars[2:length(pi0_and_fdr_vars)]
  pi0 <- unname(pi0_and_fdr_means[1])
  pi0_var <- unname(pi0_and_fdr_vars[1])

  return(list(
    fdr = fdr, fdr_var = fdr_var, pi0 = pi0, pi0_var = pi0_var
  ))
}
