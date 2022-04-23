#' Ensemble
#'
#' @param test_statistics vector of test statistics
#' @param top_grid data frame of method, row number combinations to ensemble
#' @param locfdr_grid data frame where each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data frame where each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data frame where each row is a set of hyperparameters for qvalue
#' @param df degrees of freedom of test statistics, if known
#' @param parallel_param if not NULL, processes are run in parallel
#' @param verbose if TRUE, status updates will be displayed
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
  parallel_param = NULL,
  verbose = TRUE
) {

  # nrow(top_grid) rows x len(test_statistics) + 1 columns
  pi0_and_fdr_table <- do.call(rbind, parlapply(
    X = 1:nrow(top_grid),
    parallel_param = NULL,
    FUN = function(i) {
      method = top_grid[i,'method']
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

  # take the mean of each column
  # len(test_statistics) + 1
  pi0_and_fdr_means = sapply(pi0_and_fdr_table, mean)
  pi0_and_fdr_vars = sapply(pi0_and_fdr_table, var)

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
