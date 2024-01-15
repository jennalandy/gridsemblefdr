#' Run grid row
#'
#' @param test_statistics vector, test statistics
#' @param to_pval_function function, converts test statistics vector to a
#' p-value vector.
#' @param grids list, grids named by method
#' @param method string, one of c('locfdr','fdrtool','qvalue')
#' @param row integer, which row of method grid to consider
#' @param returnFdr boolean, whether to calculate Fdr form fdr
#' @param verbose boolean
#'
#' @return
#' list of estimates
#' \itemize{
#'   \item `fdr`: vector, local false discovery rates
#'   \item `Fdr`: vector, tail-end false discovery rates if returnFdr = TRUE
#'   \item `pi0`: double, proportion of tests that are null
#' }
#' @noRd
run_row <- function(
    test_statistics,
    to_pval_function,
    grids,
    method,
    row,
    returnFdr = TRUE,
    verbose = FALSE
) {
  grid = grids[[method]]
  if (method == 'locfdr') {
    row_res <- run_locfdr_row(
      test_statistics = test_statistics,
      locfdr_grid = grid,
      row = row,
      returnFdr = returnFdr
    )
  } else if (method == 'fdrtool') {
    row_res <- run_fdrtool_row(
      test_statistics = test_statistics,
      fdrtool_grid = grid,
      row = row,
      returnFdr = returnFdr
    )
  } else if (method == 'qvalue') {
    row_res <- run_qvalue_row(
      test_statistics = test_statistics,
      to_pval_function = to_pval_function,
      qvalue_grid = grid,
      row = row,
      returnFdr = returnFdr,
      verbose = verbose
    )
  }
  return(row_res)
}
