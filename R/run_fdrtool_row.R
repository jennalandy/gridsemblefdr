#' @title Run fdrtool
#' @description Run fdrtool with a specific set of parameters
#'
#' @param t vector of test statistics
#' @param fdrtool_grid data frame where each row is a set of hyperparameters
#' @param row row of fdrtool_grid, i.e. which set of hyperparameters to run fdrtool with
#'
#' @return
#'
#' @examples
#'
#' @importFrom fdrtool fdrtool
run_fdrtool_row <- function(t, fdrtool_grid, row) {
  tryCatch(
    {
      res <- fdrtool::fdrtool(
        x = t, plot = 0, verbose = FALSE,
        cutoff.method = as.character(
          fdrtool_grid$cutoff.method[row]
        ),
        pct0 = fdrtool_grid$pct0[row]
      )

      return(list(
        'fdr' = res$lfdr,
        'Fdr' = Fdr_from_fdr(res$lfdr, t),
        'pi0' = res$param[1, 'eta0']
      ))
    },
    warning = function(cond) {
      return(NULL)
    },
    error = function(cond) {
      return(NULL)
    }
  )
}
