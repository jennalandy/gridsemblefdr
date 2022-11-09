#' @title Run fdrtool with a specific set of parameters
#'
#' @param test_statistics vector, test statistics
#' @param fdrtool_grid data.frame, each row is a set of hyperparameters
#' @param row integer, row of fdrtool_grid
#' @param returnFdr boolean, whether to calculate Fdr form fdr
#'
#' @return
#' list of estimates
#' \itemize{
#'   \item `fdr`: vector, local false discovery rates
#'   \item `Fdr`: vector, tail-end false discovery rates if returnFdr = TRUE
#'   \item `pi0`: double, proportion of tests that are null
#' }
#'
#' @importFrom fdrtool fdrtool
#' @noRd
run_fdrtool_row <- function(
  test_statistics,
  fdrtool_grid,
  row,
  returnFdr = TRUE
) {
  tryCatch(
    {
      res <- fdrtool::fdrtool(
        x = test_statistics,
        plot = 0,
        verbose = FALSE,
        cutoff.method = as.character(
          fdrtool_grid$cutoff.method[row]
        ),
        pct0 = fdrtool_grid$pct0[row]
      )

      if (returnFdr) {
        Fdr = Fdr_from_fdr(
          fdr = res$lfdr,
          test_statistics = test_statistics,
          direction = 'left'
        )
      } else {
        Fdr = 'not computed'
      }

      return(list(
        'fdr' = res$lfdr,
        'Fdr' = Fdr,
        'pi0' = unname(unlist(res$param[1, 'eta0']))
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
