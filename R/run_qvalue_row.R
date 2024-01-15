#' @title Run qvalue
#' @description Run qvalue with a specific set of parameters
#'
#' @param test_statistics vector, test statistics
#' @param to_pval_function function, converts test statistics vector to a
#' p-value vector.
#' @param qvalue_grid data.frame, each row is a set of hyperparameters
#' @param row integer, row of qvalue_grid
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
#'
#' @importFrom qvalue qvalue
#' @noRd
run_qvalue_row <- function(
  test_statistics,
  to_pval_function,
  qvalue_grid,
  row,
  returnFdr = TRUE,
  verbose = FALSE
) {
  p <- to_pval_function(
    test_statistics = test_statistics
  )

  res <- NULL

  # try running qvalue and fitting pi0
  tryCatch(
    {
      res <- qvalue::qvalue(
        p = p,
        transf = as.character(qvalue_grid$transf[row]),
        adj = qvalue_grid$adj[row],
        pi0.method = as.character(qvalue_grid$pi0.method[row]),
        smooth.log.pi0 = as.logical(qvalue_grid$smooth.log.pi0[row])
      )
    },
    warning = function(cond) {
    },
    error = function(cond) {
      # don't have enough p-values to warrant estimating pi0
      # try again with lambda = 0
    }
  )

  if (is.null(res)) {
    # if null, run without fitting pi0 (lambda = 0)
    res <- qvalue::qvalue(
      p = p,
      transf = as.character(qvalue_grid$transf[row]),
      adj = qvalue_grid$adj[row],
      pi0.method = as.character(qvalue_grid$pi0.method[row]),
      smooth.log.pi0 = as.logical(qvalue_grid$smooth.log.pi0[row]),
      lambda = 0
    )
  }

  if (is.null(res)) {
    # if still null
    return(NULL)
  } else {
    # otherwise return results
    if (returnFdr) {
      Fdr = Fdr_from_fdr(
        fdr = res$lfdr,
        test_statistics = test_statistics
      )
    } else {
      Fdr = 'not computed'
    }
    return(list(
      'fdr' = res$lfdr,
      'Fdr' = Fdr,
      'pi0' = res$pi0
    ))
  }
}
