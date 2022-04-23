#' @title P from T
#' @description get p-value from test statistic
#'
#' @param test_statistics vector of t statistics
#' @param df degrees of freedom to compute p-value from test statistics,
#' assume standard normal if NULL
#'
#' @return vector of p-values
#' @export
p_from_t <- function(test_statistics, df = NULL, sides = 'one') {

  if (is.null(df)) {
    # assume standard normal
    one_sided <- unlist(lapply(test_statistics, function(z) {
      pnorm(-1*abs(z))
    }))
  } else {
    # assume t_df
    one_sided <-unlist(lapply(test_statistics, function(z) {
      pt(-1*abs(z), df = df)
    }))
  }

  if (sides == 'two') {
    return (2*one_sided)
  } else {
    return(one_sided)
  }
}

#' @title Run qvalue
#' @description Run fdrtool with a specific set of parameters
#'
#' @param test_statistics vector of test statistics
#' @param qvalue_grid data frame where each row is a set of hyperparameters
#' @param row row of qvalue_grid, i.e. which set of hyperparameters to run qvalue with
#' @param df degrees of freedom to compute p-value from test statistics,
#' assume standard normal if NULL
#' @param verbose
#'
#' @return
#' \itemize{
#'   \item fdr - estimated local false discovery rates
#'   \item Fdr - estimated left tail-end false discovery rates
#'   \item pi0 - estimated proportion of tests that are null
#' }
#'
#' @importFrom qvalue qvalue
#' @export
run_qvalue_row <- function(
  test_statistics,
  qvalue_grid,
  row,
  df = NULL,
  returnFdr = TRUE,
  verbose = FALSE
) {

  p <- p_from_t(
    test_statistics = test_statistics,
    df = df,
    sides = 'two'
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
        test_statistics = test_statistics,
        direction = 'left'
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
