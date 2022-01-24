#' @title P from T
#' @description get p-value from test statistic
#'
#' @param t vector of t statistics
#' @param df degrees of freedom to compute p-value from test statistics,
#' assume standard normal if NULL
#'
#' @return
#'
#' @examples
p_from_t <- function(t, df = NULL, sides = 'one') {

  if (is.null(df)) {
    # assume standard normal
    one_sided <- lapply(t, function(z) {
      pnorm(-1*abs(z))
    }) %>%
      unlist()
  } else {
    # assume t_df
    one_sided <-lapply(t, function(z) {
      pt(-1*abs(z), df = df)
    }) %>%
      unlist()
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
#' @param t vector of test statistics
#' @param qvalue_grid data frame where each row is a set of hyperparameters
#' @param row row of qvalue_grid, i.e. which set of hyperparameters to run qvalue with
#' @param lambda0 if TRUE, pi0 is not estimated (necessary to avoid some errors)
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
#' @export
#'
#' @examples
run_qvalue_row <- function(
  t, qvalue_grid, row,
  lambda0 = FALSE, df = NULL,
  verbose = FALSE
) {

  p <- p_from_t(
    t = t,
    df = df,
    sides = 'one'
  )

  tryCatch(
    {
      if (lambda0) {
        res <- qvalue(
          p = p,
          transf = as.character(qvalue_grid$transf[row]),
          adj = qvalue_grid$adj[row],
          pi0.method = as.character(qvalue_grid$pi0.method[row]),
          smooth.log.pi0 = as.logical(qvalue_grid$smooth.log.pi0[row]),
          lambda = 0
        )
      } else {
        res <- qvalue(
          p = p,
          transf = qvalue_grid$transf[row],
          adj = qvalue_grid$adj[row],
          pi0.method = as.character(qvalue_grid$pi0.method[row]),
          smooth.log.pi0 = as.logical(qvalue_grid$smooth.log.pi0[row])
        )
      }

      return(list(
        'fdr' = res$lfdr,
        'Fdr' = Fdr_from_fdr(res$lfdr, t),
        'pi0' = res$pi0
      ))
    },
    warning = function(cond) {
      print(cond)
    },
    error = function(cond) {
      # don't have enough p-values to warrant estimating pi0
      # try again with lambda = 0

      if (verbose) {
        print(paste("Error.", cond, "\nTry with lambda0 = TRUE."))
      }

    }
  )
}
