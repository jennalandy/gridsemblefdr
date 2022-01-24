#' @title Reduce qvalue grid
#' @description Reduces grid to parameter combinations
#' that can run qvalue on the provided data without error
#'
#' @param t vector of test statistics
#' @param qvalue_grid data frame where each
#' @param verbose
#'
#' @return dataframe where each row is a possible set of hyperparameters for qvalue
#'
#' @examples
reduce_qvalue_grid <- function(t, qvalue_grid, verbose = FALSE) {

  ok_rows <- c()
  for (i in 1:nrow(qvalue_grid)) {
    # for each row, attempt to run qvalue
    run_i <- run_qvalue_row(
      t = t,
      qvalue_grid = qvalue_grid,
      row = i
    )

    # if null, try again with lambda = 0
    if (is.null(run_i)) {
      run_i <- run_qvalue_row(
        t = t,
        qvalue_grid = qvalue_grid,
        row = i,
        lambda0 = TRUE
      )
    }

    # keep this row in final grid if it ran without error,
    # and if it does not estimate fdrs as all 0 or all 1
    if (!is.null(run_i)) {
      if ((sum(run_i$fdr != 1) > 0) & (sum(run_i$fdr != 0) > 0)) {
        ok_rows <- c(ok_rows, i)
      }
    }
  }

  if(verbose) {
    print(paste(
      length(ok_rows),'/',nrow(qvalue_grid),
      'qvalue parameter sets are usable'
    ))
  }

  return (qvalue_grid[ok_rows,])
}

#' @title Build qvalue Grid
#' @description Build a grid of possible hyperparameters for qvalue
#' from separate vectors of hyperparameter options. Final grid only considers
#' hyperparameters that can be run on provided data without error.
#'
#' @param t vector of test statistics
#' @param transf vector of options for transf hyperparameter
#' @param adj vector of options for adj hyperparameter
#' @param pi0.method vector of options for pi0.method hyperparameter
#' @param smooth.log.pi0 vector of options for smooth.log.pi0 hyperparameter
#' @param verbose
#'
#' @return dataframe where each row is a possible set of hyperparameters for
#' qvalue on a specific set of test statistics
#' @export
#'
#' @examples
build_qvalue_grid <- function(
  t,
  transf = c('probit', 'logit'),
  adj = c(0.8,0.9,1,1.1,1.2), # default 1
  pi0.method = c('bootstrap','smoother'),
  smooth.log.pi0 = c(TRUE, FALSE),
  verbose = FALSE
) {

  qvalue_grid <- expand.grid(
    transf = transf,
    adj = adj,
    pi0.method = pi0.method,
    smooth.log.pi0 = smooth.log.pi0
  )

  # smooth.log.pi0 is only used if pi0.method = 'smoother'
  default_smooth.log.pi0 <- unique(smooth.log.pi0)[1]
  qvalue_grid <- qvalue_grid[
    qvalue_grid$pi0.method == 'smoother' |
      qvalue_grid$smooth.log.pi0 == default_smooth.log.pi0,
  ]

  qvalue_grid_reduced <- reduce_qvalue_grid(
    t = t,
    qvalue_grid = qvalue_grid,
    verbose = verbose
  )

  return(qvalue_grid_reduced)
}
