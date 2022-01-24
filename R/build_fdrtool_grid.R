#' @title Reduce fdrtool grid
#' @description Reduces grid to parameter combinations
#' that can run fdrtool on the provided data without error
#'
#' @param t vector of test statistics
#' @param fdrtool_grid data frame where each
#' @param verbose
#'
#' @return dataframe where each row is a possible set of hyperparameters for fdrtool
#'
#' @examples
reduce_fdrtool_grid <- function(t, fdrtool_grid, verbose = FALSE) {

  ok_rows <- c()
  for (i in 1:nrow(fdrtool_grid)) {
    # for each row, attempt to run fdrtool
    run_i <- run_fdrtool_row(
      t = t,
      fdrtool_grid = fdrtool_grid,
      row = i
    )

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
      length(ok_rows),'/',nrow(fdrtool_grid),
      'fdrtool parameter sets are usable'
    ))
  }

  return (fdrtool_grid[ok_rows,])
}

#' @title Build fdrtool Grid
#' @description Build a grid of possible hyperparameters for fdrtool
#' from separate vectors of hyperparameter options. Final grid only considers
#' hyperparameters that can be run on provided data without error.
#'
#' @param t vector of test statistics
#' @param cutoff.method vector of options for cutoff method hyperparameter
#' @param pct0 vector of options for pct0 hyperparameter
#'
#' @return dataframe where each row is a possible set of hyperparameters for
#' fdrtool on a specific set of test statistics
#' @export
#'
#' @examples
build_fdrtool_grid <- function(
  t,
  cutoff.method = c('fndr','pct0','locfdr'),
  pct0 = 1/c(2:10),
  verbose = FALSE
) {

  fdrtool_grid <- expand.grid(
    cutoff.method = cutoff.method,
    pct0 = pct0
  )

  #pct0 is only used when cutoff.method = 'pct0'
  default_pct0 <- unique(pct0)[1]
  fdrtool_grid <- fdrtool_grid[
    fdrtool_grid$cutoff.method == 'pct0' |
      fdrtool_grid$pct0 == default_pct0,
  ]

  fdrtool_grid_reduced <- reduce_fdrtool_grid(
    t = t,
    fdrtool_grid = fdrtool_grid,
    verbose = verbose
  )

  return(fdrtool_grid_reduced)
}
