#' @title Reduce fdrtool grid
#' @description Reduces grid to parameter combinations
#' that can run fdrtool on the provided data without error
#'
#' @param test_statistics vector of test statistics
#' @param fdrtool_grid data frame where each row is a possible set of hyperparameters for fdrtool
#' @param parallel if TRUE, processes are run in parallel
#' @param verbose if TRUE, status updates will be displayed
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#'
#' @return dataframe where each row is a possible set of hyperparameters for fdrtool
#' @export
reduce_fdrtool_grid <- function(
  test_statistics,
  fdrtool_grid,
  parallel = TRUE,
  verbose = FALSE
) {

  n.cores = ifelse(
    parallel,
    min(nrow(fdrtool_grid), parallel::detectCores() - 1),
    1
  )
  cl = parallel::makeCluster(
    n.cores,
    type = 'PSOCK'
  )
  doParallel::registerDoParallel(cl)

  ok_rows <- foreach::foreach(
    i=1:nrow(fdrtool_grid),
    .combine = c,
    .export = c(
      'run_fdrtool_row'
    )
  ) %dopar% {
    # for each row, attempt to run fdrtool
    run_i <- run_fdrtool_row(
      test_statistics = test_statistics,
      fdrtool_grid = fdrtool_grid,
      row = i
    )

    # keep this row in final grid if it ran without error,
    # and if it does not estimate fdrs as all 0 or all 1
    if (!is.null(run_i)) {
      if (
        (sum(run_i$fdr != 1) > 0) &
        (sum(run_i$fdr != 0) > 0) &
        run_i$pi0 <= 1
      ) {
        return(i)
      }
    }
  }
  parallel::stopCluster(cl)

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
#' @param test_statistics vector of test statistics
#' @param cutoff.method vector of options for cutoff method hyperparameter.
#' `cutoff.method` is one of "fndr" (default), "pct0", "locfdr".
#' @param pct0 vector of options for pct0 hyperparameter. `pct0` is the
#' fraction of data used for fitting null model - only if cutoff.method="pct0".
#' @param parallel if TRUE, processes are run in parallel
#' @param verbose if TRUE, status updates will be displayed
#'
#' @return dataframe where each row is a possible set of hyperparameters for
#' fdrtool on a specific set of test statistics
#' @export
build_fdrtool_grid <- function(
  test_statistics,
  cutoff.method = c('fndr','pct0','locfdr'),
  pct0 = 1/c(2:10),
  parallel = TRUE,
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
    test_statistics = test_statistics,
    fdrtool_grid = fdrtool_grid,
    parallel = parallel,
    verbose = verbose
  )

  return(fdrtool_grid_reduced)
}
