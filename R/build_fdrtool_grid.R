#' @title Check Fdrtool Row
#' @description For a given row of fdrtool_grid, checks whether hyperparameter combinations
#' can be used for fdrtool on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param fdrtool_grid data.frame, each row is a possible set of hyperparameters for fdrtool
#' @param row integer, row of fdrtool_grid considered
#' @param lower_pi0 double, hyperparameter combinations that give pi0 estimates below this are excluded from the grid
#'
#' @return boolean, whether the hyperparameter combination run without error
check_fdrtool_row <- function(
  test_statistics,
  fdrtool_grid,
  row,
  lower_pi0
) {
  # for each row, attempt to run fdrtool
  run_i <- run_fdrtool_row(
    test_statistics = test_statistics,
    fdrtool_grid = fdrtool_grid,
    row = row,
    returnFdr = FALSE
  )

  # keep this row in final grid if it ran without error,
  # and if it does not estimate fdrs as all 0 or all 1
  if (!is.null(run_i)) {
    if (run_i$pi0 <= 1 & run_i$pi0 >= lower_pi0) {
      return(TRUE)
    }
  }

  return(FALSE)
}

#' @title Reduce fdrtool grid
#' @description Reduces grid to parameter combinations
#' that can run fdrtool on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param fdrtool_grid data.frame where each row is a possible set of hyperparameters for fdrtool
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#' @param verbose boolean
#'
#' @return reduced fdrtool_grid: dataframe where each row is a possible set of hyperparameters for fdrtool
reduce_fdrtool_grid <- function(
  test_statistics,
  fdrtool_grid,
  lower_pi0,
  parallel_param = NULL,
  verbose = FALSE
) {

  ok_rows <- unlist(parlapply(
    X = 1:nrow(fdrtool_grid),
    parallel_param = parallel_param,
    FUN = function(i, run_fdrtool_row, fdrtool_grid, test_statistics) {
      if(check_fdrtool_row(
        test_statistics = test_statistics,
        fdrtool_grid = fdrtool_grid,
        row = i,
        lower_pi0 = lower_pi0
      )) {
        return(i)
      }
    },
    run_fdrtool_row = run_fdrtool_row,
    fdrtool_grid = fdrtool_grid,
    test_statistics = test_statistics
  ))

  if(verbose) {
    print(paste(
      length(ok_rows),'/',nrow(fdrtool_grid),
      'fdrtool parameter sets are usable'
    ))
  }

  return (fdrtool_grid[ok_rows,])
}

#' @title Build fdrtool Grid
#' @description Build a factorial grid of possible hyperparameters for fdrtool and
#' reduce to hyperparameter sets that can be run on provided data without error.
#'
#' @param test_statistics vector, test statistics
#' @param cutoff.method vector, options for cutoff method hyperparameter.
#' `cutoff.method` is one of "fndr" (default), "pct0", "locfdr".
#' @param pct0_range vector c(min, max), range for pct0 hyperparameter. `pct0` is the
#' fraction of data used for fitting null model - only if cutoff.method="pct0".
#' @param grid_size integer, desired size of grid to use. Note that this is *not a guarantee*
#' as some hyperparameter combinations may fail when run on the data.
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#' @param method string, one of c('random', 'grid'). 'random' will sample from a uniform
#' distribution within the ranges, 'grid' will select equally spaced values.
#' @param seed integer, random seed used if method = 'random'
#' @param verbose boolean
#'
#' @return data.frame, each row is a possible set of hyperparameters for fdrtool
#' @export
build_fdrtool_grid <- function(
  test_statistics,
  cutoff.method = c('fndr','pct0','locfdr'),
  pct0_range = c(0, 0.5),
  grid_size = 40,
  parallel_param = NULL,
  method = 'grid',
  lower_pi0 = 0.7,
  seed = NULL,
  verbose = FALSE
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  #pct0 parameter is only used when cutoff.method = 'pct0'
  if ('pct0' %in% cutoff.method) {
    other_cutoff_methods <- cutoff.method[cutoff.method != 'pct0']
    n_pct0 <- grid_size - length(other_cutoff_methods)
    if (method == 'random') {
      pct0_vec <- c(
        rep(0.75, grid_size - n_pct0),
        runif(n = n_pct0, min = pct0_range[1], max = pct0_range[2])
      )
    } else if (method == 'grid') {
      pct0_vec <- c(
        rep(0.75, grid_size - n_pct0),
        seq(from = pct0_range[1], to = pct0_range[2], length.out = n_pct0)
      )
    }
    fdrtool_grid <- data.frame(
      cutoff.method = c(
        other_cutoff_methods,
        rep('pct0', n_pct0)
      ),
      pct0 = pct0_vec
    )
  } else {
    fdrtool_grid <- data.frame(
      cutoff.method = cutoff.method,
      pct0 = rep(0.75, length(cutoff.method))
    )
    if (verbose & nrow(fdrtool_grid) <  grid_size) {
      message(paste0(
        "Couldn't create grid size ", grid_size,
        ". Without 'pct0' as a cutoff.method option, there are only ",
        nrow(fdrtool_grid), " possible hyperparameter combinations. Using ",
        "a grid size of ", nrow(fdrtool_grid), "."
      ))
    }
  }

  fdrtool_grid_reduced <- reduce_fdrtool_grid(
    test_statistics = test_statistics,
    fdrtool_grid = fdrtool_grid,
    lower_pi0 = lower_pi0,
    parallel_param = parallel_param,
    verbose = verbose
  )
  if (verbose & nrow(fdrtool_grid_reduced) < nrow(fdrtool_grid)) {
    message(paste0(
      nrow(fdrtool_grid) - nrow(fdrtool_grid_reduced), " hyperparameter combinations failed",
      ifelse(lower_pi0 > 0, paset(' or had pi0 below', lower_pi0), ''),
      " when run on the data. Using a grid size of ", nrow(fdrtool_grid_reduced), "."
    ))
  }

  return(fdrtool_grid_reduced)
}
