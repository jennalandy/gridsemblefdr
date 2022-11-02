#' @title Check Fdrtool Row
#' @description For a given row of fdrtool_grid, checks whether hyperparameter combinations
#' can be used for fdrtool on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param fdrtool_grid data.frame, each row is a possible set of hyperparameters for fdrtool
#' @param row integer, row of fdrtool_grid considered
#' @param lower_pi0 double, exclude hyperparameter combinations that give pi0 estimates below this threshold
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
    if (run_i$pi0 > 1) {
      return('pi0 > 1')
    } else if (run_i$pi0 < lower_pi0) {
      return('pi0 < lower_pi0')
    } else {
      return(TRUE)
    }
  }

  return('error')
}

#' @title Reduce fdrtool grid
#' @description Reduces grid to parameter combinations
#' that can run fdrtool on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param fdrtool_grid data.frame where each row is a possible set of hyperparameters for fdrtool
#' @param lower_pi0 double, exclude hyperparameter combinations that give pi0 estimates below this threshold
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

  checks <- unlist(parlapply(
    X = 1:nrow(fdrtool_grid),
    parallel_param = parallel_param,
    FUN = function(i, run_fdrtool_row, fdrtool_grid, test_statistics) {
      check = check_fdrtool_row(
        test_statistics = test_statistics,
        fdrtool_grid = fdrtool_grid,
        row = i,
        lower_pi0 = lower_pi0
      )

      if (check == TRUE) {
        return(i)
      } else {
        return(check)
      }
    },
    run_fdrtool_row = run_fdrtool_row,
    fdrtool_grid = fdrtool_grid,
    test_statistics = test_statistics
  ))

  rows = grepl('^[[:digit:]]+$', checks)
  ok_rows = as.numeric(checks[rows])
  fails = checks[!rows]

  if(verbose) {
    message(paste0(
      '\t',length(ok_rows),'/',nrow(fdrtool_grid),
      ' fdrtool thetas included'
    ))
  }

  if (verbose & length(fails) > 0) {
    fails_tab = table(fails)
    for (reason in names(fails_tab)) {
      message(paste0(
        '\t\t', fails_tab[reason],
        ' thetas dropped for ',
        reason
      ))
    }
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
#' @param grid_size integer, maximum size of grid to use. Note that this is *not the final grid size*
#' as some hyperparameter combinations may fail when run on the data.
#' @param lower_pi0 double, exclude hyperparameter combinations that give pi0 estimates below this threshold
#' @param method string, one of c('random', 'grid'). 'random' will sample from a uniform
#' distribution within the ranges, 'grid' will select equally spaced values.
#' @param seed integer, random seed used if method = 'random'
#' @param parallel boolean, whether to utilize parallelization
#' @param n_workers integer, number of cores to use if parallel, default 2 less than available.
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#' @param verbose boolean
#'
#' @importFrom stats runif
#' @return data.frame, each row is a possible set of hyperparameters for fdrtool
#' @export
build_fdrtool_grid <- function(
  test_statistics,
  cutoff.method = c('fndr','pct0','locfdr'),
  pct0_range = c(0, 0.5),
  grid_size = 40,
  lower_pi0 = 0.7,
  method = 'grid',
  seed = NULL,
  parallel = min(TRUE, n_workers > 1),
  n_workers = max(parallel::detectCores() - 2, 1),
  parallel_param = NULL,
  verbose = FALSE
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (parallel & is.null(parallel_param) & n_workers > 1) {
    parallel_param = BiocParallel::MulticoreParam(
      workers = n_workers,
      tasks = n_workers
    )
  }
  if (!is.null(parallel_param) & verbose) {
    message('Building fdrtool grid in parallel')
  } else if (verbose) {
    message('Building fdrtool grid')
  }

  #pct0 parameter is only used when cutoff.method = 'pct0'
  if ('pct0' %in% cutoff.method) {
    other_cutoff_methods <- cutoff.method[cutoff.method != 'pct0']
    n_pct0 <- grid_size - length(other_cutoff_methods)
    if (method == 'random') {
      pct0_vec <- c(
        rep(0.75, grid_size - n_pct0),
        stats::runif(n = n_pct0, min = pct0_range[1], max = pct0_range[2])
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
        "\tCouldn't create grid size ", grid_size,
        ". Without 'pct0' as a cutoff.method option, there are only ",
        nrow(fdrtool_grid), " possible thetas. Using ",
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

  return(fdrtool_grid_reduced)
}
