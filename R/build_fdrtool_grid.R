#' @title Check Fdrtool Row
#' @description For a given row of fdrtool_grid, checks whether hyperparameter
#' combinations can be used for fdrtool on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param fdrtool_grid data.frame, rows are possible hyperparameters for fdrtool
#' @param row integer, row of fdrtool_grid considered
#'
#' @return boolean, whether the hyperparameter combination run without error
#' @noRd
check_fdrtool_row <- function(
  test_statistics,
  fdrtool_grid,
  row
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
    }  else {
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
#' @param fdrtool_grid data.frame, rows are possible hyperparameters for fdrtool
#' @param parallel_param BiocParallel object
#' @param verbose boolean
#'
#' @return reduced fdrtool_grid: dataframe, rows are hyperparameters for fdrtool
#' @noRd
reduce_fdrtool_grid <- function(
  test_statistics,
  fdrtool_grid,
  parallel_param = NULL,
  verbose = FALSE
) {

  checks <- unlist(parlapply(
    X = seq_len(nrow(fdrtool_grid)),
    parallel_param = parallel_param,
    FUN = function(i, run_fdrtool_row, fdrtool_grid, test_statistics) {
      check = check_fdrtool_row(
        test_statistics = test_statistics,
        fdrtool_grid = fdrtool_grid,
        row = i
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
      ' fdrtool models included'
    ))
  }

  if (verbose & length(fails) > 0) {
    fails_tab = table(fails)
    for (reason in names(fails_tab)) {
      message(paste0(
        '\t\t', fails_tab[reason],
        ' models dropped for ',
        reason
      ))
    }
  }

  return (fdrtool_grid[ok_rows,])
}

#' @title Build fdrtool Grid
#' @description Build a grid of possible hyperparameters for fdrtool and reduce
#' to hyperparameter sets that can be run on provided data without error.
#'
#' @param test_statistics vector, test statistics
#' @param cutoff.method vector, options for cutoff method hyperparameter.
#' `cutoff.method` is one of "fndr" (default), "pct0", "locfdr".
#' @param pct0_range vector c(min, max), range of values for pct0 hyperparameter, the
#' fraction of data used for fitting null model - only if cutoff.method="pct0".
#' Default in `fdrtool` package is `pct0 = 0.75`.
#' @param grid_depth integer, number of evenly-spaced values of continuous parameters
#' considered within their respective `_range`.
#' @param parallel boolean, whether to utilize parallelization
#' @param n_workers integer, number of cores to use if parallel
#' @param parallel_param BiocParallel object
#' @param verbose boolean
#'
#' @importFrom stats runif
#' @return data.frame, each row is a possible set of hyperparameters for fdrtool
#' @export
#'
#' @examples
#' set.seed(123)
#' test_statistics = c(rnorm(800), runif(100, -10, -5), runif(100, 5, 10))
#' fdrtool_grid = build_fdrtool_grid(test_statistics)
build_fdrtool_grid <- function(
  test_statistics, cutoff.method = c('fndr','pct0','locfdr'),
  pct0_range = c(0.4, 1), grid_depth = 5,
  parallel = min(TRUE, n_workers > 1),
  n_workers = max(parallel::detectCores() - 2, 1), parallel_param = NULL,
  verbose = FALSE
) {
  pct0 = seq(pct0_range[1], pct0_range[2], length.out = grid_depth)

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

  fdrtool_grid <- expand.grid(
    cutoff.method = cutoff.method[cutoff.method != 'pct0'],
    pct0 = 0.75 # this value doesn't affect the model
  )

  if ('pct0' %in% cutoff.method) {
    fdrtool_grid = rbind(
      fdrtool_grid,
      expand.grid(
        cutoff.method = c('pct0'),
        pct0 = pct0
      )
    )
  }

  fdrtool_grid_reduced <- reduce_fdrtool_grid(
    test_statistics = test_statistics, fdrtool_grid = fdrtool_grid,
    parallel_param = parallel_param, verbose = verbose
  )

  return(fdrtool_grid_reduced)
}
