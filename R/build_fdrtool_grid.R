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
#' @description Build a grid of possible hyperparameters for fdrtool and reduce
#' to hyperparameter sets that can be run on provided data without error.
#'
#' @param test_statistics vector, test statistics
#' @param cutoff.method vector, options for cutoff method hyperparameter.
#' `cutoff.method` is one of "fndr" (default), "pct0", "locfdr".
#' @param pct0_range vector c(min, max), range for pct0 hyperparameter, the
#' fraction of data used for fitting null model - only if cutoff.method="pct0".
#' @param grid_size integer, maximum size of grid to use. Note that this is
#' *not the final grid size*, combinations may fail when run on the data.
#' @param method string, one of c('random', 'grid'). 'random' will sample
#' uniformly within the ranges, 'grid' will select equally spaced values
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
  pct0_range = c(0,0.8), grid_size = 40, method = 'grid',
  parallel = min(TRUE, n_workers > 1),
  n_workers = max(parallel::detectCores() - 2, 1), parallel_param = NULL,
  verbose = FALSE
) {
  if (parallel & is.null(parallel_param) & n_workers > 1) {
    parallel_param = BiocParallel::MulticoreParam(
      workers = n_workers,
      tasks = n_workers
    )
  }
  if (!is.null(parallel_param) & verbose) {
    message('Building fdrtool grid in parallel')
  } else if (verbose) {  message('Building fdrtool grid') }

  if (method == 'random') {
    fdrtool_grid <- data.frame(
      cutoff.method = sample(cutoff.method, size = grid_size, replace = TRUE),
      pct0 = stats::runif(
        n = grid_size, min = pct0_range[1], max = pct0_range[2]
      )
    )
  } else if (method == 'grid') {
    non <- cutoff.method[cutoff.method != 'pct0']
    # pct0 method uses pct0 parameter, others dont
    n_pct = grid_size - length(non)

    fdrtool_grid <- rbind(
      expand.grid(
        cutoff.method = c('pct0'),
        pct0 = seq(
          from = pct0_range[1], to = pct0_range[2], length.out = n_pct
        )
      ),
      expand.grid(
        cutoff.method = non,
        pct0 = 0.75
      )
    )
  } else {stop('method must be one of c("random","grid")')}

  if (verbose & nrow(fdrtool_grid) < grid_size & !('pct0' %in% cutoff.method)) {
    message(paste0(
      "\tWithout 'pct0' as a cutoff.method option, there are only ",
      nrow(fdrtool_grid), " possible thetas."
    ))
  } else if (verbose & nrow(fdrtool_grid) > grid_size & method == 'grid') {
    message(paste0(
      '\t', nrow(fdrtool_grid) - grid_size, " extra thetas needed for full grid"
    ))
  }

  fdrtool_grid_reduced <- reduce_fdrtool_grid(
    test_statistics = test_statistics, fdrtool_grid = fdrtool_grid,
    parallel_param = parallel_param, verbose = verbose
  )

  return(fdrtool_grid_reduced)
}
