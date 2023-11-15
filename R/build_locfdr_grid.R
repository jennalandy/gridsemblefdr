#' Check locfdr Row
#' @description For a given row of locfdr_grid, checks whether hyperparameter
#' combinations can be used for locfdr on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param locfdr_grid data.frame, rows are possible hyperparameters for locfdr
#' @param row integer, row of locfdr_grid considered
#'
#' @return boolean, whether the hyperparameter combination ran without error
#' @noRd
check_locfdr_row <- function(
  test_statistics,
  locfdr_grid,
  row
) {
  # for each row, attempt to run fdrtool
  run_i <- run_locfdr_row(
    test_statistics = test_statistics,
    locfdr_grid = locfdr_grid,
    row = row,
    returnFdr = FALSE
  )

  # keep this row in final grid if it ran without error,
  # and if it does not estimate fdrs as all 0 or all 1
  if (!is.null(run_i)) {
    if (run_i$pi0 > 1) {
      return('pi0 > 1')
    } else {
      return(TRUE)
    }
  }

  return('error')
}

#' @title Reduce locfdr grid
#' @description Reduces grid to parameter combinations
#' that can run locfdr on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param locfdr_grid data.frame, rows are possible hyperparameters for locfdr
#' @param parallel_param BiocParallel object
#' @param verbose boolean
#'
#' @return data.frame, each row is a possible set of hyperparameters for locfdr
#' @noRd
reduce_locfdr_grid <- function(
  test_statistics,
  locfdr_grid,
  parallel_param = NULL,
  verbose = FALSE
) {

  checks <- unlist(parlapply(
    X = seq_len(nrow(locfdr_grid)),
    parallel_param = parallel_param,
    FUN = function(i, run_locfdr_row, locfdr_grid, test_statistics) {
      check = check_locfdr_row(
        test_statistics,
        locfdr_grid,
        row = i
      )

      if (check == TRUE) {
        return(i)
      } else {
        return(check)
      }
    },
    run_locfdr_row = run_locfdr_row,
    locfdr_grid = locfdr_grid,
    test_statistics = test_statistics
  ))

  rows = grepl('^[[:digit:]]+$', checks)
  ok_rows = as.numeric(checks[rows])
  fails = checks[!rows]

  if(verbose) {
    message(paste0(
      '\t',length(ok_rows),'/',nrow(locfdr_grid),
      ' locfdr thetas included'
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

  return (locfdr_grid[ok_rows,])
}

#' @title Build locfdr Grid
#' @description Build a factorial grid of hyperparameters for locfdr and
#' reduce to hyperparameter sets that can be run on provided data without error.
#'
#' @param test_statistics vector, test statistics
#' @param pct_range vector c(min, max), range for pct hyperparameter. `pct` is
#' the excluded tail proportions of zz's when fitting f(z).
#' pct=0 includes full range of zz's.
#' @param pct0_range vector c(min, max), range for pct0 hyperparameter. `pct0`
#' is the proportion of the zz distribution used in fitting the null density
#' f0(z) by central matching.
#' @param nulltype vector, options for nulltype hyperparameter. `nulltype`
#' is the type of null hypothesis assumed in estimating f0(z), for
#' use in the fdr calculations. 0 is the theoretical null N(0,1),
#' 1 is maximum likelihood estimation, 2 is central matching estimation,
#' 3 is a split normal version of 2.
#' @param type vector, options for type hyperparameter. `type` is the type of
#' fitting used for f; 0 is a natural spline, 1 is a polynomial.
#' @param grid_size integer, maximum size of grid to use. Note that this is
#' *not the final grid size*, combinations may fail when run on the data.
#' @param method string, one of c('random', 'grid'). 'random' will sample
#' uniformly within the ranges, 'grid' will select equally spaced values
#'
#' @param parallel boolean, whether to utilize parallelization
#' @param n_workers integer, number of cores to use if parallel
#' @param parallel_param BiocParallel object
#' @param verbose boolean
#'
#' @importFrom stats runif
#' @return data.frame, each row is a possible set of hyperparameters for locfdr
#' @export
#' @examples
#' set.seed(123)
#' test_statistics = c(rnorm(800), runif(100, -10, -5), runif(100, 5, 10))
#' locfdr_grid = build_locfdr_grid(test_statistics)
build_locfdr_grid <- function(
  test_statistics, pct_range = c(0, 0.3), pct0_range = c(0, 0.45),
  nulltype = c(1,2,3), type = c(0,1), grid_size = 40,
  method = 'grid', parallel_param = NULL, parallel = min(TRUE, n_workers > 1),
  n_workers = max(parallel::detectCores() - 2, 1), verbose = FALSE
) {

  if (parallel & is.null(parallel_param) & n_workers > 1) {
    parallel_param = BiocParallel::MulticoreParam(
      workers = n_workers, tasks = n_workers
    )
  }
  if (!is.null(parallel_param) & verbose) {
    message('Building locfdr grid in parallel')
  } else if (verbose) { message('Building locfdr grid') }

  if (method == 'random') { # select uniformly within ranges
   locfdr_grid <- data.frame(
     pct = stats::runif(
       n = grid_size, min = pct_range[1], max = pct_range[2]
     ),
     pct0 = stats::runif(
       n = grid_size, min = pct0_range[1], max = pct0_range[2]
     ),
     nulltype = sample(nulltype, size = grid_size, replace = TRUE),
     type = sample(type, size = grid_size, replace = TRUE)
   )
  } else if (method == 'grid') { # select equally spaced within ranges
    # calculate number of values for pct and pct0 to each expand grid with
    nG <- ceiling(sqrt(grid_size/(length(nulltype)*length(type))))
    locfdr_grid <- expand.grid(
      pct = seq(from = pct_range[1], to = pct_range[2], length.out = nG),
      pct0 = seq(from = pct0_range[1], to = pct0_range[2], length.out = nG),
      nulltype = nulltype,
      type = type
    )
  }

  if (verbose & nrow(locfdr_grid) > grid_size & method == 'grid') {
    message(paste0(
      '\t', nrow(locfdr_grid) - grid_size, " extra thetas needed for full grid"
    ))
  }

  locfdr_grid_reduced <- reduce_locfdr_grid(
    test_statistics = test_statistics,
    locfdr_grid = locfdr_grid,
    parallel_param = parallel_param,
    verbose = verbose
  )

  return(locfdr_grid_reduced)
}
