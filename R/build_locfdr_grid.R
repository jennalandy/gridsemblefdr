#' Check locfdr Row
#' @description For a given row of locfdr_grid, checks whether hyperparameter
#' combinations can be used for locfdr on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param locfdr_grid data.frame, rows are possible hyperparameters for locfdr
#' @param row integer, row of locfdr_grid considered
#' @param drop_pi0_1 boolean, whether to discard models that estimate pi0 = 1
#'
#' @return boolean, whether the hyperparameter combination ran without error
#' @noRd
check_locfdr_row <- function(
  test_statistics,
  locfdr_grid,
  row,
  drop_pi0_1 = TRUE
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
    } else if (drop_pi0_1 & run_i$pi0 == 1) {
      return('pi0 == 1')
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
#' @param drop_pi0_1 boolean, whether to discard models that estimate pi0 = 1
#' @param parallel_param BiocParallel object
#' @param verbose boolean
#'
#' @return data.frame, each row is a possible set of hyperparameters for locfdr
#' @export
reduce_locfdr_grid <- function(
  test_statistics,
  locfdr_grid,
  drop_pi0_1 = TRUE,
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
        row = i,
        drop_pi0_1 = drop_pi0_1
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
      ' locfdr models included'
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
#' @param grid_depth integer, number of evenly-spaced values of continuous parameters
#' considered within their respective `_range`.
#' @param drop_pi0_1 boolean, whether to discard models that estimate pi0 = 1
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
  test_statistics, pct_range = c(0, 0.2), pct0_range = c(0, 0.3),
  nulltype = c(1,2,3), type = c(0,1), grid_depth = 5,
  parallel_param = NULL, parallel = min(TRUE, n_workers > 1),
  drop_pi0_1 = TRUE,
  n_workers = max(parallel::detectCores() - 2, 1), verbose = FALSE
) {
  pct = seq(pct_range[1], pct_range[2], length.out = grid_depth)
  pct0 = seq(pct0_range[1], pct0_range[2], length.out = grid_depth)

  if (parallel & is.null(parallel_param) & n_workers > 1) {
    parallel_param = BiocParallel::MulticoreParam(
      workers = n_workers, tasks = n_workers
    )
  }

  if (!is.null(parallel_param) & verbose) {
    message('Building locfdr grid in parallel')
  } else if (verbose) {
    message('Building locfdr grid')
  }

  locfdr_grid <- expand.grid(
    pct = pct,
    pct0 = pct0,
    nulltype = nulltype,
    type = type
  )

  locfdr_grid_reduced <- reduce_locfdr_grid(
    test_statistics = test_statistics,
    locfdr_grid = locfdr_grid,
    drop_pi0_1 = drop_pi0_1,
    parallel_param = parallel_param,
    verbose = verbose
  )

  return(locfdr_grid_reduced)
}
