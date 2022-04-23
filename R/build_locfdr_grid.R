#' @title Reduce locfdr grid
#' @description Reduces grid to parameter combinations
#' that can run locfdr on the provided data without error
#'
#' @param test_statistics vector of test statistics
#' @param locfdr_grid data frame where each row is a possible set of hyperparameters for locfdr
#' @param parallel_param if not NULL, processes are run in parallel with this parameter
#' @param verbose if TRUE, status updates will be displayed
#'
#' @return dataframe where each row is a possible set of hyperparameters for locfdr
#'
#' @export
reduce_locfdr_grid <- function(
  test_statistics,
  locfdr_grid,
  parallel_param = NULL,
  verbose = FALSE
) {

  ok_rows <- unlist(parlapply(
    X = 1:nrow(locfdr_grid),
    parallel_param = parallel_param,
    FUN = function(i, run_locfdr_row, locfdr_grid, test_statistics) {
      # for each row, attempt to run locfdr
      run_i <- run_locfdr_row(
        test_statistics = test_statistics,
        locfdr_grid = locfdr_grid,
        row = i,
        returnFdr = FALSE
      )

      # keep this row in final grid if it ran without error,
      # and if it does not estimate fdrs as all 0 or all 1
      if (!is.null(run_i)) {
        if (run_i$pi0 <= 1) {
          return(i)
        }
      }
    },
    run_locfdr_row = run_locfdr_row,
    locfdr_grid = locfdr_grid,
    test_statistics = test_statistics
  ))

  if(verbose) {
    print(paste(
      length(ok_rows),'/',nrow(locfdr_grid),
      'locfdr parameter sets are usable'
    ))
  }

  return (locfdr_grid[ok_rows,])
}

#' @title Build locfdr Grid
#' @description Build a grid of possible hyperparameters for locfdr
#' from separate vectors of hyperparameter options. Final grid only considers
#' hyperparameters that can be run on provided data without error.
#'
#' @param test_statistics vector of test statistics
#' @param pct vector of options for pct hyperparameter. `pct` is the
#' excluded tail proportions of zz's when fitting f(z).
#' pct=0 includes full range of zz's.
#' @param pct0 vector of options for pct0 hyperparameter. `pct0` is the
#' proportion of the zz distribution used in fitting the null density
#' f0(z) by central matching. Scalar, range [pct0, 1-pct0] is used.
#' @param nulltype vector of options for nulltype hyperparameter. `nulltype`
#' is the type of null hypothesis assumed in estimating f0(z), for
#' use in the fdr calculations. 0 is the theoretical null N(0,1),
#' 1 is maximum likelihood estimation, 2 is central matching estimation,
#' 3 is a split normal version of 2.
#' @param type vector of options for type hyperparameter. `type` is the type of
#' fitting used for f; 0 is a natural spline, 1 is a polynomial.
#' @param parallel_param if not NULL, processes are run in parallel with this parameter
#' @param verbose if TRUE, status updates will be displayed
#'
#' @return dataframe where each row is a possible set of hyperparameters for
#' locfdr on a specific set of test statistics
#' @export
build_locfdr_grid <- function(
  test_statistics,
  pct = c(0,2,4,6)/1000,
  pct0 = 1/c(3:10),
  nulltype = 1:3,
  type = 0:1,
  parallel_param = NULL,
  verbose = FALSE
) {

  locfdr_grid <- expand.grid(
    pct = pct,
    pct0 = pct0,
    nulltype = nulltype,
    type = type
  )

  # pct0 is only used when central matching (nulltype = 2 or 3) is used
  # limit to default_pct0 value if it won't be used
  default_pct0 <- ifelse(1/4 %in% unique(pct0), 1/4, unique(pct0)[1] )
  locfdr_grid <- locfdr_grid[
    locfdr_grid$nulltype == 2 |
      locfdr_grid$nulltype == 3 |
      locfdr_grid$pct0 == default_pct0,
  ]

  locfdr_grid_reduced <- reduce_locfdr_grid(
    test_statistics = test_statistics,
    locfdr_grid = locfdr_grid,
    parallel_param = parallel_param,
    verbose = verbose
  )

  return(locfdr_grid_reduced)
}
