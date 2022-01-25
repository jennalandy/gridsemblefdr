#' @title Reduce locfdr grid
#' @description Reduces grid to parameter combinations
#' that can run locfdr on the provided data without error
#'
#' @param t vector of test statistics
#' @param locfdr_grid data frame where each
#' @param verbose
#'
#' @return dataframe where each row is a possible set of hyperparameters for locfdr
reduce_locfdr_grid <- function(t, locfdr_grid, verbose = FALSE) {

  ok_rows <- c()
  for (i in 1:nrow(locfdr_grid)) {
    # for each row, attempt to run locfdr
    run_i <- run_locfdr_row(
      t = t,
      locfdr_grid = locfdr_grid,
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
#' @param t vector of test statistics
#' @param pct vector of options for pct hyperparameter
#' @param pct0 vector of options for pct0 hyperparameter
#' @param nulltype vector of options for nulltype hyperparameter
#' @param type vector of options for type hyperparameter
#' @param verbose
#'
#' @return dataframe where each row is a possible set of hyperparameters for
#' locfdr on a specific set of test statistics
#' @export
build_locfdr_grid <- function(
  t,
  pct = c(0,2,4,6)/1000,
  pct0 = 1/c(3:10),
  nulltype = 1:3,
  type = 0:1,
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
    t = t,
    locfdr_grid = locfdr_grid,
    verbose = verbose
  )

  return(locfdr_grid_reduced)
}
