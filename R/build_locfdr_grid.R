#' Check locfdr Row
#' For a given row of locfdr_grid, checks whether hyperparameter combinations
#' can be used for locfdr on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param locfdr_grid data.frame, each row is a possible set of hyperparameters for locfdr
#' @param row integer, row of locfdr_grid considered
#' @param lower_pi0 double, exclude hyperparameter combinations that give pi0 estimates below this threshold
#'
#' @return boolean, whether the hyperparameter combination ran without error
check_locfdr_row <- function(
  test_statistics,
  locfdr_grid,
  row,
  lower_pi0
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
    if (run_i$pi0 <= 1 & run_i$pi0 >= lower_pi0) {
      return(TRUE)
    }
  }

  return(FALSE)
}

#' @title Reduce locfdr grid
#' @description Reduces grid to parameter combinations
#' that can run locfdr on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param locfdr_grid data.frame, each row is a possible set of hyperparameters for locfdr
#' @param lower_pi0 double, exclude hyperparameter combinations that give pi0 estimates below this threshold
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#' @param verbose boolean
#'
#' @return data.frame, each row is a possible set of hyperparameters for locfdr
reduce_locfdr_grid <- function(
  test_statistics,
  locfdr_grid,
  lower_pi0,
  parallel_param = NULL,
  verbose = FALSE
) {

  ok_rows <- unlist(parlapply(
    X = 1:nrow(locfdr_grid),
    parallel_param = parallel_param,
    FUN = function(i, run_locfdr_row, locfdr_grid, test_statistics) {
      if (check_locfdr_row(
        test_statistics,
        locfdr_grid,
        row = i,
        lower_pi0 = lower_pi0
      )) {
        return(i)
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
#' @description Build a factorial grid of possible hyperparameters for locfdr and
#' reduce to hyperparameter sets that can be run on provided data without error.
#'
#' @param test_statistics vector, test statistics
#' @param pct_range vector c(min, max), range for pct hyperparameter. `pct` is the
#' excluded tail proportions of zz's when fitting f(z).
#' pct=0 includes full range of zz's.
#' @param pct0_range vector c(min, max), range for pct0 hyperparameter. `pct0` is the
#' proportion of the zz distribution used in fitting the null density
#' f0(z) by central matching.
#' @param nulltype vector, options for nulltype hyperparameter. `nulltype`
#' is the type of null hypothesis assumed in estimating f0(z), for
#' use in the fdr calculations. 0 is the theoretical null N(0,1),
#' 1 is maximum likelihood estimation, 2 is central matching estimation,
#' 3 is a split normal version of 2.
#' @param type vector, options for type hyperparameter. `type` is the type of
#' fitting used for f; 0 is a natural spline, 1 is a polynomial.
#' @param grid_size integer, maximum size of grid to use. Note that this is *not the final grid size*
#' as some hyperparameter combinations may fail when run on the data.
#' @param lower_pi0 double, exclude hyperparameter combinations that give pi0 estimates below this threshold
#' @param method string, one of c('random', 'grid'). 'random' will sample from a uniform
#' distribution within the ranges, 'grid' will select equally spaced values.
#' @param seed integer, random seed used if method = 'random'
#'
#' @param parallel boolean, whether to utilize parallelization
#' @param n_workers integer, number of cores to use if parallel, default 2 less than available.
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#' @param verbose boolean
#'
#' @importFrom stats runif
#' @return data.frame, each row is a possible set of hyperparameters for locfdr
#' @export
build_locfdr_grid <- function(
  test_statistics,
  pct_range = c(0, 0.1),
  pct0_range = c(0, 1/3),
  nulltype = 1:3,
  type = 0:1,
  grid_size = 40,
  lower_pi0 = 0.7,
  method = 'grid',
  seed = NULL,
  parallel_param = NULL,
  parallel = TRUE,
  n_workers = max(parallel::detectCores() - 2, 1),
  verbose = FALSE
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (parallel & is.null(parallel_param)) {
    parallel_param = BiocParallel::MulticoreParam(
      workers = n_workers,
      tasks = n_workers
    )
  }

  # if nulltypes 2 or 3 are options,
  # make sure that parameter combinations are not redundant
  # i.e. nulltype 1, pct0 = 0.3 vs 0.4 doesn't matter
  nulltypes_use_pct0 = nulltype[nulltype == 2 | nulltype == 3]
  nulltype_no_pct0 = nulltype[nulltype != 2 & nulltype != 3]
  prop_nulltype_no_pct0 = length(nulltype_no_pct0) / length(nulltype)
  n_nulltype_no_pct0 = round(grid_size*prop_nulltype_no_pct0)
  n_nulltypes_use_pct0 = grid_size - n_nulltype_no_pct0

  # if random method, select from uniform distribution within ranges
  if (method == 'random') {

    # pct0 parameter is only used if nulltype is 2 or 3
    if (2 %in% nulltype | 3 %in% nulltype) {

      locfdr_grid <- rbind(
        # n_nulltypes_use_pct0 combinations
        # nulltype is selected from nulltypes_use_pct0
        # pct0 has multiple options
        data.frame(
          pct = stats::runif(n = n_nulltypes_use_pct0, min = pct_range[1], max = pct_range[2]),
          pct0 = stats::runif(n = n_nulltypes_use_pct0, min = pct0_range[1], max = pct0_range[2]),
          nulltype = sample(nulltypes_use_pct0, size = n_nulltypes_use_pct0, replace = TRUE),
          type = sample(type, size = n_nulltypes_use_pct0, replace = TRUE)
        ),
        # n_nulltype_no_pct0 combinations
        # nulltype is selected from nulltype_no_pct0
        # pct0 has default value
        data.frame(
          pct = stats::runif(n = n_nulltype_no_pct0, min = pct_range[1], max = pct_range[2]),
          pct0 = rep(1/4, n_nulltype_no_pct0), # won't be used, keep as default
          nulltype = sample(nulltype_no_pct0, size = n_nulltype_no_pct0, replace = TRUE),
          type = sample(type, size = n_nulltype_no_pct0, replace = TRUE)
        )
      )

    } else {

      # if nulltypes 2 or 3 are not options,
      # pct0 never matters, keep as default value
      locfdr_grid <- data.frame(
        pct = stats::runif(n = grid_size, min = pct_range[1], max = pct_range[2]),
        pct0 = rep(1/4, grid_size), # won't be used, keep as default
        nulltype = sample(nulltype, size = grid_size, replace = TRUE),
        type = sample(type, size = grid_size, replace = TRUE)
      )
    }

    # if grid method, select equally spaced values within ranges
  } else if (method == 'grid') {

    # pct0 parameter is only used if nulltype is 2 or 3
    if (2 %in% nulltype | 3 %in% nulltype) {

      nulltypes_use_pct0 = nulltype[nulltype == 2 | nulltype == 3]
      nulltype_no_pct0 = nulltype[nulltype != 2 & nulltype != 3]
      prop_nulltype_no_pct0 = length(nulltype_no_pct0) / length(nulltype)
      n_nulltype_no_pct0 = round(grid_size*prop_nulltype_no_pct0)
      n_nulltypes_use_pct0 = grid_size - n_nulltype_no_pct0

      # make sure expanded grid length is as close as possible to grid_size
      # n_combos_* = number combinations before numeric vars are samples from uniform in group *
      # seq_length_* = length of remaining sequences to make expanded grid correct size

      n_combos_use_pct0 = length(type)*length(nulltypes_use_pct0)
      # sqrt is because there are two sequences that need to be added to this one
      seq_length_use_pct0 = ceiling(sqrt(n_nulltypes_use_pct0/n_combos_use_pct0))

      n_combos_no_pct0 = length(type)*length(nulltype_no_pct0)
      # no sqrt because there is only one sequence that needs to be added to this one
      seq_length_no_pct0 = ceiling(n_nulltype_no_pct0/n_combos_no_pct0)

      locfdr_grid <- rbind(
        expand.grid(
          pct = seq(from = pct_range[1], to = pct_range[2], length.out = seq_length_use_pct0),
          pct0 = seq(from = pct0_range[1], to = pct0_range[2], length.out = seq_length_use_pct0),
          nulltype = nulltypes_use_pct0,
          type = type
        ),
        expand.grid(
          pct = seq(from = pct_range[1], to = pct_range[2], length.out = seq_length_no_pct0),
          pct0 = 1/4, # won't be used, keep as default
          nulltype = nulltype_no_pct0,
          type = type
        )
      )
    } else {
      # make sure expanded grid length is as close as possible to grid_size
      # n_combos = number combinations before numeric vars are samples from uniform
      # seq_length = length of remaining sequences to make expanded grid correct size

      n_combos = length(type)*length(nulltype)
      seq_length = ceiling(grid_size/n_combos)

      locfdr_grid <- expand.grid(
        pct = seq(from = pct_range[1], to = pct_range[2], length.out = seq_length),
        pct0 = 1/4, # won't be used, keep as default
        nulltype = nulltype,
        type = type
      )
    }
  }

  if (verbose & nrow(locfdr_grid) > grid_size & method == 'grid') {
    message(paste0(
      nrow(locfdr_grid) - grid_size, " extra hyperparameter combinations",
      " to have a fully expanded grid. Using a grid size of ", nrow(locfdr_grid), "."
    ))
  }

  locfdr_grid_reduced <- reduce_locfdr_grid(
    test_statistics = test_statistics,
    locfdr_grid = locfdr_grid,
    lower_pi0 = lower_pi0,
    parallel_param = parallel_param,
    verbose = verbose
  )

  if (verbose & nrow(locfdr_grid_reduced) < nrow(locfdr_grid)) {
    message(paste0(
      nrow(locfdr_grid) - nrow(locfdr_grid_reduced), " hyperparameter combinations failed",
      ifelse(lower_pi0 > 0, paste(' or had pi0 below', lower_pi0), ''),
      " when run on the data. Using a grid size of ", nrow(locfdr_grid_reduced), "."
    ))
  }

  return(locfdr_grid_reduced)
}
