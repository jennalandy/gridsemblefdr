#' Check qvalue Row
#' @description For a given row of qvalue_grid, checks whether hyperparameter combinations
#' can be used for qvalue on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param qvalue_grid data.frame, each row is a possible set of hyperparameters for locfdr
#' @param lower_pi0 double, exclude hyperparameter combinations that give pi0 estimates below this threshold
#' @param row integer, row of qvalue_grid considered
#' @param df integer, degrees of freedom of test statistics, if known
#'
#' @return boolean, whether the hyperparameter combination ran without error
check_qvalue_row <- function(
  test_statistics,
  qvalue_grid,
  lower_pi0,
  row,
  df
) {
  # for each row, attempt to run qvalue
  run_i <- run_qvalue_row(
    test_statistics = test_statistics,
    qvalue_grid = qvalue_grid,
    row = row,
    df = df,
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

#' @title Reduce qvalue grid
#' @description Reduces grid to parameter combinations
#' that can run qvalue on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param qvalue_grid data.frame, each row is a possible set of hyperparameters for qvalue
#' @param lower_pi0 double, exclude hyperparameter combinations that give pi0 estimates below this threshold
#' @param df integer, degrees of freedom of test statistics, if known
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#' @param verbose boolean
#'
#' @importFrom magrittr %>%
#'
#' @return data.frame, each row is a possible set of hyperparameters for qvalue
reduce_qvalue_grid <- function(
  test_statistics,
  qvalue_grid,
  lower_pi0,
  df = NULL,
  parallel_param = NULL,
  verbose = FALSE
) {

  checks <- unlist(parlapply(
    X = 1:nrow(qvalue_grid),
    parallel_param = parallel_param,
    FUN = function(i, run_qvalue_row, qvalue_grid, test_statistics, df) {
      check = check_qvalue_row(
        test_statistics = test_statistics,
        qvalue_grid = qvalue_grid,
        row = i,
        df = df,
        lower_pi0 = lower_pi0
      )
      if (check == TRUE) {
        return(i)
      } else {
        return(check)
      }
    },
    run_qvalue_row = run_qvalue_row,
    qvalue_grid = qvalue_grid,
    test_statistics = test_statistics,
    df = df
  ))

  rows = grepl('^[[:digit:]]+$', checks)
  ok_rows = as.numeric(checks[rows])
  fails = checks[!rows]

  if(verbose) {
    message(paste0(
      '\t', length(ok_rows),'/',nrow(qvalue_grid),
      ' qvalue thetas included'
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

  return (qvalue_grid[ok_rows,])
}

#' @title Build qvalue Grid
#' @description Build a factorial grid of possible hyperparameters for qvalue and
#' reduce to hyperparameter sets that can be run on provided data without error.
#'
#' @param test_statistics vector, test statistics
#' @param transf vector, options for transf hyperparameter. `transf` is
#' a transformation is applied to the p-values so that a local FDR estimate can
#' be formed that does not involve edge effects of the \[0,1\] interval in which
#' the p-values lie, either "probit" or "logit".
#' @param adj_range vector c(min, max), range for adj hyperparameter. `adj` is a numeric value
#' that is applied as a multiple of the smoothing bandwidth used in the density
#' estimation.
#' @param pi0.method vector, options for pi0.method hyperparameter. `pi0.method`
#' is the method for automatically choosing tuning parameter in the estimation of
#' pi_0, the proportion of true null hypotheses, either "smoother" or "bootstrap".
#' @param smooth.log.pi0 vector of options for smooth.log.pi0 hyperparameter.
#' If `smooth.log.pi0` is TRUE and pi0.method = "smoother", pi_0 will be estimated
#' by applying a smoother to a scatterplot of log(pi_0) estimates against the
#' tuning parameter lambda.
#' @param df integer, degrees of freedom of test statistics, if known
#'
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
#' @return data.frame, each row is a possible set of hyperparameters for qvalue
#' @export
build_qvalue_grid <- function(
  test_statistics,
  transf = c('probit', 'logit'),
  adj_range = c(0.5, 1.5),
  pi0.method = c('bootstrap','smoother'),
  smooth.log.pi0 = c(TRUE, FALSE),
  df = NULL,
  grid_size = 40,
  lower_pi0 = 0.7,
  method = 'grid',
  seed = NULL,
  parallel_param = NULL,
  parallel = min(TRUE, n_workers > 1),
  n_workers = max(parallel::detectCores() - 2, 1),
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
    message('Building qvalue grid in parallel')
  } else if (verbose) {
    message('Building qvalue grid')
  }

  # smooth.log.pi0 is only used if pi0.method = 'smoother'
  # if smoother is a pi0.method option,
  # make sure that parameter combinations are not redundant
  # i.e. pi0.method = 'bootstrap', smooth.log.pi0 = T vs F doesn't matter
  non_smoother = pi0.method[pi0.method != 'smoother']
  prop_non_smoother = length(non_smoother)/length(pi0.method)
  n_non_smoother = round(grid_size*prop_non_smoother)
  n_smoother = grid_size - n_non_smoother

  # if random method, select from uniform distribution within ranges
  if (method == 'random') {

    # smooth.log.pi0 is only used if pi0.method = 'smoother'
    if ('smoother' %in% pi0.method) {

      qvalue_grid <- rbind(
        # n_smoother combinations
        # pi0.method = 'smoother'
        # smooth.log.pi0 has multiple options
        data.frame(
          transf = sample(transf, size = n_smoother, replace = TRUE),
          adj = stats::runif(n = n_smoother, min = adj_range[1], max = adj_range[2]),
          pi0.method = rep('smoother', n_smoother), #all smoother
          smooth.log.pi0 = sample(smooth.log.pi0, size = n_smoother, replace = TRUE)
        ),
        # n_non_smoother combinations
        # pi0.method sampled from non_smoother
        # smooth.log.pi0 has default value
        data.frame(
          transf = sample(transf, size = n_non_smoother, replace = TRUE),
          adj = stats::runif(n = n_non_smoother, min = adj_range[1], max = adj_range[2]),
          pi0.method = sample(non_smoother, size = n_non_smoother, replace = TRUE), # sample from non-smoother
          smooth.log.pi0 = rep(FALSE, n_non_smoother) #won't be used, keep as default
        )
      )
    } else {

      # if pi0.method = 'smoother' is not an option,
      # smooth.log.pi0 never matters, keep as default value
      qvalue_grid <- data.frame(
        transf = sample(transf, size = grid_size, replace = TRUE),
        adj = stats::runif(n = grid_size, min = adj_range[1], max = adj_range[2]),
        pi0.method = sample(pi0.method, size = grid_size, replace = TRUE),
        smooth.log.pi0 = rep(FALSE, grid_size) #won't be used, keep as default
      )
    }

    # if grid method, select equally spaced values within ranges
  } else if (method == 'grid') {

    # smooth.log.pi0 is only used if pi0.method = 'smoother'
    if ('smoother' %in% pi0.method) {
      non_smoother = pi0.method[pi0.method != 'smoother']
      prop_non_smoother = length(non_smoother)/length(pi0.method)
      n_non_smoother = round(grid_size*prop_non_smoother)
      n_smoother = grid_size - n_non_smoother

      n_combos_smoother = length(transf)*length(smooth.log.pi0)
      seq_length_smoother = ceiling(n_smoother/n_combos_smoother)

      n_combos_non_smoother = length(transf)*length(non_smoother)
      seq_length_non_smoother = ceiling(n_non_smoother/n_combos_non_smoother)

      qvalue_grid <- rbind(
        expand.grid(
          transf = transf,
          adj = seq(from = adj_range[1], to = adj_range[2], length.out = seq_length_smoother),
          pi0.method = "smoother",
          smooth.log.pi0 = smooth.log.pi0
        ),
        expand.grid(
          transf = transf,
          adj = seq(from = adj_range[1], to = adj_range[2], length.out = seq_length_non_smoother),
          pi0.method = non_smoother,
          smooth.log.pi0 = FALSE
        )
      )
    } else {
      n_combos = length(transf)*length(pi0.method)
      seq_length = ceiling(grid_size/n_combos)

      qvalue_grid <- data.frame(
        transf = transf,
        adj = seq(from = adj_range[1], to = adj_range[2], length.out = seq_length),
        pi0.method = pi0.method,
        smooth.log.pi0 = FALSE
      )
    }
  }

  if (verbose & nrow(qvalue_grid) > grid_size & method == 'grid') {
    message(paste0(
      "\t", nrow(qvalue_grid) - grid_size, " extra thetas",
      " needed for full grid"
    ))
  }

  qvalue_grid_reduced <- reduce_qvalue_grid(
    test_statistics = test_statistics,
    qvalue_grid = qvalue_grid,
    df = df,
    lower_pi0 = lower_pi0,
    parallel_param = parallel_param,
    verbose = verbose
  )

  return(qvalue_grid_reduced)
}
