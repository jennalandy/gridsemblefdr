#' Check qvalue Row
#' @description For a given row of locfdr_grid, checks whether hyperparameter
#' combinations can be used for locfdr on the provided data without error
#'
#' @param test_statistics vector, test statistics
#' @param qvalue_grid data.frame, rows are possible hyperparameters for locfdr
#' @param row integer, row of qvalue_grid considered
#' @param df integer, degrees of freedom of test statistics, if known
#'
#' @return boolean, whether the hyperparameter combination ran without error
#' @noRd
check_qvalue_row <- function(
  test_statistics,
  qvalue_grid,
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
    }  else {
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
#' @param qvalue_grid data.frame, rows are possible hyperparameters for qvalue
#' @param df integer, degrees of freedom of test statistics, if known
#' @param parallel_param BiocParallel object
#' @param verbose boolean
#'
#' @importFrom magrittr %>%
#'
#' @return data.frame, each row is a possible set of hyperparameters for qvalue
#' @noRd
reduce_qvalue_grid <- function(
  test_statistics,
  qvalue_grid,
  df = NULL,
  parallel_param = NULL,
  verbose = FALSE
) {

  checks <- unlist(parlapply(
    X = seq_len(nrow(qvalue_grid)),
    parallel_param = parallel_param,
    FUN = function(i, run_qvalue_row, qvalue_grid, test_statistics, df) {
      check = check_qvalue_row(
        test_statistics = test_statistics,
        qvalue_grid = qvalue_grid,
        row = i,
        df = df
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
      ' qvalue models included'
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

  return (qvalue_grid[ok_rows,])
}

#' @title Build qvalue Grid
#' @description Build a factorial grid of hyperparameters for qvalue and
#' reduce to hyperparameter sets that can be run on provided data without error.
#'
#' @param test_statistics vector, test statistics
#' @param transf vector, options for transf hyperparameter. `transf` is
#' a transformation is applied to the p-values so that a local FDR estimate can
#' be formed that does not involve edge effects of the \[0,1\] interval in which
#' the p-values lie, either "probit" or "logit".
#' @param adj_range vector c(min, max), range for adj hyperparameter. `adj` is a
#' numeric value that is applied as a multiple of the smoothing bandwidth used
#' in the density estimation.
#' @param pi0.method vector, options for pi0.method hyperparameter. `pi0.method`
#' is the method for automatically choosing tuning parameter in the estimation
#' of pi_0, the proportion of true null tests, one of c("smoother","bootstrap").
#' @param smooth.log.pi0 vector of options for smooth.log.pi0 hyperparameter.
#' If `smooth.log.pi0` is TRUE and pi0.method = "smoother", pi_0 will be e
#' stimated by applying a smoother to a scatterplot of log(pi_0) estimates
#' against the tuning parameter lambda.
#' @param df integer, degrees of freedom of test statistics, if known
#'
#' @param grid_depth integer, number of evenly-spaced values of continuous parameters
#' considered within their respective `_range`.
#' @param method string, one of c('random', 'grid'). 'random' will sample
#' uniformly within the ranges, 'grid' will select equally spaced values
#'
#' @param parallel boolean, whether to utilize parallelization
#' @param n_workers integer, number of cores to use if parallel
#' @param parallel_param BiocParallel object
#' @param verbose boolean
#'
#' @importFrom stats runif
#' @return data.frame, each row is a possible set of hyperparameters for qvalue
#' @export
#'
#' @examples
#' set.seed(123)
#' test_statistics = c(rnorm(800), runif(100, -10, -5), runif(100, 5, 10))
#' qvalue_grid = build_qvalue_grid(test_statistics)
build_qvalue_grid <- function(
  test_statistics, transf = c('probit', 'logit'), adj_range = c(0.5, 2),
  pi0.method = c('bootstrap','smoother'), smooth.log.pi0 = c(TRUE, FALSE),
  df = NULL, grid_depth = 20,
  parallel_param = NULL, parallel = min(TRUE, n_workers > 1),
  n_workers = max(parallel::detectCores() - 2, 1), verbose = FALSE
) {

  adj = seq(adj_range[1], adj_range[2], length.out = grid_depth)

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

  non = pi0.method[pi0.method != 'smoother']
  qvalue_grid <- expand.grid(
    transf = transf,
    adj = adj,
    pi0.method = non,
    smooth.log.pi0 = FALSE # only affects model if pi0.method = "smoother"
  )

  if ('smoother' %in% pi0.method) {
    qvalue_grid <- rbind(
      qvalue_grid,
      expand.grid(
        transf = transf,
        adj = adj,
        pi0.method = "smoother",
        smooth.log.pi0 = smooth.log.pi0
      )
    )
  }

  qvalue_grid_reduced <- reduce_qvalue_grid(
    test_statistics = test_statistics, qvalue_grid = qvalue_grid, df = df,
    parallel_param = parallel_param, verbose = verbose
  )

  return(qvalue_grid_reduced)
}
