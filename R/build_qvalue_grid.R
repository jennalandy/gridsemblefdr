#' @title Reduce qvalue grid
#' @description Reduces grid to parameter combinations
#' that can run qvalue on the provided data without error
#'
#' @param test_statistics vector of test statistics
#' @param qvalue_grid data frame where each row is a possible set of hyperparameters for qvalue
#' @param parallel if TRUE, processes are run in parallel
#' @param verbose if TRUE, status updates will be displayed
#'
#' @return dataframe where each row is a possible set of hyperparameters for qvalue
#'
#' @importFrom magrittr %>%
#' @export
reduce_qvalue_grid <- function(
  test_statistics,
  qvalue_grid,
  df = NULL,
  parallel = TRUE,
  verbose = FALSE
) {

  n.cores = ifelse(
    parallel,
    min(nrow(qvalue_grid), parallel::detectCores() - 1),
    1
  )
  cl = parallel::makeCluster(
    n.cores,
    type = 'PSOCK'
  )
  doParallel::registerDoParallel(cl)

  ok_rows <- foreach::foreach(
    i=1:nrow(qvalue_grid),
    .combine = c,
    .export = c(
      'run_qvalue_row'
    )
  ) %dopar% {
    # for each row, attempt to run qvalue
    run_i <- run_qvalue_row(
      test_statistics = test_statistics,
      qvalue_grid = qvalue_grid,
      row = i,
      df = df
    )

    # keep this row in final grid if it ran without error,
    # and if it does not estimate fdrs as all 0 or all 1
    if (!is.null(run_i)) {
      if (
        (sum(run_i$fdr != 1) > 0) &
        (sum(run_i$fdr != 0) > 0) &
        run_i$pi0 <= 1
      ) {
        return(i)
      }
    }
  }
  parallel::stopCluster(cl)

  if(verbose) {
    print(paste(
      length(ok_rows),'/',nrow(qvalue_grid),
      'qvalue parameter sets are usable'
    ))
  }

  return (qvalue_grid[ok_rows,])
}

#' @title Build qvalue Grid
#' @description Build a grid of possible hyperparameters for qvalue
#' from separate vectors of hyperparameter options. Final grid only considers
#' hyperparameters that can be run on provided data without error.
#'
#' @param test_statistics vector of test statistics
#' @param transf vector of options for transf hyperparameter. `transf` is
#' a transformation is applied to the p-values so that a local FDR estimate can
#' be formed that does not involve edge effects of the [0,1] interval in which
#' the p-values lie, either "probit" or "logit".
#' @param adj vector of options for adj hyperparameter. `adj` is a numeric value
#' that is applied as a multiple of the smoothing bandwidth used in the density
#' estimation.
#' @param pi0.method vector of options for pi0.method hyperparameter. `pi0.method`
#' is the method for automatically choosing tuning parameter in the estimation of
#' pi_0, the proportion of true null hypotheses, either "smoother" or "bootstrap".
#' @param smooth.log.pi0 vector of options for smooth.log.pi0 hyperparameter.
#' If `smooth.log.pi0` is TRUE and pi0.method = "smoother", pi_0 will be estimated
#' by applying a smoother to a scatterplot of log(pi_0) estimates against the
#' tuning parameter lambda.
#' @param parallel if TRUE, processes are run in parallel
#' @param verbose if TRUE, status updates will be displayed
#'
#' @return dataframe where each row is a possible set of hyperparameters for
#' qvalue on a specific set of test statistics
#' @export
build_qvalue_grid <- function(
  test_statistics,
  transf = c('probit', 'logit'),
  adj = c(0.8,0.9,1,1.1,1.2),
  pi0.method = c('bootstrap','smoother'),
  smooth.log.pi0 = c(TRUE, FALSE),
  df = NULL,
  parallel = TRUE,
  verbose = FALSE
) {

  qvalue_grid <- expand.grid(
    transf = transf,
    adj = adj,
    pi0.method = pi0.method,
    smooth.log.pi0 = smooth.log.pi0
  )

  # smooth.log.pi0 is only used if pi0.method = 'smoother'
  default_smooth.log.pi0 <- unique(smooth.log.pi0)[1]
  qvalue_grid <- qvalue_grid[
    qvalue_grid$pi0.method == 'smoother' |
      qvalue_grid$smooth.log.pi0 == default_smooth.log.pi0,
  ]

  qvalue_grid_reduced <- reduce_qvalue_grid(
    test_statistics = test_statistics,
    qvalue_grid = qvalue_grid,
    df = df,
    parallel = parallel,
    verbose = verbose
  )

  return(qvalue_grid_reduced)
}
