#' Ensemble
#'
#' @param test_statistics vector of test statistics
#' @param top_grid data frame of method, row number combinations to ensemble
#' @param locfdr_grid data frame where each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data frame where each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data frame where each row is a set of hyperparameters for qvalue
#' @param df degrees of freedom of test statistics, if known
#' @param verbose
#'
#' @return
#' \itemize{
#'   \item fdr - estimated local false discovery rates
#'   \item pi0 - estimated proportion of tests that are null
#' }
#' @export
ensemble <- function(
  test_statistics,
  top_grid,
  locfdr_grid,
  fdrtool_grid,
  qvalue_grid,
  df = NULL,
  parallel = TRUE,
  verbose = TRUE
) {

  n.cores = ifelse(
    parallel,
    min(nrow(top_grid), parallel::detectCores() - 1),
    1
  )
  cl = parallel::makeCluster(
    n.cores,
    type = 'PSOCK'
  )
  doParallel::registerDoParallel(cl)

  # iterate through models to ensemble over
  vec_mean = function(vec1, vec2) {
    (vec1 + vec2)/2
  }

  pi0_and_fdr <- foreach (
    i = 1:nrow(top_grid),
    .combine = vec_mean,
    .export = c(
      "run_locfdr_row",
      "run_qvalue_row",
      "run_fdrtool_row"
    )
  ) %dopar% {
    method = top_grid[i,'method']
    if (method == 'locfdr') {
      i_fdr <- run_locfdr_row(
        test_statistics = test_statistics,
        locfdr_grid = locfdr_grid,
        row = as.numeric(top_grid[i, 'row'])
      )
    } else if (method == 'fdrtool') {
      i_fdr <- run_fdrtool_row(
        test_statistics = test_statistics,
        fdrtool_grid = fdrtool_grid,
        row = as.numeric(top_grid[i, 'row'])
      )
    } else if (method == 'qvalue') {
      i_fdr <- run_qvalue_row(
        test_statistics = test_statistics,
        qvalue_grid = qvalue_grid,
        row = as.numeric(top_grid[i, 'row']),
        df = df
      )
    }

    if (!is.null(i_fdr)) {
      if (
          (sum(i_fdr$fdr != 1) > 0) &
          (sum(i_fdr$fdr != 0) > 0)
      ) {
        return(c(unname(i_fdr$pi0), i_fdr$fdr))
      }
    }
  }

  parallel::stopCluster(cl)

  topn_fdr <- pi0_and_fdr[2:length(pi0_and_fdr)]
  pi0_est <- pi0_and_fdr[1]

  return(list(
    fdr = topn_fdr,
    pi0 = pi0_est
  ))
}
