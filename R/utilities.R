#' general lapply function with specification of parallel or not
#'
#' @param X a vector (atomic or list) or an expression object.
#' @param FUN the function to be applied to each element of X
#' @param parallel_param BiocParallel object
#' @param ... additional parameters to `BiocParallel::bplapply` or `lapply`
#'
#' @return result of lapply
#'
#' @importFrom BiocParallel bplapply
#' @noRd
parlapply <- function(X, FUN, parallel_param, ...) {
  if (!is.null(parallel_param)) {
    BiocParallel::bplapply(X, FUN, BPPARAM = parallel_param, ...)
  } else {
    lapply(X, FUN, ...)
  }
}

#' @title Fdr from fdr
#' @description Calculate tail-end false discovery rate (Fdr)
#' from local false discovery rate (fdr)
#'
#' @param fdr vector, local false discovery rate estimates
#' @param test_statistics vector, test statistics
#'#'
#' @return vector, tail-end false discovery rates
#' @export
Fdr_from_fdr <- function(fdr, test_statistics) {

  Fdr = vector(length = length(test_statistics))
  for (i in seq_len(length(test_statistics))) {
    Fdr[i] = mean(fdr[
      abs(test_statistics) >= abs(test_statistics[i])
    ])
  }

  return(Fdr)
}

#' P-value from t-statistic
#' @description Get one-sided p-value from test statistic
#'
#' @param test_statistics vector, test statistics
#' @param df integer, degrees of freedom of test statistics t-distribution,
#' otherwise assumed standard normal
#'
#' @importFrom stats pnorm pt
#' @return vector of 2-sided p-values
#' @noRd
p_from_t <- function(test_statistics, df = NULL) {

  if (is.null(df)) {
    # assume standard normal
    one_sided <- unlist(lapply(test_statistics, function(z) {
      stats::pnorm(-1*abs(z))
    }))
  } else {
    # assume t_df
    one_sided <-unlist(lapply(test_statistics, function(z) {
      stats::pt(-1*abs(z), df = df)
    }))
  }

  return (2*one_sided)
}
