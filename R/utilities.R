#' general lapply function with specification of parallel or not
#'
#' @param X a vector (atomic or list) or an expression object.
#' @param FUN the function to be applied to each element of X
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#' @param ... additional parameters to `BiocParallel::bplapply` or `lapply`
#'
#' @return result of lapply
#'
#' @importFrom BiocParallel bplapply
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
#' @param direction string, one of c('left','right'), direction of tail-end false discovery rate
#'
#' @importFrom dplyr cummean last arrange mutate group_by summarize pull
#' @importFrom magrittr %>%
#'
#' @return vector, tail-end false discovery rates
Fdr_from_fdr <- function(fdr, test_statistics, direction = 'left') {

  Fdr = vector(length = length(test_statistics))
  for (i in 1:length(test_statistics)) {
    Fdr[i] = mean(fdr[test_statistics <= test_statistics[i]])
  }

  return(Fdr)
}

#' P-value from t-statistic
#' @description Get one-sided p-value from test statistic
#'
#' @param test_statistics vector, test statistics
#' @param df integer, degrees of freedom to compute p-value from test statistics,
#' assume standard normal if NULL
#'
#' @importFrom stats pnorm pt
#' @return vector of 2-sided p-values
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
