#' general lapply function with specification of parallel or not
#'
#' @param X a vector (atomic or list) or an expression object.
#' @param FUN the function to be applied to each element of X
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
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
