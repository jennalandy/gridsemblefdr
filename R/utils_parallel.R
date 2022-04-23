#' Title
#'
#' @param X a vector (atomic or list) or an expression object.
#' @param FUN the function to be applied to each element of X
#' @param parallel whether to run in parallel with `bplapply` or regularly with `lapply`
#'
#' @return result of lapply
#'
#' @importFrom BiocParallel bplapply
#'
#' @export
parlapply <- function(X, FUN, parallel_param, ...) {
  if (!is.null(parallel_param)) {
    BiocParallel::bplapply(X, FUN, BPPARAM = parallel_param, ...)
  } else {
    lapply(X, FUN, ...)
  }
}
