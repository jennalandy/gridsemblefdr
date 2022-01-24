#' @title Fdr from fdr
#' @description Calculate tail-end false discovery rate (Fdr)
#' from local false discovery rate (fdr)
#'
#' @param fdr vector of local false discovery rate estimates
#' @param t vector of test statistics
#' @param direction direction of tail-end false discovery rate, c('left','right')
#'
#' @return vector of tail-end false discovery rates
#'
#' @examples
#'
Fdr_from_fdr <- function(fdr, t, direction = 'left') {

  sapply(t, function(x) {
    if (direction == 'left') {
      selection  <- (t <= x)
    } else if (direction == 'right') {
      selection  <- (t >= x)
    } else {
      return(NULL)
    }
    mean(fdr[selection])
  })
}
