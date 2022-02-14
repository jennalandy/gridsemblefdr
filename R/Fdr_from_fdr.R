#' @title Fdr from fdr
#' @description Calculate tail-end false discovery rate (Fdr)
#' from local false discovery rate (fdr)
#'
#' @param fdr vector of local false discovery rate estimates
#' @param test_statistics vector of test statistics
#' @param direction direction of tail-end false discovery rate, c('left','right')
#'
#' @return vector of tail-end false discovery rates
Fdr_from_fdr <- function(fdr, test_statistics, direction = 'left') {

  sapply(test_statistics, function(t) {
    if (direction == 'left') {
      selection  <- (test_statistics <= t)
    } else if (direction == 'right') {
      selection  <- (test_statistics >= t)
    } else {
      return(NULL)
    }
    mean(fdr[selection])
  })
}
