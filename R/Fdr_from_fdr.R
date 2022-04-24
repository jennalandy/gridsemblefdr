#' @title Fdr from fdr
#' @description Calculate tail-end false discovery rate (Fdr)
#' from local false discovery rate (fdr)
#'
#' @param fdr vector of local false discovery rate estimates
#' @param test_statistics vector of test statistics
#' @param direction direction of tail-end false discovery rate, c('left','right')
#'
#' @importFrom dplyr cummean last arrange mutate group_by summarize pull
#' @importFrom magrittr %>%
#'
#' @return vector of tail-end false discovery rates
#' @export
Fdr_from_fdr <- function(fdr, test_statistics, direction = 'left') {

  Fdr = vector(length = length(test_statistics))
  for (i in 1:length(test_statistics)) {
    Fdr[i] = mean(fdr[test_statistics <= test_statistics[i]])
  }

  return(Fdr)
}
