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
