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

  # set up a dataframe with fdr, test statistics, and indices
  # sort by test statistics
  dat <- data.frame(
    'fdr' = fdr,
    'test_statistics' = test_statistics,
    'i' = 1:length(test_statistics)
  ) %>%
    dplyr::arrange(test_statistics)

  # cumulative mean and then choosing the last observation
  # of each fdr will get the average fdr among
  # t-statistics <= this t-statistic
  Fdrs <- dat %>%
    dplyr::mutate(
      'Fdr' = dplyr::cummean(fdr)
    ) %>%
    dplyr::group_by(fdr) %>%
    dplyr::summarize(
      'Fdr' = dplyr::last(Fdr)
    )

  # merge back with fdrs (may be multiple observations)
  # sort back in original order, and return Fdr
  Fdr <- merge(dat, Fdrs, by = 'fdr') %>%
    dplyr::arrange(i) %>%
    dplyr::pull(Fdr)

  return(Fdr)
}
