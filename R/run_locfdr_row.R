#' @title Run locfdr
#' @description Run locfdr with a specific set of parameters
#'
#' @param test_statistics vector, test statistics
#' @param locfdr_grid data.frame, each row is a set of hyperparameters
#' @param row integer, row of locfdr_grid
#' @param returnFdr boolean, whether to calculate Fdr form fdr
#'
#' @return
#' list of estimates
#' \itemize{
#'   \item `fdr`: vector, local false discovery rates
#'   \item `Fdr`: vector, tail-end false discovery rates if returnFdr = TRUE
#'   \item `pi0`: double, proportion of tests that are null
#' }
#'
#' @importFrom locfdr locfdr
#' @importFrom dplyr case_when
#' @noRd
run_locfdr_row <- function(
  test_statistics,
  locfdr_grid,
  row,
  returnFdr = TRUE
) {
  # pi0 estimation method to matches desired nulltype
  est_method <- dplyr::case_when(
    locfdr_grid[row,'nulltype'] == 0 ~ 'thest',
    locfdr_grid[row,'nulltype'] == 1 ~ 'mlest',
    TRUE ~ 'cmest' # nulltypes 2 and 3 both use central matcing
  )

  tryCatch(
    {
      res <- locfdr::locfdr(
        test_statistics,
        pct = locfdr_grid$pct[row], pct0 = locfdr_grid$pct0[row],
        nulltype = locfdr_grid$nulltype[row], type = locfdr_grid$type[row],
        plot = 0
      )

      if (returnFdr) {
        Fdr = Fdr_from_fdr(fdr = res$fdr,test_statistics = test_statistics)
      } else { Fdr = 'not computed' }

      return(list(
        'fdr' = res$fdr, 'Fdr' = Fdr, 'pi0' = res$fp0[est_method,'p0']
      ))
    },
    warning = function(w) {
      if (
        grepl('misfit', w) | grepl('Discrepancy between central matching', w) |
        grepl('CM estimation failed', w)
      ) {
        res <- suppressWarnings(locfdr::locfdr(
          test_statistics,
          pct = locfdr_grid$pct[row], pct0 = locfdr_grid$pct0[row],
          nulltype = locfdr_grid$nulltype[row], type = locfdr_grid$type[row],
          plot = 0
        ))

        if (returnFdr) {
          Fdr = Fdr_from_fdr(fdr = res$fdr,test_statistics = test_statistics)
        } else { Fdr = 'not computed' }

        return(list(
          'fdr' = res$fdr, 'Fdr' = Fdr, 'pi0' = unlist(res$fp0[est_method,'p0'])
        ))
      } else { return(NULL) }
    },
    error = function(e) { return(NULL) }
  )
}

# Note: 'CM estimation failed' is only as a warning for pi0 estimations
# if CM fails when the type/nulltype uses central matching, it will be an error
