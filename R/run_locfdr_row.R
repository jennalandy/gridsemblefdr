#' @title Run locfdr
#' @description Run locfdr with a specific set of parameters
#'
#' @param t vector of test statistics
#' @param locfdr_grid data frame where each row is a set of hyperparameters
#' @param row row of locfdr_grid, i.e. which set of hyperparameters to run locfdr with
#'
#' @return
#' \itemize{
#'   \item fdr - estimated local false discovery rates
#'   \item Fdr - estimated left tail-end false discovery rates
#'   \item pi0 - estimated proportion of tests that are null
#' }
#'
#' @importFrom locfdr locfdr
#' @importFrom dplyr case_when
run_locfdr_row <- function(t, locfdr_grid, row) {

  # consider pi0 estimation method that matches desired nulltype
  est_method <- dplyr::case_when(
    locfdr_grid[row,'nulltype'] == 0 ~ 'thest',
    locfdr_grid[row,'nulltype'] == 1 ~ 'mlest',
    TRUE ~ 'cmest' # nulltypes 2 and 3 both use central matcing
  )

  tryCatch(
    {
      # try to run locfdr with desired hyperparameters
      res <- locfdr::locfdr(
        t,
        pct = locfdr_grid$pct[row],
        pct0 = locfdr_grid$pct0[row],
        nulltype = locfdr_grid$nulltype[row],
        type = locfdr_grid$type[row],
        plot = 0
      )

      return(list(
        'fdr' = res$fdr,
        'Fdr' = Fdr_from_fdr(res$fdr, t),
        'pi0' = res$fp0[est_method,'p0']
      ))
    },
    warning = function(w) {
      if (
        # ignore if one of these three warnings, run again
        grepl('misfit', w) |
        grepl('Discrepancy between central matching', w) |

        # if nulltype = 3 (i.e. using cm) this will be an error
        # as a warning it's just for the pi0 estimation
        grepl('CM estimation failed', w)
      ) {
        res <- suppressWarnings(locfdr::locfdr(
          t,
          pct = locfdr_grid$pct[row],
          pct0 = locfdr_grid$pct0[row],
          nulltype = locfdr_grid$nulltype[row],
          type = locfdr_grid$type[row],
          plot = 0
        ))

        return(list(
          'fdr' = res$fdr,
          'Fdr' = Fdr_from_fdr(res$fdr, t),
          'pi0' = res$fp0[est_method,'p0']
        ))
      } else {
        # for other warnings, don't consider this set of hyperparameters
        return(NULL)
      }
    },
    error = function(e) {
      # for errors, don't consider this set of hyperparameters
      return(NULL)
    }
  )
}
