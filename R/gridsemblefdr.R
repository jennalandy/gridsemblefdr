#' Gridsemble
#' @description use ensemble methods to estimate local (fdr) and tail-end (Fdr)
#' false discovery rates.
#'
#' @param test_statistics vector of test statistics
#' @param nsim Number of datasets to simulate and grid search over.
#' If 0, no datasets are simulated and model(s) are randomly selected.
#' @param topn Number of models chosen from each simulation to ensemble over.
#' If `nsim = 0`, number of random models to ensemble over.
#' @param locfdr_grid data frame where each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data frame where each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data frame where each row is a set of hyperparameters for qvalue
#' @param df degrees of freedom of test statistics, if known. Otherwise assumed
#' to be from N(0, 1).
#' @param methods vector of methods to be used. Default c('locfdr','fdrtool','qvalue')
#' @param asym boolean, whether to consider an asymmetric simulation model.
#' @param focus_metric which metric to prioritize in the grid search.
#' Must be one of c('pr','roc','brier','Fdrerror')
#' @param large_abs_metric if TRUE, only consider focus_metric looking at the
#' large absolute value test statistics (specifically, top quartile of abs(t))
#' @param parallel if TRUE, process is run in parallel
#' @param verbose
#'
#' @return
#' \itemize{
#'      \item fdr local false discovery rates
#'      \item Fdr left tail false discovery rates
#'      \item pi0 proportion of null test statistics
#'      \item top_grid dataframe containing hyperparameter sets
#'      that were ensembled over and their metrics on simulated daa
#'      \item all_grids dataframe containing all hyperparameter sets considered
#'      and their metrics on simulated data
#'      \item fit parameters used for simulation step
#' }
#' @export
gridsemble <- function(
  test_statistics,
  nsim = 10,
  topn = 1,
  locfdr_grid = NULL,
  fdrtool_grid = NULL,
  qvalue_grid = NULL,
  df = NULL,
  methods = c('locfdr','fdrtool','qvalue'),
  asym = FALSE,
  focus_metric = 'pr',
  large_abs_metric = TRUE,
  parallel = TRUE,
  verbose = TRUE
) {

  if (!(focus_metric %in% c('Fdrerror','roc','pr','brier'))) {
    cat("focus_metric must be one of c('Fdrerror','roc','pr','brier'), using default (Fdrerror)")
    focus_metric = 'Fdrerror'
  }

  method_list = c()
  row_list = c()
  if ('locfdr' %in% methods) {
    if (is.null(locfdr_grid)) {
      locfdr_grid <- build_locfdr_grid(test_statistics, parallel = parallel)
    }
    if (nrow_null0(locfdr_grid) > 0) {
      method_list = c(method_list, rep('locfdr', nrow(locfdr_grid)))
      row_list = c(row_list, 1:nrow(locfdr_grid))
    }
  }
  if ('fdrtool' %in% methods) {
    if (is.null(fdrtool_grid)) {
      fdrtool_grid <- build_fdrtool_grid(test_statistics, parallel = parallel)
    }
    if (nrow_null0(fdrtool_grid) > 0) {
      method_list = c(method_list, rep('fdrtool', nrow(fdrtool_grid)))
      row_list = c(row_list, 1:nrow(fdrtool_grid))
    }
  }
  if ('qvalue' %in% methods) {
    if (is.null(qvalue_grid)) {
      qvalue_grid <- build_qvalue_grid(test_statistics, df = df, parallel = parallel)
    }
    if (nrow_null0(qvalue_grid) > 0) {
      method_list = c(method_list, rep('qvalue', nrow(qvalue_grid)))
      row_list = c(row_list, 1:nrow(qvalue_grid))
    }
  }

  default_locfdr <- locfdr::locfdr(test_statistics, plot = 0)
  default_fdrtool <- fdrtool::fdrtool(test_statistics, pct = 0, plot = 0, verbose = 0)
  p_values = p_from_t(
    test_statistics = test_statistics,
    df = df,
    sides = 'two'
  )
  default_qvalue <- NULL
  tryCatch({
    default_qvalue <- qvalue::qvalue(p_values, plot = 0)
  }, error = function(e) {
  })
  if (is.null(default_qvalue)) {
    default_qvalue <- qvalue::qvalue(p_values, plot = 0, lambda = 0)
  }

  if (nsim == 0) {

    # just ensemble, no simulation / grid search
    top_grid = data.frame(
      'method' = method_list,
      'row' = row_list
    )

    top_grid = top_grid[
      sample(1:nrow(top_grid), size = topn),
    ]

    ensemble_res = ensemble(
      test_statistics = test_statistics,
      top_grid = top_grid,
      locfdr_grid = locfdr_grid,
      fdrtool_grid = fdrtool_grid,
      qvalue_grid = qvalue_grid,
      df = df,
      verbose = verbose
    )

    fit = NULL
    all_grids = NULL

  } else {
    asym_type = ifelse(asym, 'asymmetric', 'symmetric')

    fit <- fit_sim(test_statistics, type = asym_type)
    fit$parameters$pi0 = 0.8 # todo delete this

    grid_res <- grid_search(
      n = length(test_statistics),
      nsim = nsim,
      topn = topn,
      fit = fit,
      method_list = method_list,
      row_list = row_list,
      df = df,
      locfdr_grid = locfdr_grid,
      qvalue_grid = qvalue_grid,
      fdrtool_grid = fdrtool_grid,
      focus_metric = focus_metric,
      large_abs_metric = large_abs_metric,
      params_type = asym_type,
      parallel = parallel,
      verbose = verbose
    )

    top_grid = grid_res$top_grid
    all_grids = grid_res$all_grids
  }

  ensemble_res = ensemble(
    test_statistics = test_statistics,
    top_grid = top_grid,
    locfdr_grid = locfdr_grid,
    fdrtool_grid = fdrtool_grid,
    qvalue_grid = qvalue_grid,
    df = df,
    verbose = verbose
  )

  Fdr <- Fdr_from_fdr(
    fdr = ensemble_res$fdr,
    test_statistics = test_statistics
  )

  out <- list(
    'fdr' = ensemble_res$fdr,
    'Fdr' = Fdr,
    'pi0' = ensemble_res$pi0,
    'top_grid' = top_grid,
    'default_locfdr' = default_locfdr,
    'default_fdrtool' = default_fdrtool,
    'default_qvalue' = default_qvalue,
    'all_grids' = all_grids,
    'fit' = fit
  )

  return(out)
}