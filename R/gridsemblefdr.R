#' Gridsemble
#' @description use ensemble methods to estimate local (fdr) and tail-end (Fdr)
#' false discovery rates.
#'
#' @param test_statistics vector, test statistics
#' @param nsim integer, number of datasets to simulate and grid search over.
#' If 0, no datasets are simulated and model(s) are randomly selected.
#' @param ensemble_size integer, number of models chosen to ensemble over.
#' @param locfdr_grid data.frame, each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data.frame, each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data.frame, each row is a set of hyperparameters for qvalue
#' @param df integer, degrees of freedom of test statistics, if known. Otherwise assumed
#' to be from N(0, 1).
#' @param methods vector, methods to be used. Default c('locfdr','fdrtool','qvalue')
#' @param focus_metric string, one of one of c('fdrerror','Fdrerror','pr','roc','brier'),
#' which metric to optimize in the grid search
#' @param large_abs_metric boolean, if TRUE, only consider focus_metric looking at the
#' large absolute value test statistics (top quartile of abs(t))
#' @param parallel boolean
#' @param paralllel_param BiocParallel object, specified to run in parallel or NULL
#' @param sim_n integer, size of simulated datasets, default number of test statistics
#' @param verbose boolean
#'
#' @importFrom BiocParallel DoparParam
#' @importFrom parallel detectCores
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
  ensemble_size = 10,
  locfdr_grid = NULL,
  fdrtool_grid = NULL,
  qvalue_grid = NULL,
  df = NULL,
  methods = c('locfdr','fdrtool','qvalue'),
  focus_metric = 'fdrerror',
  lower_pi0 = 0.7,
  n_workers = NULL,
  large_abs_metric = TRUE,
  parallel = TRUE,
  parallel_param = NULL,
  sim_n = NULL,
  verbose = TRUE
) {

  if (parallel & is.null(parallel_param)) {
    if (is.null(n_workers)) {
      n_workers = parallel::detectCores() - 2
    }

    parallel_param <- BiocParallel::MulticoreParam(
      workers = n_workers,
      tasks = n_workers
    )
  }
  if (!parallel) {
    parallel_param <- NULL
  }

  p_values = p_from_t(
    test_statistics = test_statistics,
    df = df,
    sides = 'two'
  )

  if (!(focus_metric %in% c('Fdrerror','fdrerror','roc','pr','brier','accuracy','recall','precision','specificity','f1'))) {
    cat("focus_metric must be one of c('Fdrerror','fdrerror','roc','pr','brier'), using default (fdrerror)")
    focus_metric = 'fdrerror'
  }

  if ('locfdr' %in% methods) {
    if (is.null(locfdr_grid)) {
      locfdr_grid <- build_locfdr_grid(test_statistics, lower_pi0=lower_pi0, parallel_param = parallel_param)
    }
  }
  if ('fdrtool' %in% methods) {
    if (is.null(fdrtool_grid)) {
      fdrtool_grid <- build_fdrtool_grid(test_statistics, lower_pi0=lower_pi0, parallel_param = parallel_param)
    }
  }
  if ('qvalue' %in% methods) {
    if (is.null(qvalue_grid)) {
      qvalue_grid <- build_qvalue_grid(test_statistics, df = df, lower_pi0=lower_pi0, parallel_param = parallel_param)
    }
  }

  grid_size = nrow_null0(locfdr_grid) +
    nrow_null0(fdrtool_grid) +
    nrow_null0(qvalue_grid)

  if (grid_size < ensemble_size) {
    ensemble_size = grid_size
    warning(paste("Grid size of ", grid_size, " is too small for ensemble size ", ensemble_size,
                  ". Setting ensemble size to ", grid_size, ".", sep = ''))
  }

  # store results from each implementation with default parameters
  # leave null if unable to run implementation

  default_locfdr <- NULL
  default_fdrtool <- NULL
  default_qvalue <- NULL

  tryCatch({
    default_locfdr <- locfdr::locfdr(test_statistics, plot = 0)
  }, error = function(e) {})

  tryCatch({
    default_fdrtool <- fdrtool::fdrtool(test_statistics, pct = 0, plot = 0, verbose = 0)
  }, error = function(e) {})

  tryCatch({
    default_qvalue <- qvalue::qvalue(p_values, plot = 0)
  }, error = function(e) {})

  if (is.null(default_qvalue)) {
    default_qvalue <- qvalue::qvalue(p_values, plot = 0, lambda = 0)
  }

  # fit simulation model to test statistics
  if (nsim > 0) {
    fit <- fit_sim(test_statistics, type = 'symmetric')
  } else {
    fit <- NULL
  }

  # perform grid search on simulated datasets
  grid_res <- grid_search(
    n = length(test_statistics),
    nsim = nsim,
    ensemble_size = ensemble_size,
    fit = fit,
    df = df,
    locfdr_grid = locfdr_grid,
    qvalue_grid = qvalue_grid,
    fdrtool_grid = fdrtool_grid,
    focus_metric = focus_metric,
    large_abs_metric = large_abs_metric,
    params_type = 'symmetric',
    parallel_param = parallel_param,
    sim_n = sim_n,
    verbose = verbose
  )

  top_grid = grid_res$top_grid
  all_grids = grid_res$all_grids

  # ensemble over top performing grid search methods
  ensemble_res = ensemble(
    test_statistics = test_statistics,
    focus_metric = focus_metric,
    large_abs_metric = large_abs_metric,
    top_grid = top_grid,
    locfdr_grid = locfdr_grid,
    fdrtool_grid = fdrtool_grid,
    qvalue_grid = qvalue_grid,
    df = df,
    parallel_param = parallel_param,
    verbose = verbose
  )

  # compute tail end Fdr from local fdr
  Fdr <- Fdr_from_fdr(
    fdr = ensemble_res$fdr,
    test_statistics = test_statistics
  )

  out <- list(
    'fdr' = ensemble_res$fdr,
    'fdr_var' = ensemble_res$fdr_var,
    'Fdr' = Fdr,
    'pi0' = ensemble_res$pi0,
    'pi0_var' = ensemble_res$pi0_var,
    'top_grid' = top_grid,
    'default_locfdr' = default_locfdr,
    'default_fdrtool' = default_fdrtool,
    'default_qvalue' = default_qvalue,
    'all_grids' = all_grids,
    'fit' = fit,
    'locfdr_grid' = locfdr_grid,
    'fdrtool_grid' = fdrtool_grid,
    'qvalue_grid' = qvalue_grid
  )
  return(out)
}
