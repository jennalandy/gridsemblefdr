#' Gridsemble
#' @description use ensemble methods to estimate local (fdr) and tail-end (Fdr)
#' false discovery rates.
#'
#' @param test_statistics vector, test statistics
#' @param df integer, degrees of freedom of test statistics, if known. Otherwise assumed
#' to be from N(0, 1).
#'
#' @param locfdr_grid data.frame, each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data.frame, each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data.frame, each row is a set of hyperparameters for qvalue
#'
#' @param nsim integer, number of datasets to simulate and grid search over.
#' If 0, no datasets are simulated and model(s) are randomly selected.
#' @param ensemble_size integer, number of models chosen to ensemble over.
#' @param lower_pi0 double, lower bound cutoff of pi0 for model consideration.
#' Default 0.7 (assume majority null tests).
#' @param sim_size integer, size of simulated datasets, default number of test statistics
#'
#' @param parallel boolean, whether to utilize parallelization
#' @param n_workers integer, number of cores to use if parallel, default 2 less than available.
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
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
#'      \item generating_model parameters used for simulation step
#' }
#' @export
gridsemble <- function(
  test_statistics,
  df = NULL,
  locfdr_grid = build_locfdr_grid(test_statistics, lower_pi0 = lower_pi0, parallel = parallel,n_workers = n_workers, verbose = verbose),
  fdrtool_grid = build_fdrtool_grid(test_statistics, lower_pi0 = lower_pi0, parallel = parallel, n_workers = n_workers, verbose = verbose),
  qvalue_grid = build_qvalue_grid(test_statistics, df = df, lower_pi0 = lower_pi0, parallel = parallel, n_workers = n_workers, verbose = verbose),
  nsim = 10,
  ensemble_size = 10,
  lower_pi0 = 0.7,
  sim_size = length(test_statistics),
  parallel = TRUE,
  n_workers = parallel::detectCores() - 2,
  parallel_param = NULL,
  verbose = TRUE
) {

  focus_metric = 'fdrerror'
  large_abs_metric = FALSE

  if (parallel & is.null(parallel_param)) {
    parallel_param = BiocParallel::MulticoreParam(
      workers = n_workers,
      tasks = n_workers
    )
  }

  p_values = p_from_t(
    test_statistics = test_statistics,
    df = df
  )

  methods = c()
  if (!is.null(locfdr_grid)) {
    methods = c(methods, 'locfdr')
  }
  if (is.null(fdrtool_grid)) {
    methods = c(methods, 'fdrtool')
  }
  if (is.null(qvalue_grid)) {
    methods = c(methods, 'qvalue')
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

  # fit generating model to test statistics
  if (nsim > 0) {
    generating_model <- fit_generating_model(test_statistics)
  } else {
    generating_model <- NULL
  }

  # perform grid search on simulated datasets
  grid_res <- grid_search(
    generating_model = generating_model,
    nsim = nsim,
    sim_size = sim_size,
    ensemble_size = ensemble_size,
    df = df,
    locfdr_grid = locfdr_grid,
    qvalue_grid = qvalue_grid,
    fdrtool_grid = fdrtool_grid,
    focus_metric = focus_metric,
    large_abs_metric = large_abs_metric,
    parallel = parallel,
    n_workers =  n_workers,
    parallel_param = parallel_param,
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
    'generating_model' = generating_model,
    'locfdr_grid' = locfdr_grid,
    'fdrtool_grid' = fdrtool_grid,
    'qvalue_grid' = qvalue_grid
  )
  return(out)
}
