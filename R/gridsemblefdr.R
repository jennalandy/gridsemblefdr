#' Gridsemble
#' @description use ensemble methods to estimate local (fdr) and tail-end (Fdr)
#' false discovery rates.
#'
#' @param test_statistics vector, test statistics
#' @param df integer, degrees of freedom of test statistics t-distribution,
#' otherwise assumed standard normal
#'
#' @param locfdr_grid data.frame or 'default', rows are possible hyperparameters
#' for locfdr, or `build_locfdr_grid` will be run with default values, NULL to
#' exclude locfdr from gridsemble
#' @param fdrtool_grid data.frame or 'default', rows are possible
#' hyperparameters for fdrtool, or `build_fdrtool_grid` will be run with default
#' values, NULL to exclude fdrtool from gridsemble
#' @param qvalue_grid data.frame or 'default', rows are possible hyperparameters
#' for qvalue, or `build_qvalue_grid` will be run with default values, NULL to
#' exclude qvalue from gridsemble
#'
#' @param ensemble_size integer, number of models to ensemble
#' @param nsim integer, number of datasets to simulate and perform grid search
#' over, or 0 for models to be randomly selected.
#' @param lower_pi0 double, lower bound cutoff of `pi0` for model consideration,
#' efault 0.7 (assume majority null tests)
#' @param sim_size integer, size of simulated datasets
#'
#' @param n_workers integer, number of cores to use if parallel
#' @param parallel boolean, whether to utilize parallelization
#' @param parallel_param `BiocParallel` object
#' @param verbose boolean
#'
#' @importFrom BiocParallel DoparParam
#' @importFrom parallel detectCores
#'
#' @return
#' \itemize{
#'      \item `fdr`: vector, local false discovery rates
#'      \item `Fdr`: vector, left tail false discovery rates
#'      \item `pi0`: double, proportion of null test statistics
#'      \item `top_grid`: data.frame, hyperparameter sets in ensemble and
#'      their metrics on each simulated dataset
#'      \item `all_grids`: data.frame, all hyperparameter sets considered
#'      and their metrics on simulated data
#'      \item `generating_model`: list, parameters of generating model
#' }
#' @export
#' @examples
#' set.seed(123)
#' test_statistics = c(rnorm(800), runif(100, -10, -5), runif(100, 5, 10))
#' res = gridsemble(test_statistics)
#' res$pi0
#' res$fdr
gridsemble <- function(
  test_statistics,
  df = NULL,
  locfdr_grid = 'default',
  fdrtool_grid = 'default',
  qvalue_grid = 'default',
  ensemble_size = 10,
  nsim = 10,
  lower_pi0 = 0.7,
  sim_size = length(test_statistics),
  n_workers = max(parallel::detectCores() - 2, 1),
  parallel = min(TRUE, n_workers > 1),
  parallel_param = NULL,
  verbose = TRUE
) {

  if (typeof(locfdr_grid) == "character") {if (locfdr_grid == 'default') {
    locfdr_grid = build_locfdr_grid(
      test_statistics, lower_pi0 = lower_pi0, parallel = parallel,
      n_workers = n_workers, verbose = verbose
    )
  }}
  if (typeof(fdrtool_grid) == "character") {if (fdrtool_grid == 'default') {
    fdrtool_grid = build_fdrtool_grid(
      test_statistics, lower_pi0 = lower_pi0, parallel = parallel,
      n_workers = n_workers, verbose = verbose
    )
  }}
  if (typeof(qvalue_grid) == "character") {if (qvalue_grid == 'default') {
    qvalue_grid = build_qvalue_grid(
      test_statistics, df = df, lower_pi0 = lower_pi0, parallel = parallel,
      n_workers = n_workers, verbose = verbose
    )
  }}

  focus_metric = 'fdrerror'
  large_abs_metric = FALSE

  n_workers = min(n_workers, nsim)
  if (parallel & is.null(parallel_param) & n_workers > 1) {
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
    warning(paste0(
      "Grid size of ", grid_size, " is too small for ensemble size ",
      ensemble_size, ". Setting ensemble size to ", grid_size, ".", sep = ''
    ))
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
    default_fdrtool <- fdrtool::fdrtool(
      test_statistics, pct = 0, plot = 0, verbose = 0
    )
  }, error = function(e) {})

  tryCatch({
    default_qvalue <- qvalue::qvalue(p_values, plot = 0)
  }, error = function(e) {})

  if (is.null(default_qvalue)) {
    default_qvalue <- qvalue::qvalue(p_values, plot = 0, lambda = 0)
  }

  # fit generating model to test statistics
  if (nsim > 0) {
    generating_model <- fit_generating_model(test_statistics, verbose = verbose)
  } else {
    generating_model <- NULL
    if (verbose) {
      message('No generting model fit with nsim = 0')
    }
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
