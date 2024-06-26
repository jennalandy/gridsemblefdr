#' Gridsemble
#' @description use ensemble methods to estimate local (fdr) and tail-end (Fdr)
#' false discovery rates.
#'
#' @param test_statistics vector, test statistics
#' @param df integer, degrees of freedom of test statistics null t-distribution,
#' otherwise assumed standard normal. Ignored if `to_pval_function` is provided.
#' @param to_pval_function function, converts test statistics vector to a
#' p-value vector. Default assumes t-distribution with given df under the null.
#' @param ensemble_size integer, number of models to ensemble
#' @param n_synthetic integer, number of datasets to simulate and perform grid search
#' over, or 0 for models to be randomly selected.
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
#' @param synthetic_size integer, size of synthetic datasets. Defaults to
#' `length(test_statistics)`. Can be reduced for time/memory constraints.
#'
#' @param n_workers integer, number of cores to use if parallel
#' @param parallel boolean, whether to utilize parallelization
#' @param verbose boolean, whether to report status
#'
#' @importFrom BiocParallel DoparParam
#' @importFrom parallel detectCores
#'
#' @return
#' \itemize{
#'      \item `fdr`: vector, estimated local false discovery rates
#'      \item `Fdr`: vector, estimated left tail-end false discovery rates
#'      \item `pi0`: double, estimated proportion of null test statistics
#'      \item `top_grid`: data.frame, hyperparameter sets in ensemble and
#'      their estimated metrics on each simulated dataset
#'      \item `all_grids`: data.frame, all hyperparameter sets considered
#'      and their estimated metrics on simulated data
#'      \item `working_model`: list, parameters of working model
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
  to_pval_function = function(test_statistics) {p_from_t(test_statistics, df = df)},
  ensemble_size = 10,
  n_synthetic = 10,
  locfdr_grid = 'default',
  fdrtool_grid = 'default',
  qvalue_grid = 'default',
  synthetic_size = length(test_statistics),
  n_workers = max(parallel::detectCores() - 2, 1),
  parallel = min(TRUE, n_workers > 1),
  verbose = TRUE,
  drop_pi0_1 = TRUE
) {

  focus_metric = 'fdrerror'
  p_values = to_pval_function(
    test_statistics = test_statistics
  )

  if (typeof(locfdr_grid) == "character") {if (locfdr_grid == 'default') {
    locfdr_grid = build_locfdr_grid(
      test_statistics, drop_pi0_1 = drop_pi0_1, parallel = parallel,
      n_workers = n_workers, verbose = verbose
    )
  }}
  if (typeof(fdrtool_grid) == "character") {if (fdrtool_grid == 'default') {
    fdrtool_grid = build_fdrtool_grid(
      test_statistics, drop_pi0_1 = drop_pi0_1, parallel = parallel,
      n_workers = n_workers, verbose = verbose
    )
  }}
  if (typeof(qvalue_grid) == "character") {if (qvalue_grid == 'default') {
    qvalue_grid = build_qvalue_grid(
      test_statistics, drop_pi0_1 = drop_pi0_1, to_pval_function = to_pval_function,
      parallel = parallel, n_workers = n_workers, verbose = verbose
    )
  }}

  n_workers = min(n_workers, n_synthetic)
  if (parallel & n_workers > 1) {
    parallel_param = BiocParallel::MulticoreParam(
      workers = n_workers,
      tasks = n_workers
    )
  } else {
    parallel_param = NULL
  }

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

  # fit working model to test statistics
  if (n_synthetic > 0) {
    working_model <- fit_working_model(
      test_statistics,
      verbose = verbose
    )
  } else if (n_synthetic == 0) {
    working_model <- NULL
    if (verbose) {
      message('No generting model fit with n_synthetic = 0')
    }
  }

  # perform grid search on simulated datasets
  grid_res <- grid_search(
    working_model = working_model,
    nsim = n_synthetic,
    synthetic_size = synthetic_size,
    ensemble_size = ensemble_size,
    to_pval_function = to_pval_function,
    locfdr_grid = locfdr_grid,
    qvalue_grid = qvalue_grid,
    fdrtool_grid = fdrtool_grid,
    focus_metric = focus_metric,
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
    to_pval_function = to_pval_function,
    focus_metric = focus_metric,
    top_grid = top_grid,
    locfdr_grid = locfdr_grid,
    fdrtool_grid = fdrtool_grid,
    qvalue_grid = qvalue_grid,
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
    'working_model' = working_model,
    'locfdr_grid' = locfdr_grid,
    'fdrtool_grid' = fdrtool_grid,
    'qvalue_grid' = qvalue_grid
  )
  return(out)
}
