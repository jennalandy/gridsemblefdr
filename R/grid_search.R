#' Number of rows, 0 if null
#'
#' @param dataframe data.frame or NULL
#'
#' @return integer, number of rows in the dataframe, or 0 if a
#' NULL value is passed in.
nrow_null0 <- function(dataframe) {
  ifelse(is.null(dataframe), 0, nrow(dataframe))
}

#' Get true tail end Fdr from test statistics and truth labels
#'
#' @param test_statistics vector, test statistics
#' @param truth vector, truth values
#'
#' @return vector, tail end Fdr values
get_true_Fdr <- function(test_statistics, truth)  {
  out <- rep(NA, length(test_statistics))
  for (i in 1:length(test_statistics)) {
    t = test_statistics[i]
    # Pr(null | T <= t) = Pr(truth = 0 | T <= t)
    out[i] <- mean(1 - truth[test_statistics <= t])
  }
  return(out)
}

#' Single Grid Search
#' @description run grid search on a given dataset
#'
#' @param this_dat list, result of simulate_from_generating_model()
#' @param df integer, degrees of freedom of test statistics, if known, or NULL
#'
#' @param method_list vector, methods to consider
#' @param row_list vector, rows of each method grid to consider
#'
#' @param fdrtool_grid data.frame, each row is a set of hyperparameters for fdrtool
#' @param locfdr_grid data.frame, each row is a set of hyperparameters for locfdr
#' @param qvalue_grid data.frame, each row is a set of hyperparameters for qvalue
#'
#' @param verbose boolean
#'
#' @return data.frame where each row gives metrics for a model
#'
#' @importFrom dplyr arrange summarise_all group_by select
#' @importFrom rlang sym
single_grid_search <- function(
  this_dat,
  df,

  method_list,
  row_list,

  fdrtool_grid,
  locfdr_grid,
  qvalue_grid,

  verbose
) {

  # run a grid search for best set of parameters
  # record PR AUC, ROC AUC, and Brier Score
  # on all data and on top absolute quantile

  this_score <- do.call(rbind, lapply(
    1:length(row_list),
    function(i) {
      row_res = run_row(
        test_statistics = this_dat$t,
        grids = list(
          'locfdr' = locfdr_grid,
          'fdrtool' = fdrtool_grid,
          'qvalue' = qvalue_grid
        ),
        method = method_list[i],
        row = row_list[i],
        df = df
      )
      if (!is.null(row_res)) {
        # if not null, record pi0 estimate and metrics
        this_metrics <- metrics(
          fdr = row_res$fdr,
          true_fdr = this_dat$true_fdr,
          topq = this_dat$topq
        )
        this_metrics$pi0 <- row_res$pi0
      } else {
        # o.w. record placeholder (NA) metrics and pi0 estimate
        this_metrics <- metrics(
          fdr = NULL,
          true_fdr = this_dat$true_fdr,
          topq = this_dat$topq
        )
        this_metrics$pi0 <- NA
      }

      # record method and grid row
      this_metrics$method <- method_list[i]
      this_metrics$row <- row_list[i]

      return(data.frame(this_metrics))
    }
  ))
  return(this_score)
}

#' Grid Search
#' @description simulate data and run grid search to determine which
#' models to ensemble over on the real dataset
#'
#' @param generating_model list, result of fit_generating_model()
#' @param nsim integer, number of datasets to simulate.
#' If 0, no datasets are simulated and model(s) are randomly selected.
#' @param sim_size integer, sample size of each simulation.
#' @param ensemble_size integer, number of models chosen to ensemble over OR double,
#' proportion of grid size to ensemble over
#'
#' @param df integer, degrees of freedom of test statistics, if known, or NULL
#' @param locfdr_grid data.frame, each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data.frame, each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data.frame, each row is a set of hyperparameters for qvalue
#'
#' @param focus_metric string, one of one of c('fdrerror'),
#' which metric to optimize in the grid search
#' @param large_abs_metric boolean, if TRUE, only consider focus_metric looking at the
#' large absolute value test statistics (top quartile of abs(t))
#'
#' @param parallel boolean, whether to utilize parallelization
#' @param n_workers integer, number of cores to use if parallel, default 2 less than available.
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#'
#' @param verbose boolean
#'
#' @return
#' \itemize{
#'    \item generating_model - list, result of generating_model_sim()
#'    \item top_grid - data.frame, each row corresponds to a model
#'    to ensemble over. Each row references a method (locfdr, fdrtool,
#'    or qvalue) and row numbers in the respective grid
#'    \item all_grids - data.frame, simulation metrics across all
#'    simulations and all models.
#' }
#'
#' @importFrom dplyr arrange summarise_all group_by select
#' @importFrom rlang sym
grid_search <- function(
  generating_model,
  nsim,
  sim_size,
  ensemble_size,

  df = NULL,
  locfdr_grid = NULL,
  fdrtool_grid = NULL,
  qvalue_grid = NULL,

  focus_metric = 'fdrerror',
  large_abs_metric = FALSE,

  parallel = TRUE,
  n_workers =  parallel::detectCores() - 2,
  parallel_param = NULL,

  verbose = TRUE
) {

  if (ensemble_size < 1) {
    ensemble_size = max(1, round(length(method_list)*ensemble_size))
  }

  if (parallel & is.null(parallel_param)) {
    parallel_param <- BiocParallel::MulticoreParam(
      workers = n_workers,
      tasks = n_workers
    )
  }

  method_list = c()
  row_list = c()
  if (nrow_null0(locfdr_grid) > 0) {
    method_list = c(method_list, rep('locfdr', nrow(locfdr_grid)))
    row_list = c(row_list, 1:nrow(locfdr_grid))
  }
  if (nrow_null0(fdrtool_grid) > 0) {
    method_list = c(method_list, rep('fdrtool', nrow(fdrtool_grid)))
    row_list = c(row_list, 1:nrow(fdrtool_grid))
  }
  if (nrow_null0(qvalue_grid) > 0) {
    method_list = c(method_list, rep('qvalue', nrow(qvalue_grid)))
    row_list = c(row_list, 1:nrow(qvalue_grid))
  }

  if (nsim == 0) {
    # no simulation/grid_search
    # randomly sample `ensemble_size` rows to use
    all_grids = data.frame(
      'method' = method_list,
      'row' = row_list
    )
    top_grid = all_grids[
      sample(1:nrow(all_grids), size = ensemble_size),
    ]
    return(list(
      'generating_model' = NA,
      'top_grid' = top_grid,
      'all_grids' = all_grids
    ))
  }

  # grid_size is the number of possible hyperparameter
  # combinations we're searching over
  grid_size = nrow_null0(locfdr_grid) +
    nrow_null0(fdrtool_grid) +
    nrow_null0(qvalue_grid)

  if (large_abs_metric) {
    focus_metric = paste(focus_metric, '_topq', sep = '')
  }

  # simulate and perform grid search `nsim` times
  all_grids <- do.call(rbind, parlapply(
    X = 1:nsim,
    parallel_param = parallel_param,
    FUN = function(
      sim,
      sim_size,
      focus_metric,
      generating_model,
      df,
      method_list,
      row_list,
      fdrtool_grid,
      locfdr_grid,
      qvalue_grid,
      verbose,
      # all functions and objects that must be exported to
      # each of the parallel sessions
      simulate_from_generating_model,
      p_from_t,
      get_true_Fdr,
      metrics,
      run_fdrtool_row,
      run_locfdr_row,
      run_qvalue_row
    ){
      if(verbose) {
        print(paste('Simulation',sim))
      }

      # simulate data from generating_model
      this_dat <- simulate_from_generating_model(sim_size, generating_model)
      this_score <- single_grid_search(
        this_dat = this_dat,
        df = df,
        method_list = method_list,
        row_list = row_list,
        fdrtool_grid = fdrtool_grid,
        locfdr_grid = locfdr_grid,
        qvalue_grid = qvalue_grid,
        verbose = verbose
      )
      this_score$sim <- sim
      return(this_score)
    },
    sim_size = sim_size,
    simulate_from_generating_model = simulate_from_generating_model,
    generating_model = generating_model,
    focus_metric = focus_metric,
    p_from_t = p_from_t,
    df = df,
    get_true_Fdr = get_true_Fdr,
    method_list = method_list,
    row_list = row_list,
    metrics = metrics,
    run_fdrtool_row = run_fdrtool_row,
    run_locfdr_row = run_locfdr_row,
    run_qvalue_row = run_qvalue_row,
    fdrtool_grid = fdrtool_grid,
    locfdr_grid = locfdr_grid,
    qvalue_grid = qvalue_grid,
    verbose = verbose
  ))

  method <- sim <- NULL
  grid_avgs <- all_grids %>%
    dplyr::select(-sim) %>%
    dplyr::group_by(method, row) %>%
    dplyr::summarise_all(list( ~ mean(., na.rm = TRUE)))

  sorted <- dplyr::arrange(grid_avgs, !!rlang::sym(focus_metric))

  # smaller is better for fdrerror
  best_rows = (1:ensemble_size)

  sorted[,'best'] = rep(FALSE, nrow(sorted))
  sorted[best_rows, 'best'] = TRUE

  return(list(
    'generating_model' = generating_model,
    'top_grid' = sorted[sorted$best,],
    'avg_grid' = sorted,
    'all_grids' = all_grids
  ))
}
