#' Number of rows, 0 if null
#'
#' @param dataframe data.frame or NULL
#'
#' @return integer, number of rows in the dataframe, or 0 if a
#' NULL value is passed in.
#' @noRd
nrow_null0 <- function(dataframe) {
  ifelse(is.null(dataframe), 0, nrow(dataframe))
}

#' Get true tail end Fdr from test statistics and truth labels
#'
#' @param test_statistics vector, test statistics
#' @param truth vector, truth values
#'
#' @return vector, tail end Fdr values
#' @noRd
get_true_Fdr <- function(test_statistics, truth)  {
  out <- rep(NA, length(test_statistics))
  for (i in seq_len(length(test_statistics))) {
    t = test_statistics[i]
    # Pr(null | T <= t) = Pr(truth = 0 | T <= t)
    out[i] <- mean(1 - truth[test_statistics <= t])
  }
  return(out)
}

#' Single Grid Search
#' @description run grid search on a given dataset
#'
#' @param this_dat list, result of simulate_from_working_model()
#' @param df integer, degrees of freedom of test statistics, if known, or NULL
#'
#' @param method_list vector, methods to consider
#' @param row_list vector, rows of each method grid to consider
#'
#' @param fdrtool_grid data.frame, rows are possible hyperparameters for fdrtool
#' @param locfdr_grid data.frame, rows are possible hyperparameters for locfdr
#' @param qvalue_grid data.frame, rows are possible hyperparameters for qvalue
#'
#' @param verbose boolean
#'
#' @return data.frame where each row gives metrics for a model
#'
#' @importFrom dplyr arrange summarise_all group_by select
#' @importFrom rlang sym
#' @noRd
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
    seq_len(length(row_list)),
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
          true_fdr = this_dat$true_fdr
        )
        this_metrics$pi0 <- row_res$pi0
      } else {
        # o.w. record placeholder (NA) metrics and pi0 estimate
        this_metrics <- metrics(
          fdr = NULL,
          true_fdr = this_dat$true_fdr
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

#' Get Methods and Rows Lists
#'
#' @param locfdr_grid data.frame, rows are possible hyperparameters for locfdr
#' @param fdrtool_grid data.frame, rows are possible hyperparameters for fdrtool
#' @param qvalue_grid data.frame, rows are possible hyperparameters for qvalue
#'
#' @return named list with methods and rows lists
#' @noRd
get_method_row_list <- function(
    locfdr_grid,
    fdrtool_grid,
    qvalue_grid
) {

  method_list = c(
    rep('locfdr', nrow_null0(locfdr_grid)),
    rep('fdrtool', nrow_null0(fdrtool_grid)),
    rep('qvalue', nrow_null0(qvalue_grid))
  )
  row_list = unlist(lapply(
    unique(method_list),
    function(m) {
      seq_len(sum(method_list == m))
    }
  ))

  return(list(
    'methods' = method_list,
    'rows' = row_list
  ))
}

#' Grid Search
#' @description simulate data and run grid search to determine which
#' models to ensemble over on the real dataset
#'
#' @param working_model list, result of fit_working_model()
#' @param nsim integer, number of datasets to simulate.
#' If 0, no datasets are simulated and model(s) are randomly selected.
#' @param synthetic_size integer, sample size of each synthetic dataset
#' @param ensemble_size integer, number of models chosen to ensemble over
#'
#' @param df integer, degrees of freedom of test statistics, if known, or NULL
#' @param locfdr_grid data.frame, rows are possible hyperparameters for locfdr
#' @param fdrtool_grid data.frame, rows are possible hyperparameters for fdrtool
#' @param qvalue_grid data.frame, rows are possible hyperparameters for qvalue
#'
#' @param focus_metric string, one of one of c('fdrerror'),
#' which metric to optimize in the grid search
#'
#' @param parallel boolean, whether to utilize parallelization
#' @param n_workers integer, number of cores to use if parallel
#' @param parallel_param `BiocParallel` object
#'
#' @param verbose boolean
#'
#' @return
#' \itemize{
#'    \item top_grid - data.frame, each row corresponds to a model
#'    to ensemble over. Each row references a method (locfdr, fdrtool,
#'    or qvalue) and row numbers in the respective grid
#'    \item all_grids - data.frame, simulation metrics across all
#'    simulations and all models.
#' }
#'
#' @importFrom dplyr arrange summarise_all group_by select
#' @importFrom rlang sym
#' @noRd
grid_search <- function(
  working_model, nsim, synthetic_size, ensemble_size,
  df = NULL,
  locfdr_grid = NULL,
  fdrtool_grid = NULL,
  qvalue_grid = NULL,
  focus_metric = 'fdrerror',
  parallel = min(TRUE, n_workers > 1),
  n_workers =  max(parallel::detectCores() - 2, 1),
  parallel_param = NULL,
  verbose = TRUE
) {

  n_workers = min(n_workers, nsim)
  if (parallel & is.null(parallel_param) & n_workers > 1) {
    parallel_param <- BiocParallel::MulticoreParam(
      workers = n_workers,
      tasks = n_workers
    )
  }
  if (!is.null(parallel_param) & verbose) {
    message('Running grid search in parallel')
  } else if (verbose) { message('Running grid search') }

  methods_rows <- get_method_row_list(
    locfdr_grid = locfdr_grid,
    fdrtool_grid = fdrtool_grid,
    qvalue_grid = qvalue_grid
  )

  if (nsim == 0) {  # no simulation/grid_search, sample rows to ensemble
    all_grids = data.frame(
      'method' = methods_rows$methods,
      'row' = methods_rows$rows
    )
    top_grid = all_grids[
      sample(seq_len(nrow(all_grids)), size = ensemble_size),
    ]
    return(list(
      'working_model' = NA,
      'top_grid' = top_grid,
      'all_grids' = all_grids
    ))
  }

  # simulate data from working_model
  generated_dat <- list()
  for (sim in seq_len(nsim)) {
    generated_dat[[sim]] <- simulate_from_working_model(
      synthetic_size, working_model
    )
  }

  # perform grid search `nsim` times
  all_grids <- do.call(rbind, parlapply(
    X = seq_len(nsim), parallel_param = parallel_param,
    FUN = function(
      generated_dat, sim, nsim, synthetic_size, focus_metric, working_model, df,
      method_list, row_list, fdrtool_grid, locfdr_grid, qvalue_grid, verbose,
      simulate_from_working_model, p_from_t, get_true_Fdr, metrics,
      run_fdrtool_row, run_locfdr_row, run_qvalue_row
    ){
      if(verbose) { message(paste0('\tSimulation ',sim, '/', nsim)) }

      # metrics on a single generated dataset
      this_score <- single_grid_search(
        this_dat = generated_dat[[sim]], df = df,
        method_list = method_list, row_list = row_list,
        fdrtool_grid = fdrtool_grid, locfdr_grid = locfdr_grid,
        qvalue_grid = qvalue_grid, verbose = verbose
      )
      this_score$sim <- sim
      return(this_score)
    },
    generated_dat = generated_dat, synthetic_size = synthetic_size, nsim = nsim,
    simulate_from_working_model = simulate_from_working_model,
    working_model = working_model, focus_metric = focus_metric,
    p_from_t = p_from_t, df = df, get_true_Fdr = get_true_Fdr,
    method_list = methods_rows$methods, row_list = methods_rows$rows,
    metrics = metrics, run_fdrtool_row = run_fdrtool_row,
    run_locfdr_row = run_locfdr_row, run_qvalue_row = run_qvalue_row,
    fdrtool_grid = fdrtool_grid, locfdr_grid = locfdr_grid,
    qvalue_grid = qvalue_grid, verbose = verbose
  ))

  method <- sim <- NULL
  grid_avgs <- all_grids %>%
    dplyr::select(-sim) %>%
    dplyr::group_by(method, row) %>%
    dplyr::summarise_all(list( ~ mean(., na.rm = TRUE)))

  sorted <- dplyr::arrange(grid_avgs, !!rlang::sym(focus_metric))
  best_rows = seq_len(ensemble_size)

  sorted[,'best'] = rep(FALSE, nrow(sorted))
  sorted[best_rows, 'best'] = TRUE

  return(list(
    'top_grid' = sorted[sorted$best,],
    'all_grids' = all_grids
  ))
}
