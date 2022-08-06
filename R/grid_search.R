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
#' @param direction string, one of c('left','right')
#'
#' @return vector, tail end Fdr values
get_true_Fdr <- function(test_statistics, truth, direction = 'left')  {
  out <- rep(NA, length(test_statistics))
  if (direction == 'left') {
    # Pr(null | T <= t) = Pr(truth = 0 | T <= t)
    for (i in 1:length(test_statistics)) {
      t = test_statistics[i]
      out[i] <- mean(1 - truth[test_statistics <= t])
    }
  } else {
    # Pr(null | T >= t)
    for (i in 1:length(test_statistics)) {
      t = test_statistics[i]
      out[i] <- mean(1 - truth[test_statistics >= t])
    }
  }
  return(out)
}

#' Grid Search
#' @description simulate data and run grid search to determine which
#' models to ensemble over on the real dataset
#'
#' @param n integer, number of test statistics in each simulated dataset
#' @param nsim integer, number of datasets to simulate.
#' If 0, no datasets are simulated and model(s) are randomly selected.
#' @param ensemble_size integer, number of models chosen to ensemble over OR double,
#' proportion of grid size to ensemble over
#' @param fit list, result of fit_sim()
#' @param df integer, degrees of freedom of test statistics, if known, or NULL
#' @param locfdr_grid data.frame, each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data.frame, each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data.frame, each row is a set of hyperparameters for qvalue
#' @param focus_metric string, one of one of c('fdrerror','Fdrerror','pr','roc','brier'),
#' which metric to optimize in the grid search
#' @param large_abs_metric boolean, if TRUE, only consider focus_metric looking at the
#' large absolute value test statistics (top quartile of abs(t))
#' @param params_type string, type of simulation model fit, one of c('symmetric')
#' @param parallel_param BiocParallel object, specified to run in parallel or NULL
#' @param verbose boolean
#'
#' @return
#' \itemize{
#'    \item fit - list, result of fit_sim()
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
  n, nsim, ensemble_size, fit,
  method_list,
  row_list,
  df = NULL,
  locfdr_grid = NULL,
  fdrtool_grid = NULL,
  qvalue_grid = NULL,
  focus_metric = 'pr',
  large_abs_metric = TRUE,
  params_type = 'symmetric',
  parallel_param = NULL,
  sim_subset = NULL,
  verbose = TRUE
) {

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

  if (ensemble_size < 1) {
    ensemble_size = max(1, round(length(method_list)*ensemble_size))
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
      'fit' = NA,
      'top_grid' = top_grid,
      'all_grids' = all_grids
    ))

  }

  # grid_size is the number of possible hyperparameter
  # combinations we're searching over
  grid_size = nrow_null0(locfdr_grid) +
    nrow_null0(fdrtool_grid) +
    nrow_null0(qvalue_grid)

  if (large_abs_metric & focus_metric %in% c('brier','Fdrerror','fdrerror','pr','roc')) {
    focus_metric = paste(focus_metric, '_topq', sep = '')
  }

  # simulate and perform grid search `nsim` times
  sample_n = ifelse(!is.null(sim_subset), sim_subset, n)
  all_grids <- do.call(rbind, parlapply(
    X = 1:nsim,
    parallel_param = parallel_param,
    FUN = function(
      sim,
      # all functions and objects that must be exported to
      # each of the parallel sessions
      sample_n,
      focus_metric,
      simulate_from_fit,
      fit,
      p_from_t,
      df,
      get_true_Fdr,
      method_list,
      row_list,
      metrics,
      run_fdrtool_row,
      run_locfdr_row,
      run_qvalue_row,
      fdrtool_grid,
      locfdr_grid,
      qvalue_grid,
      verbose
    ) {
      if(verbose) {
        print(paste('Simulation',sim))
      }

      # simulate data from fit
      this_dat <- simulate_from_fit(sample_n, fit)
      this_dat$p = p_from_t(
        test_statistics = this_dat$t,
        df = df,
        sides = 'two'
      )

      this_dat$true_Fdr = get_true_Fdr(
        test_statistics = this_dat$t,
        truth = this_dat$truth
      )

      # run a grid search for best set of parameters
      # record PR AUC, ROC AUC, and Brier Score
      # on all data and on top absolute quantile
      topq = abs(this_dat$t) > quantile(abs(this_dat$t))['75%']

      this_score <- do.call(rbind, lapply(
        1:length(row_list),
        function(i) {
          if (method_list[i] == 'locfdr') {
            row_res <- run_locfdr_row(
              test_statistics = this_dat$t,
              locfdr_grid = locfdr_grid,
              row = row_list[i]
            )
          } else if (method_list[i] == 'fdrtool') {
            row_res <- run_fdrtool_row(
              test_statistics = this_dat$t,
              fdrtool_grid = fdrtool_grid,
              row = row_list[i]
            )
          } else if (method_list[i] == 'qvalue') {
            row_res <- run_qvalue_row(
              test_statistics = this_dat$t,
              qvalue_grid = qvalue_grid,
              row = row_list[i],
              df = df
            )
          }

          # if not null, record pi0 estimate and metrics
          if (!is.null(row_res)) {

            this_metrics <- metrics(
              fdr = row_res$fdr,
              Fdr = row_res$Fdr,
              test_statistics = this_dat$t,
              truth = this_dat$truth,
              true_Fdr = this_dat$true_Fdr,
              true_fdr = this_dat$true_fdr,
              topq = topq
            )

            this_metrics$method <- method_list[i]
            this_metrics$row <- row_list[i]
            this_metrics$pi0 <- row_res$pi0
            this_metrics$sim <- sim

          } else {

            # number of metrics that would've been returned
            this_metrics <- as.list(rep(NA, 15))
            names(this_metrics) <- c(
              "pr", "roc", "brier", "Fdrerror","fdrerror",
              "pr_topq", "roc_topq", "brier_toq", "Fdrerror_topq","fdrerror_topq",
              "accuracy", "precision", "recall", "specificity", "f1"
            )

            this_metrics$method <- method_list[i]
            this_metrics$row <- row_list[i]
            this_metrics$pi0 <- NA
            this_metrics$sim <- sim
          }

          return(data.frame(this_metrics))
        }
      ))

      return(this_score)
    },
    sample_n = sample_n,
    simulate_from_fit = simulate_from_fit,
    fit = fit,
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

  grid_avgs <- all_grids %>%
    dplyr::select(-sim) %>%
    dplyr::group_by(method, row) %>%
    dplyr::summarise_all(list( ~ mean(., na.rm = TRUE)))

  # order by focus_metric, lowest to highest
  sorted <- dplyr::arrange(grid_avgs, !!rlang::sym(focus_metric))

  # if Fdrerror or brier, we want low values
  # if ROC or PR AUC, we want high values

  if(any(startsWith(
    # smaller is better for these three
    focus_metric, c('Fdrerror','brier','fdrerror')
  ))) {
    best_rows = (1:ensemble_size)
  } else {
    # get last idx where metric is not NA
    last = sum(!is.na(sorted[,focus_metric]))
    # larger is better for the rest
    best_rows = (last-ensemble_size+1):last
  }


  sorted[,'best'] = rep(FALSE, nrow(sorted))
  sorted[best_rows, 'best'] = TRUE

  return(list(
    'fit' = fit,
    'top_grid' = sorted[sorted$best,],
    'avg_grid' = sorted,
    'all_grids' = all_grids
  ))
}
