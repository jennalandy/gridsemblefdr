#' Number of rows, 0 if null
#'
#' @param dataframe data frame or NULL
#'
#' @return the number of rows in the dataframe, or 0 if a
#' NULL value is passed in.
#'
#' @export
nrow_null0 <- function(dataframe) {
  ifelse(is.null(dataframe), 0, nrow(dataframe))
}

#' Get true tail end Fdr from test statistics and truth labels
#'
#' @param test_statistics vector of test statistics
#' @param truth bector of truth values
#' @param direction one of c('left','right')
#'
#' @return vector of tail end Fdr values
#'
#' @export
get_true_Fdr <- function(test_statistics, truth, direction = 'left')  {
  out <- rep(NA, length(test_statistics))
  if (direction == 'left') {
    for (i in 1:length(test_statistics)) {
      t = test_statistics[i]
      out[i] <- mean(truth[test_statistics <= t])
    }
  } else {
    for (i in 1:length(test_statistics)) {
      t = test_statistics[i]
      out[i] <- mean(truth[test_statistics >= t])
    }
  }
  return(out)
}

#' Grid Search
#' @description simulate data and run grid search(es) to determine which
#' models to ensemble over on the real dataset
#'
#' @param n number of test statistics in the real dataset = number of
#' test statistics to simulate
#' @param nsim Number of datasets to simulate and grid search over.
#' If 0, no datasets are simulated and model(s) are randomly selected.
#' @param topn Number of models chosen from each simulation to ensemble over.
#' If `nsim = 0`, number of random models to ensemble over.
#' @param fit result of fit_sim()
#' @param method_list
#' @param row_list
#' @param df degrees of freedom of test statistics, if known
#' @param locfdr_grid data frame where each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data frame where each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data frame where each row is a set of hyperparameters for qvalue
#' @param focus_metric which metric to prioritize in the grid search.
#' Must be one of c('pr','roc','brier','Fdrerror')
#' @param large_abs_metric if TRUE, only consider focus_metric looking at the
#' large absolute value test statistics (specifically, top quartile of abs(t))
#' @param params_type type of simulation model fit, one of c('symmetric','asymmetric')
#' @param parallel if TRUE, process is run in parallel
#' @param verbose
#'
#' @return
#' \itemize{
#'    \item fit - result of fit_sim()
#'    \item top_grid - data frame where each row corresponds to a model
#'    to ensemble over. Each row references a method (locfdr, fdrtool,
#'    or qvalue) and row numbers in the respective grid
#'    \item all_grids - data frame with simulation metrics across all
#'    simulations and all models.
#' }
#'
#' @importFrom dplyr arrange
#' @importFrom rlang sym
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#'
#' @export
grid_search <- function(
  n, nsim, topn, fit,
  method_list,
  row_list,
  df = NULL,
  locfdr_grid = NULL,
  fdrtool_grid = NULL,
  qvalue_grid = NULL,
  focus_metric = 'pr',
  large_abs_metric = TRUE,
  params_type = 'symmetric',
  parallel = TRUE,
  verbose = TRUE
) {

  # grid_size is the number of possible hyperparameter
  # combinations we're searching over
  grid_size = nrow_null0(locfdr_grid) +
    nrow_null0(fdrtool_grid) +
    nrow_null0(qvalue_grid)

  if (large_abs_metric) {
    focus_metric = paste(focus_metric, '_topq', sep = '')
  }

  # simulate and perform grid search `nsim` times
  n.cores = ifelse(
    parallel,
    min(nsim, parallel::detectCores() - 1),
    1
  )
  cl = parallel::makeCluster(
    n.cores,
    type = 'PSOCK'
  )
  doParallel::registerDoParallel(cl)

  all_grids <- foreach(
    sim=1:nsim,
    .packages = c(
      'dplyr',
      'magrittr',
      'foreach'
    ),
    .export = c(
      "simulate_from_fit"
    ),
    .combine = rbind
  ) %dopar% {

    if(verbose) {
      print(paste('Simulation',sim))
    }

    # simulate data from fit
    this_dat <- simulate_from_fit(n, fit)
    this_dat$p = p_from_t(
      test_statistics = this_dat$t,
      df = df,
      sides = 'two'
    )
    this_dat$true_Fdr = get_true_Fdr(
      t = this_dat$t,
      truth = this_dat$truth
    )

    # run a grid search for best set of parameters
    # record PR AUC, ROC AUC, and Brier Score
    # on all data and on top absolute quantile
    topq = abs(this_dat$t) > quantile(abs(this_dat$t))['75%']

    this_score <- foreach (
      i = 1:length(row_list),
      .combine = rbind,
      .packages = c(
        "magrittr"
      )
    ) %dopar% {
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
          truth = this_dat$truth,
          true_Fdr = this_dat$true_Fdr,
          topq = topq
        )

        this_metrics$method <- method_list[i]
        this_metrics$row <- row_list[i]
        this_metrics$pi0 <- row_res$pi0
        this_metrics$sim <- sim

      } else {

        this_metrics <- list(
          "pr"=NA,
          "roc"= NA,
          "brier"=NA,
          "Fdrerror"=NA,
          "pr_topq"=NA,
          "roc_topq"=NA,
          "brier_toq"=NA,
          "Fdrerror_topq"=NA,
          "method"=method_list[i],
          "row"=row_list[i],
          "pi0"=NA,
          "sim"=sim
        )

      }
      return(data.frame(this_metrics))
    }

    # order by focus_metric, lowest to highest
    sorted <- dplyr::arrange(this_score, !!rlang::sym(focus_metric))

    return(sorted)
  }
  parallel::stopCluster(cl)

  # if Fdrerror or brier, we want low values
  # if ROC or PR AUC, we want high values
  if(any(startsWith(
      focus_metric, c('Fdrerror','brier')
  ))) {
    best_rows = (1:topn)
  } else {
    best_rows = (grid_size-topn+1):grid_size
  }

  top_grid = data.frame(do.call(
    # concatenate the topn within each sim
    rbind,
    lapply(
      1:nsim,
      # already sorted within sim, take best n from each
      function(i) {all_grids[all_grids$sim == i,][best_rows,]}
    )
  ))

  return(list(
    'fit' = fit,
    'top_grid' = top_grid,
    'all_grids' = all_grids
  ))
}
