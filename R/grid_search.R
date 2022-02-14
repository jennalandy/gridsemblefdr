#' Number of rows, 0 if null
#'
#' @param dataframe data frame or NULL
#'
#' @return the number of rows in the dataframe, or 0 if a
#' NULL value is passed in.
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
#' @param df degrees of freedom of test statistics, if known
#' @param locfdr_grid data frame where each row is a set of hyperparameters for locfdr
#' @param fdrtool_grid data frame where each row is a set of hyperparameters for fdrtool
#' @param qvalue_grid data frame where each row is a set of hyperparameters for qvalue
#' @param focus_metric which metric to prioritize in the grid search.
#' Must be one of c('pr','roc','brier','Fdrerror')
#' @param large_abs_metric if TRUE, only consider focus_metric looking at the
#' large absolute value test statistics (specifically, top quartile of abs(t))
#' @param params_type type of simulation model fit, one of c('symmetric','asymmetric')
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
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
grid_search <- function(
  n, nsim, topn, fit,
  methods,
  df = NULL,
  locfdr_grid = NULL,
  fdrtool_grid = NULL,
  qvalue_grid = NULL,
  focus_metric = 'pr',
  large_abs_metric = TRUE,
  params_type = 'symmetric',
  verbose = T
) {

  # grid_size is the number of possible hyperparameter
  # combinations we're searching over
  grid_size = nrow_null0(locfdr_grid) +
    nrow_null0(fdrtool_grid) +
    nrow_null0(qvalue_grid)

  all_grids = data.frame(
    'sim' = rep(1:nsim, each = grid_size),
    'method' = rep(NA, grid_size*nsim),
    'row' = rep(NA, grid_size*nsim),
    'pr'  = rep(NA, grid_size*nsim),
    'roc' = rep(NA, grid_size*nsim),
    'brier' = rep(NA, grid_size*nsim),
    'Fdrerror' = rep(NA, grid_size*nsim),
    'pr_topq'  = rep(NA, grid_size*nsim),
    'roc_topq' = rep(NA, grid_size*nsim),
    'brier_topq' = rep(NA, grid_size*nsim),
    'Fdrerror_topq' = rep(NA, grid_size*nsim),
    'pi0' = rep(NA, grid_size*nsim)
  )

  top_grid = data.frame(
    'sim' = rep(1:topn, each =  nsim),
    'method' = rep(NA, nsim*topn),
    'row' = rep(NA, nsim*topn),
    'pr'  = rep(NA, nsim*topn),
    'roc' = rep(NA, nsim*topn),
    'brier' = rep(NA, nsim*topn),
    'Fdrerror' = rep(NA, nsim*topn),
    'pr_topq'  = rep(NA, nsim*topn),
    'roc_topq' = rep(NA, nsim*topn),
    'brier_topq' = rep(NA, nsim*topn),
    'Fdrerror_topq' = rep(NA, nsim*topn),
    'pi0' = rep(NA, nsim*topn)
  )

  if (large_abs_metric) {
    focus_metric = paste(focus_metric, '_topq', sep = '')
  }

  rows <- c()
  if (!is.null(locfdr_grid)) {rows <- c(rows, 1:nrow(locfdr_grid))}
  if (!is.null(fdrtool_grid)) {rows <- c(rows, 1:nrow(fdrtool_grid))}
  if (!is.null(qvalue_grid)) {rows <- c(rows, 1:nrow(qvalue_grid))}

  # simulate and perform grid search `nsim` times
  # backend = parallel::makeCluster(nsim)
  # doParallel::registerDoParallel(backend)

  # foreach(
  #   sim=1:nsim,
  #   .export = c(
  #     'simulate_from_fit',
  #     'p_from_t',
  #     'get_true_Fdr',
  #     'nrow_null0',
  #     'run_locfdr_row',
  #     'run_fdrtool_row',
  #     'run_qvalue_row',
  #     'metrics'
  #   ),
  #   .packages = c('dplyr')
  # ) %dopar% {

  for (sim in 1:nsim) {

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

    # set up empty data frame to record method, grid row, and metrics
    this_score = data.frame(
      'method' = c(
        rep('locfdr', nrow_null0(locfdr_grid)),
        rep('fdrtool', nrow_null0(fdrtool_grid)),
        rep('qvalue', nrow_null0(qvalue_grid))
      ),
      'row' = rows,
      'pr'  = rep(NA, grid_size),
      'roc' = rep(NA, grid_size),
      'brier' = rep(NA, grid_size),
      'Fdrerror' = rep(NA, grid_size),
      'pr_topq'  = rep(NA, grid_size),
      'roc_topq' = rep(NA, grid_size),
      'brier_topq' = rep(NA, grid_size),
      'Fdrerror_topq' = rep(NA, grid_size),
      'pi0' = rep(NA, grid_size)
    )

    # run a grid search for best set of parameters
    # record PR AUC, ROC AUC, and Brier Score
    # on all data and on top absolute quantile
    topq = abs(this_dat$t) > quantile(abs(this_dat$t))['75%']

    for (i in 1:nrow(this_score)) {
      if (this_score$method[i] == 'locfdr') {
        row_res <- run_locfdr_row(
          test_statistics = this_dat$t,
          locfdr_grid = locfdr_grid,
          row = this_score$row[i]
        )
      } else if (this_score$method[i] == 'fdrtool') {
        row_res <- run_fdrtool_row(
          test_statistics = this_dat$t,
          fdrtool_grid = fdrtool_grid,
          row = this_score$row[i]
        )
      } else if (this_score$method[i] == 'qvalue') {
        row_res <- run_qvalue_row(
          test_statistics = this_dat$t,
          qvalue_grid = qvalue_grid,
          row = this_score$row[i],
          df = df
        )
      }

      # if not null, record pi0 estimate and metrics
      this_score$pi0[i] <- row_res$pi0

      this_score[
        i, c('pr','roc','brier','Fdrerror',
          'pr_topq','roc_topq','brier_topq','Fdrerror_topq')
      ] <- metrics(
        fdr = row_res$fdr,
        Fdr = row_res$Fdr,
        truth = this_dat$truth,
        true_Fdr = this_dat$true_Fdr,
        topq = topq
      )
    }
    # parallel::stopCluster(backend)

    # save full grid of metrics to all_grids
    all_grids[
      (grid_size*(sim-1) + 1):(grid_size*sim),
    ] <- cbind(
      rep(sim, nrow(this_score)),
      this_score
    )

    # order by focus_metric, choose topn best parameter combinations
    sorted <- this_score %>%
      dplyr::arrange(!!rlang::sym(focus_metric))

    if (any(startsWith(focus_metric, c('Fdrerror','brier')))) {
      # low value is best
      best <- sorted %>%
        head(topn)
    } else if (any(startsWith(focus_metric, c('pr','roc')))) {
      # high value is best
      best <- sorted %>%
        tail(topn)
    }

    # store top options in top_grid -> to be ensembled over!
    top_grid[
        (topn*(sim-1) + 1):(topn*sim),
      ] <- cbind(
        rep(sim, topn),
        best
      )
  }

  return(list(
    'fit' = fit,
    'top_grid' = top_grid,
    'all_grids' = all_grids
  ))
}
