#' @title Fit generating model using EM algorithm
#'
#' @param test_statistics vector, test statistics
#'
#' @param sigmasq0 double, initial value for sigmasq0
#' @param sigmasq1 double, initial value for sigmasq1
#' @param pi0 double, initial value for pi0
#'
#' @param maxiter integer, highest number of iterations of EM algorithm allowed
#' @param tol double, tolerance for change in sum of squared differences in
#' parameters in order to stop algorithm
#'
#' @return list, named parameters for the generating_model densities and values
#' across iterations
#' @noRd
fit_generating_model <- function(
  test_statistics,

  sigmasq0 = 2,
  sigmasq1 = 4,
  pi0 = 0.9,

  maxiter = 500,
  tol = 0.0001,

  verbose = TRUE
) {
  if (verbose) {
    message('Fitting generating model')
  }

  type = 'symmetric'

  # option to add other working model options later
  if (type == 'symmetric') {
    return(
      fit_generating_model_symmetric(
        test_statistics = test_statistics,

        # initialize
        sigmasq0 = sigmasq0,
        sigmasq1 = sigmasq1,
        pi0 = pi0,

        # learning parameters
        maxiter = maxiter,
        tol = tol
      )
    )
  } else {
    return(NULL)
  }
}

#' @title Fit symmetric generating model using EM algorithm
#'
#' @param test_statistics vector, test statistics
#'
#' @param sigmasq0 double, initial value for sigmasq0
#' @param sigmasq1 double, initial value for sigmasq1
#' @param pi0 double, initial value for pi0
#'
#' @param maxiter integer, max number of iterations of EM algorithm allowed
#' @param tol double, tolerance for change in sum of squared differences in
#' parameters in order to stop algorithm
#'
#' @return list, named parameters for the generating_model densities and values
#' across iterations
#' @noRd
fit_generating_model_symmetric <- function(
  test_statistics, sigmasq0, sigmasq1, pi0, maxiter, tol
) {

  prob_y1_given_t <- function(t, pi0, sigmasq0, sigmasq1) {
    (1-pi0)*(t^2/sigmasq1)*(1/sqrt(2*pi*sigmasq1))*exp(-(1/(2*sigmasq1))*t^2)/
      mix(t, pi0, sigmasq0, sigmasq1)
  }
  prob_y0_given_t <- function(t, pi0, sigmasq0, sigmasq1) {
    pi0*(1/sqrt(2*pi*sigmasq0))*exp(-(1/(2*sigmasq0))*t^2)/
      mix(t, pi0, sigmasq0, sigmasq1)
  }

  thetas <- matrix(nrow = maxiter, ncol = 3)
  i <- 1
  while( i < maxiter ) {
    thetas[i,] <- c(pi0, sigmasq0, sigmasq1)

    P_y1_given_t <- prob_y1_given_t(test_statistics, pi0, sigmasq0, sigmasq1)
    P_y0_given_t <- prob_y0_given_t(test_statistics, pi0, sigmasq0, sigmasq1)

    sigmasq0 <- sum(P_y0_given_t*test_statistics^2, na.rm = TRUE)/
      sum(P_y0_given_t, na.rm = TRUE)
    sigmasq1 <- sum(P_y1_given_t*test_statistics^2, na.rm = TRUE)/
      (3*sum(P_y1_given_t, na.rm = TRUE))

    pi0 <- mean(P_y0_given_t, na.rm = TRUE)

    i <- i + 1
    diff = sum((thetas[i,] - thetas[i-1,])^2)
    if (is.na(diff)) { break } else if (diff < tol) { break }
  }
  thetas[i,] <- c(pi0, sigmasq0, sigmasq1)
  thetas <- thetas[!is.na(thetas[,1]),]
  colnames(thetas) = c('pi0', 'sigmasq0', 'sigmasq1')

  return(list(
    'parameters' = list(sigmasq0 = sigmasq0, sigmasq1 = sigmasq1, pi0 = pi0),
    'thetas' = thetas,
    'type' = 'symmetric',
    'iters' = i
  ))
}

#' @title Simulate from generating model
#'
#' @param n integer, sample size
#' @param generating_model list, result of fit_generating_model()
#' @param df integer, degrees of freedom of test statistics, if known
#'
#' @importFrom stats quantile
#' @return data.frame, n simulated t-statistic and truth pairs
#' @noRd
simulate_from_generating_model <- function(n, generating_model, df = NULL) {

  n0 <- round(generating_model$parameters$pi0*n)

  if (generating_model$type == 'symmetric'){
    test_statistics = c(
      sample_null(
        n = n0,
        sigmasq0 = generating_model$parameters$sigmasq0
      ),
      sample_alternative(
        N = n-n0,
        sigmasq1 = generating_model$parameters$sigmasq1
      )
    )
    this_dat <- list(
      t = test_statistics,
      true_fdr = generating_model$parameters$pi0*null(
        t = test_statistics,
        sigmasq0 = generating_model$parameters$sigmasq0
      ) / mix(
        t = test_statistics,
        pi0 = generating_model$parameters$pi0,
        sigmasq0 = generating_model$parameters$sigmasq0,
        sigmasq1 = generating_model$parameters$sigmasq1
      ),
      truth = c(
        rep(0, n0),
        rep(1, n-n0)
      ),
      topq = abs(test_statistics) > stats::quantile(abs(test_statistics))['75%']
    )
    this_dat$p = p_from_t(
      test_statistics = this_dat$t,
      df = df
    )

    this_dat$true_Fdr = get_true_Fdr(
      test_statistics = this_dat$t,
      truth = this_dat$truth
    )

    this_dat$mixture = function(t) {
      mix(
        t,
        pi0 = generating_model$parameters$pi0,
        sigmasq0 = generating_model$parameters$sigmasq0,
        sigmasq1 = generating_model$parameters$sigmasq1
      )
    }

    return(this_dat)
  } else {
    return(NULL)
  }
}

# -------- Null ~ N(0, sigmasq0) --------

#' @title Null pdf
#'
#' @param t double, test statistic
#' @param sigmasq0 double, variance of null density
#'
#' @return double, value of null density pdf at value t
#' @noRd
null <- function(t, sigmasq0) {
  (2*pi*sigmasq0)^(-1/2) * exp(-t^2/(2*sigmasq0))
}

#' @title Sample from null density
#'
#' @param n integer, sample size
#' @param sigmasq0 double, variance of null distribution
#'
#' @importFrom stats rnorm
#' @return vector, test statistics sampled from null density N(0, sigmasq0)
#' @noRd
sample_null <- function(n, sigmasq0) {
  stats::rnorm(n, mean = 0, sd = sqrt(sigmasq0))
}

# -------- Alt ~ t^2/sigmasq1 N(0, sigmasq1) --------

#' @title Alternative pdf
#'
#' @param t double, test statistic
#' @param sigmasq1 double, variance of alternative distribution if symmetric
#'
#' @return double, value of alternative density pdf at value t
#' @noRd
alt <- function(t, sigmasq1) {
  (t^2/sigmasq1) * (2*pi*sigmasq1)^(-1/2) * exp(-t^2/(2*sigmasq1))
}

#' @title Sample from alternative pdf using metropolis hastings
#'
#' @param N integer, number of test statistics to sample
#' @param sigmasq1 double, variance of alternative distribution
#' @param burn_in integer, number of iterations to use as burn in and discard
#' @param t_init double, starting point for metropolis hastings
#'
#' @return vector, test statistics
#' @noRd
sample_alternative <- function(
  N,
  sigmasq1,
  burn_in = 1000,
  t_init = 1
) {
  t_vec = vector(length = N)
  t = t_init
  for (i in seq_len(N + burn_in)) {
    t = metropolis_hastings_step(t, sigmasq1)

    # use first 1000 as burn-in
    if (i > burn_in) {
      t_vec[i-burn_in] = t
    }
  }
  return (t_vec)
}

#' @title Proposal pdf
#'
#' @param of double, proposed value
#' @param given double, previous value
#' @param sd double, standard deviation of proposal distribution
#'
#' @importFrom stats dnorm
#' @return double, proposal pdf value
#' @noRd
proposal_pdf <- function(of, given, sd) {
  stats::dnorm(of, mean = given, sd = sd)
}

#' @title One Metropolis Hastings step to update t
#'
#' @param t double, current t value
#' @param sigmasq1 double, variance of alternative distribution
#'
#' @importFrom stats rnorm runif
#' @return double, next t value
#' @noRd
metropolis_hastings_step <- function(t, sigmasq1) {
  # sample t_next|t ~ q
  t_next <- stats::rnorm(1, mean = t, sd = 4)

  # compute acceptance probability
  acceptance <- min(
    1,
    alt(t_next, sigmasq1)*proposal_pdf(of = t_next, given = t, sd = 4)/
      (alt(t, sigmasq1)*proposal_pdf(of = t, given = t_next, sd = 4))
  )

  # accept t_next with acceptance probability, otherwise keep t
  u <- stats::runif(1, 0, 1)
  if (u < acceptance) {
    t <- t_next
  }

  return(t)
}

# -------- Mixture = pi0*null + (1-pi0)*alt --------

#' @title Mixture pdf
#'
#' @param t double, test statistic
#' @param pi0 double, proportion null
#' @param sigmasq0 double, variance of null distribution
#' @param sigmasq1 double, variance of alternative distribution
#'
#' @return value of mixture density pdf at value t
#' @noRd
mix <- function(t, pi0, sigmasq0, sigmasq1) {
  pi0*null(t, sigmasq0) + (1-pi0)*alt(t, sigmasq1)
}
