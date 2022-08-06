#' @title Fit working model using EM algorithm
#' @description Fit a null and alternative distribution in order to
#' simulate labeled data for the grid search.
#'
#' @param test_statistics vector, test statistics
#' @param type string, one of c('symmetric'), set up this way to add other working model options later
#' @param sigmasq0 double, initial value for sigmasq0
#' @param sigmasq1 double, initial value for sigmasq1
#' @param pi0 double, initial value for pi0
#' @param maxiter integer, highest number of iterations of EM algorithm allowed
#' @param tol double, tolerance for change in sum of squared differences in parameters
#' in order to stop algorithm
#'
#' @return vector of named parameters for the fit densities and values across iterations
fit_sim <- function(
  test_statistics,
  type = 'symmetric',
  # initialize
  sigmasq0 = 2,
  sigmasq1 = 4,
  sigmasq1n = 4,
  sigmasq1p = 4,
  pi1n = 0.5,
  pi0 = 0.9,
  # learning parameters
  maxiter = 500,
  tol = 0.0001
) {
  # set up this way to add other working model options later
  if (type == 'symmetric') {
    return(
      fit_sim_symmetric(
        test_statistics,
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

#' @title Fit Symmetric Working Model using EM algorithm
#' @description Fit a null and alternative distribution in order to
#' simulate labeled data for the grid search.
#'
#' @param test_statistics vector, test statistics
#' @param sigmasq0 double, initial value for sigmasq0
#' @param sigmasq1 double, initial value for sigmasq1
#' @param pi0 double, initial value for pi0
#' @param maxiter integer, max number of iterations of EM algorithm allowed
#' @param tol double, tolerance for change in sum of squared differences in parameters
#' in order to stop algorithm
#'
#' @return list, named parameters for the fit densities and values across iterations
fit_sim_symmetric <- function(
  test_statistics,
  # initialize
  sigmasq0 = 2,
  sigmasq1 = 4,
  pi0 = 0.9,
  # learning parameters
  maxiter = 500,
  tol = 0.0001
) {
  prob_y1_given_t <- function(t, pi0, sigmasq0, sigmasq1) {
    (1-pi0)*(t^2/sigmasq1)*(1/sqrt(2*pi*sigmasq1))*exp(-(1/(2*sigmasq1))*t^2)/mix(
      t, pi0, sigmasq0, sigmasq1
    )
  }
  prob_y0_given_t <- function(t, pi0, sigmasq0, sigmasq1) {
    pi0*(1/sqrt(2*pi*sigmasq0))*exp(-(1/(2*sigmasq0))*t^2)/mix(
      t, pi0, sigmasq0, sigmasq1
    )
  }

  thetas <- matrix(nrow = maxiter, ncol = 3)
  thetas[1,] <- c(pi0, sigmasq0, sigmasq1)

  for (i in 2:maxiter) {
    # estimate
    P_y1_given_t <- prob_y1_given_t(test_statistics, pi0, sigmasq0, sigmasq1)
    P_y0_given_t <- prob_y0_given_t(test_statistics, pi0, sigmasq0, sigmasq1)

    # maximize
    # MLEs, with means weighted by probabities
    sigmasq0 <- sum(P_y0_given_t*test_statistics^2, na.rm = TRUE)/sum(P_y0_given_t, na.rm = TRUE)
    sigmasq1 <- sum(P_y1_given_t*test_statistics^2, na.rm = TRUE)/(3*sum(P_y1_given_t, na.rm = TRUE))

    # definition of pi0
    pi0 <- mean(P_y0_given_t, na.rm = TRUE)

    thetas[i,] <- c(pi0, sigmasq0, sigmasq1)

    # check if we can break early
    diff = sum((thetas[i,] - thetas[i-1,])^2)
    if (is.na(diff)) {
      break
    } else if (diff < tol) {
      break
    }
  }

  thetas <- thetas[!is.na(thetas[,1]),]

  colnames(thetas) = c('pi0', 'sigmasq0', 'sigmasq1')

  parameters_list = list(
    sigmasq0 = sigmasq0,
    sigmasq1 = sigmasq1,
    pi0 = pi0
  )

  return(list(
    'parameters' = parameters_list,
    'thetas' = thetas,
    'type' = 'symmetric',
    'iters' = i
  ))
}

#' @title Null pdf
#'
#' @param t double, test statistic
#' @param sigmasq0 double, variance of null density
#'
#' @return double, value of null density pdf at value t
null <- function(t, sigmasq0) {
  (2*pi*sigmasq0)^(-1/2) * exp(-t^2/(2*sigmasq0))
}

#' @title Sample from null density
#'
#' @param n integer, sample size
#' @param sigmasq0 double, variance of null distribution
#'
#' @return vector, test statistics sampled from null density N(0, sigmasq0)
sample_null <- function(n, sigmasq0) {
  rnorm(n, mean = 0, sd = sqrt(sigmasq0))
}

#' @title Alternative pdf
#'
#' @param t double, test statistic
#' @param sigmasq1 double, variance of alternative distribution if symmetric
#'
#' @return double, value of alternative density pdf at value t
alt <- function(t, sigmasq1) {
  (t^2/sigmasq1) * (2*pi*sigmasq1)^(-1/2) * exp(-t^2/(2*sigmasq1))
}

#' @title Sample from alternative pdf using metropolis hastings
#'
#' @param N integer, number of test statistics to sample
#' @param sigmasq1 double, variance of alternative distribution
#' @param burn_in integer, number of iterations to use as burn in and discard
#'
#' @return vector, test statistics
sample_alternative <- function(
  N,
  sigmasq1,
  burn_in = 1000,
  t_init = 1
) {
  t_vec = vector(length = N)
  t = t_init
  for (i in 1:(N + burn_in)) {
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
#' @return double, proposal pdf value
proposal_pdf <- function(of, given, sd) {
  dnorm(of, mean = given, sd = sd)
}

#' @title One Metropolis Hastings step to update t
#'
#' @param t double, current t value
#' @param sigmasq1 double, variance of alternative distribution
#'
#' @return double, next t value
metropolis_hastings_step <- function(t, sigmasq1) {
  # sample t_next|t ~ q
  t_next <- rnorm(1, mean = t, sd = 4)

  # compute acceptance probability
  acceptance <- min(
    1,
    alt(t_next, sigmasq1)*proposal_pdf(of = t_next, given = t, sd = 4)/
      (alt(t, sigmasq1)*proposal_pdf(of = t, given = t_next, sd = 4))
  )

  # accept t_next with acceptance probability, otherwise keep t
  u <- runif(1, 0, 1)
  if (u < acceptance) {
    t <- t_next
  }

  return(t)
}

#' @title Mixture pdf
#'
#' @param t double, test statistic
#' @param pi0 double, proportion null
#' @param sigmasq0 double, variance of null distribution
#' @param type string, one of c('symmetric')
#'
#' @return value of mixture density pdf at value t
mix <- function(t, pi0, sigmasq0, sigmasq1) {
  pi0*null(t, sigmasq0) + (1-pi0)*alt(t, sigmasq1)
}

#' @title Simulate from fit
#' @description simulate a dataset of size n from the mixture and record truth label
#' of each value (null vs alternative)
#'
#' @param n integer, sample size
#' @param fit list, result of fit_sim()
#'
#' @return data.frame, n simulated t-statistic and truth pairs
simulate_from_fit <- function(n, fit) {

  n0 <- round(fit$parameters$pi0*n)

  if (fit$type == 'symmetric'){
    test_statistics = c(
      sample_null(
        n = n0,
        sigmasq0 = fit$parameters$sigmasq0
      ),
      sample_alternative(
        N = n-n0,
        sigmasq1 = fit$parameters$sigmasq1
      )
    )
    return(data.frame(
      t = test_statistics,
      true_fdr = fit$parameters$pi0*null(
        t = test_statistics,
        sigmasq0 = fit$parameters$sigmasq0
      ) / mix(
        t = test_statistics,
        pi0 = fit$parameters$pi0,
        sigmasq0 = fit$parameters$sigmasq0,
        sigmasq1 = fit$parameters$sigmasq1
      ),
      truth = c(
        rep(0, n0),
        rep(1, n-n0)
      )
    ))
  } else {
    return(NULL)
  }
}
