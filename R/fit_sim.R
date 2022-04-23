#' @title Fit for Simulation using EM algorithm
#' @description Fit a null and alternative distribution in order to
#' simulate labeled data for the grid search.
#'
#' @param test_statistics vector of test statistics
#' @param type one of c('symmetric','asymmetric')
#' @param sigmasq0 initial value for sigmasq0
#' @param sigmasq1 initial value for sigmasq1
#' @param pi0 initial value for pi0
#' @param maxiter highest number of iterations of EM algorithm allowed
#' @param tol tolerance for change in sum of squared differences in parameters
#' in order to stop algorithm
#'
#' @return vector of named parameters for the fit densities and values across iterations
#' @export
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
  } else if (type == 'asymmetric') {
    return(
      fit_sim_asymmetric(
        test_statistics,
        # initialize
        sigmasq0 = sigmasq0,
        sigmasq1n = sigmasq1n,
        sigmasq1p = sigmasq1p,
        pi1n = pi1n,
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

#' @title Fit for Simulation using EM algorithm, Symmetric
#' @description Fit a null and alternative distribution in order to
#' simulate labeled data for the grid search.
#'
#' @param test_statistics vector of test statistics
#' @param sigmasq0 initial value for sigmasq0
#' @param sigmasq1 initial value for sigmasq1
#' @param pi0 initial value for pi0
#' @param maxiter highest number of iterations of EM algorithm allowed
#' @param tol tolerance for change in sum of squared differences in parameters
#' in order to stop algorithm
#'
#' @return vector of named parameters for the fit densities and values across iterations
#' @export
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
    sigmasq0 <- sum(P_y0_given_t*test_statistics^2)/sum(P_y0_given_t)
    sigmasq1 <- sum(P_y1_given_t*test_statistics^2)/(3*sum(P_y1_given_t))

    pi0 <- mean(P_y0_given_t)

    thetas[i,] <- c(pi0, sigmasq0, sigmasq1)

    # check if we can break early
    diff = sum((thetas[i,] - thetas[i-1,])^2)
    if (diff < tol) {
      break
    }
  }

  thetas <- thetas[!is.na(thetas[,1]),]

  thetas <- data.frame(thetas) %>%
    mutate(i = 1:nrow(thetas))

  colnames(thetas) = c('pi0', 'sigmasq0', 'sigmasq1', 'i')

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

#' @title Fit for Simulation using EM algorithm, Asymmetric
#' @description Fit a null and alternative distribution in order to
#' simulate labeled data for the grid search.
#'
#' @param test_statistics vector of test statistics
#' @param sigmasq0 initial value for sigmasq0
#' @param sigmasq1n initial value for sigmasq1n
#' @param sigmasq1p initial value for sigmasq1p
#' @param pi1n initial value for pi1n
#' @param pi0 initial value for pi0
#' @param maxiter highest number of iterations of EM algorithm allowed
#' @param tol tolerance for change in sum of squared differences in parameters
#' in order to stop algorithm
#'
#' @return vector of named parameters for the fit densities and values across iterations
#' @export
fit_sim_asymmetric <- function(
  test_statistics,
  # initialize
  sigmasq0 = 2,
  sigmasq1n = 4,
  sigmasq1p = 4,
  pi1n = 0.5,
  pi0 = 0.9,
  # learning parameters
  maxiter = 500,
  tol = 0.0001
) {
  prob_y1_given_t <- function(t, pi0, sigmasq0, pi1n, sigmasq1n, sigmasq1p) {
    (1-pi0)*(
      ifelse(
        t < 0,
        pi1n * (t^2/sigmasq1n) * (1/sqrt(2*pi*sigmasq1n) * exp(-t^2/(2*sigmasq1n))),
        0
      ) + ifelse(
        t > 0,
        (1-pi1n) * (t^2/sigmasq1p) * (1/sqrt(2*pi*sigmasq1p) * exp(-t^2/(2*sigmasq1p))),
        0
      )
    )/mix_asymmetric(
      t, pi0, sigmasq0, pi1n, sigmasq1n, sigmasq1p
    )
  }

  prob_y0_given_t <- function(t, pi0, sigmasq0, pi1n, sigmasq1n, sigmasq1p) {
    pi0*(1/sqrt(2*pi*sigmasq0))*exp(-(1/(2*sigmasq0))*t^2)/mix_asymmetric(
      t, pi0, sigmasq0, pi1n, sigmasq1n, sigmasq1p
    )
  }

  thetas <- matrix(nrow = maxiter, ncol = 5)
  thetas[1,] <- c(pi0, sigmasq0, pi1n, sigmasq1n, sigmasq1p)

  for (i in 2:maxiter) {
    # estimate
    P_y1_given_t <- prob_y1_given_t(test_statistics, pi0, sigmasq0, pi1n, sigmasq1n, sigmasq1p)
    P_y0_given_t <- prob_y0_given_t(test_statistics, pi0, sigmasq0, pi1n, sigmasq1n, sigmasq1p)

    # maximize
    sigmasq0 <- sum(P_y0_given_t*test_statistics^2)/sum(P_y0_given_t)
    sigmasq1n <- sum((P_y1_given_t*test_statistics^2)[test_statistics < 0])/(3*sum(P_y1_given_t[test_statistics < 0]))
    sigmasq1p <- sum((P_y1_given_t*test_statistics^2)[test_statistics > 0])/(3*sum(P_y1_given_t[test_statistics > 0]))

    pi0 <- mean(P_y0_given_t)
    pi1n <- mean(P_y1_given_t[test_statistics < 0])*mean(test_statistics < 0)/mean(P_y1_given_t)

    # store
    thetas[i,] <- c(pi0, sigmasq0, pi1n, sigmasq1n, sigmasq1p)

    # check if we can break early
    diff = sum((thetas[i,] - thetas[i-1,])^2)
    if (diff < tol) {
      break
    }
  }

  # drop extra rows if iterations ended early
  thetas <- thetas[!is.na(thetas[,1]),]

  thetas <- data.frame(thetas) %>%
    mutate(i = 1:nrow(thetas))

  colnames(thetas) = c('pi0', 'sigmasq0', 'pi1n', 'sigmasq1n', 'sigmasq1p', 'i')

  parameters_list = list(
    sigmasq0 = sigmasq0,
    sigmasq1n = sigmasq1n,
    sigmasq1p = sigmasq1p,
    pi1n = pi1n,
    pi0 = pi0
  )

  return(list(
    'parameters' = parameters_list,
    'thetas' = thetas,
    'type' = 'asymmetric',
    'iters' = i
  ))
}

#' @title Null pdf
#'
#' @param t test statistic
#' @param sigmasq0 variance of null density
#'
#' @return value of null density pdf at value t
null <- function(t, sigmasq0) {
  (2*pi*sigmasq0)^(-1/2) * exp(-t^2/(2*sigmasq0))
}

#' @title Sample from null density
#'
#' @param n sample size
#' @param sigmasq0 variance of null distribution
#'
#' @return vector of test statistics sampled from null density N(0, sigmasq0)
sample_null <- function(n, sigmasq0) {
  rnorm(n, mean = 0, sd = sqrt(sigmasq0))
}

#' @title Alternative pdf
#'
#' @param t test statistic
#' @param sigmasq1 variance of alternative distribution if symmetric
#'
#' @return value of alternative density pdf at value t
alt <- function(t, sigmasq1) {
  (t^2/sigmasq1) * (2*pi*sigmasq1)^(-1/2) * exp(-t^2/(2*sigmasq1))
}

alt_neg <- function(t, sigmasq1n) {
  ifelse(
    t < 0,
    2 * (t^2/sigmasq1n) * (2*pi*sigmasq1n)^(-1/2) * exp(-t^2/(2*sigmasq1n)),
    0
  )
}
alt_pos <- function(t, sigmasq1p) {
  ifelse(
    t > 0,
    2 * (t^2/sigmasq1p) * (2*pi*sigmasq1p)^(-1/2) * exp(-t^2/(2*sigmasq1p)),
    0
  )
}
alt_asymmetric <- function(t, pi1n, sigmasq1n, sigmasq1p) {
  pi1n*alt_neg(t, sigmasq1n) + (1-pi1n)*alt_pos(t, sigmasq1p)
}

#' @title Sample from alternative pdf using metropolis hastings
#'
#' @param N number of test statistics to sample
#' @param sigmasq1 variance of alternative distribution
#' @param burn_in number of iterations to use as burn in and discard
#'
#' @return test_statistics vector of test statistics
#' @export
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

sample_alternative_asymmetric <- function(
  N,
  sigmasq1n,
  sigmasq1p,
  pi1n,
  burn_in = 1000,
  t_init = 1
) {
  t_vec = vector(length = N)
  t = t_init
  for (i in 1:(N + burn_in)) {
    t = metropolis_hastings_step_asymmetric(t, pi1n, sigmasq1n, sigmasq1p)

    # use first 1000 as burn-in
    if (i > burn_in) {
      t_vec[i-burn_in] = t
    }
  }
  return (t_vec)
}

#' @title Proposal pdf
#'
#' @param of returns
#' @param given
#' @param sd
#'
#' @return
#' @export
#'
#' @examples
proposal_pdf <- function(of, given, sd) {
  dnorm(of, mean = given, sd = sd)
}

#' @title One Metropolis Hastings step to update t
#'
#' @param t current t value
#' @param sigmasq1 variance of alternative distribution
#'
#' @return t next t value
#' @export
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

metropolis_hastings_step_asymmetric <- function(t, pi1n, sigmasq1n, sigmasq1p) {
  # sample beta_next|beta ~ q
  t_next <- rnorm(1, mean = t, sd = 4)

  # compute acceptance probability
  acceptance <- min(
    1,
    alt_asymmetric(t_next, pi1n, sigmasq1n, sigmasq1p)*proposal_pdf(of = t_next, given = t, sd = 4)/
      (alt_asymmetric(t, pi1n, sigmasq1n, sigmasq1p)*proposal_pdf(of = t, given = t_next, sd = 4))
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
#' @param t test statistic
#' @param pi0 proportion null
#' @param sigmasq0 variance of null distribution
#' @param type one of c('symmetric','asymmetric')
#'
#' @return value of mixture density pdf at value t
mix <- function(t, pi0, sigmasq0, sigmasq1) {
  pi0*null(t, sigmasq0) + (1-pi0)*alt(t, sigmasq1)
}

mix_asymmetric <- function(t, pi0, sigmasq0, pi1n, sigmasq1n, sigmasq1p) {
  pi0*null(t, sigmasq0) + (1-pi0)*alt_asymmetric(t, pi1n, sigmasq1n, sigmasq1p)
}

#' @title Simulate from fit
#' @description simulate a dataset of size n from the mixture and record truth label
#' of each value (null vs alternative)
#'
#' @param n sample size
#' @param fit result of fit_sim()
#'
#' @return data frame of n simulated t-statistic truth pairs
#' @export
simulate_from_fit <- function(n, fit) {

  n0 <- round(fit$parameters$pi0*n)

  if (fit$type == 'symmetric'){
    return(data.frame(
      t = c(
        sample_null(
          n = n0,
          sigmasq0 = fit$parameters$sigmasq0
        ),
        sample_alternative(
          N = n-n0,
          sigmasq1 = fit$parameters$sigmasq1
        )
      ),
      truth = c(
        rep(0, n0),
        rep(1, n-n0)
      )
    ))

  } else {
    return(data.frame(
      t = c(
        sample_null(
          n = n0,
          sigmasq0 = fit$parameters$sigmasq0
        ),
        sample_alternative_asymmetric(
          N = n-n0,
          sigmasq1n = fit$parameters$sigmasq1n,
          sigmasq1p = fit$parameters$sigmasq1p,
          pi1n = fit$parameters$pi1n
        )
      ),
      truth = c(
        rep(0, n0),
        rep(1, n-n0)
      )
    ))
  }
}
