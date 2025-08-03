#' @title Fit synthetic generator using EM algorithm
#'
#' @param test_statistics vector, test statistics
#' @param type string, type of null distribution, one of c("Normal","t")
#' @param standardize logical, whether to divide by standard deviation before fitting
#' @param standardize_by
#'
#' @param sigmasq0 double, initial value for sigmasq0
#' @param sigmasq1 double, initial value for sigmasq1
#' @param pi0 double, initial value for pi0
#'
#' @param sigmasq0_fixed logical, whether sigmasq0 should be learned or kept at sigmasq0
#' @param sigmasq1_fixed logical, whether sigmasq1 should be learned or kept at sigmasq1
#' @param pi0_fixed logical, whether pi0 should be learned or kept at pi0
#'
#' @param maxiter integer, highest number of iterations of EM algorithm allowed
#' @param tol double, tolerance for change in sum of squared differences in
#' parameters in order to stop algorithm
#'
#' @return list, named parameters for the synthetic_generator densities and values
#' across iterations
#' @export
fit_synthetic_generator <- function(
  test_statistics,

  sigmasq0 = 1,
  sigmasq1 = 2,
  pi0 = 0.9,

  sigmasq0_fixed = FALSE,
  sigmasq1_fixed = FALSE,
  pi0_fixed = FALSE,

  maxiter = 500,
  tol = 0.0001,

  verbose = TRUE,
  type = 'symmetric'
) {
  if (verbose) {
    message('Fitting synthetic generator')
  }

  # option to add other synthetic generator options later
  if (type == 'symmetric') {
    return(
      fit_synthetic_generator_z(
        test_statistics = test_statistics,

        # initialize
        sigmasq0 = sigmasq0,
        sigmasq1 = sigmasq1,
        pi0 = pi0,


        sigmasq0_fixed = sigmasq0_fixed,
        sigmasq1_fixed = sigmasq1_fixed,
        pi0_fixed = pi0_fixed,

        # learning parameters
        maxiter = maxiter,
        tol = tol
      )
    )
  } else if (type == 'asymmetric') {
    # init <- fit_synthetic_generator_z(
    #   test_statistics = test_statistics,

    #   # initialize
    #   sigmasq0 = sigmasq0,
    #   sigmasq1 = sigmasq1,
    #   pi0 = pi0,


    #   sigmasq0_fixed = sigmasq0_fixed,
    #   sigmasq1_fixed = sigmasq1_fixed,
    #   pi0_fixed = pi0_fixed,

    #   # learning parameters
    #   maxiter = maxiter,
    #   tol = tol
    # )
    # return(
    #   fit_synthetic_generator_asym(
    #     test_statistics = test_statistics,

    #     # initialize based on symmetric parameters
    #     sigmasq0 = init$parameters$sigmasq0,
    #     sigmasq1n = init$parameters$sigmasq1,
    #     sigmasq1p = init$parameters$sigmasq1,
    #     pi0 = init$parameters$pi0,
    #     pi1n = 0.5,

    #     # learning parameters
    #     maxiter = maxiter,
    #     tol = tol
    #   )
    # )

    return(
      fit_synthetic_generator_asym(
        test_statistics = test_statistics,

        # initialize based on symmetric parameters
        sigmasq0 = sigmasq0,
        sigmasq1n = sigmasq1,
        sigmasq1p = sigmasq1,
        pi0 = pi0,
        pi1n = 0.5,

        # learning parameters
        maxiter = maxiter,
        tol = tol
      )
    )
  } else {
    stop(paste("type", type, "not defined."))
  }
}

#' @title Fit normal-based synthetic generator using EM algorithm
#'
#' @param test_statistics vector, test statistics
#'
#' @param sigmasq0 double, initial value for sigmasq0
#' @param sigmasq1 double, initial value for sigmasq1
#' @param pi0 double, initial value for pi0
#'
#' @param sigmasq0_fixed logical, whether sigmasq0 should be learned or kept at sigmasq0
#' @param sigmasq1_fixed logical, whether sigmasq1 should be learned or kept at sigmasq1
#' @param pi0_fixed logical, whether pi0 should be learned or kept at pi0
#'
#' @param maxiter integer, max number of iterations of EM algorithm allowed
#' @param tol double, tolerance for change in sum of squared differences in
#' parameters in order to stop algorithm
#'
#' @return list, named parameters for the synthetic_generator densities and values
#' across iterations
#' @noRd
fit_synthetic_generator_z <- function(
  test_statistics, sigmasq0, sigmasq1, pi0,
  sigmasq0_fixed, sigmasq1_fixed, pi0_fixed,
  maxiter, tol
) {

  prob_y1_given_t <- function(t, pi0, sigmasq0, sigmasq1) {
    (1-pi0)*alt(t, sigmasq1)/
      mix(t, pi0, sigmasq0, sigmasq1)
  }
  prob_y0_given_t <- function(t, pi0, sigmasq0, sigmasq1) {
    pi0*null(t, sigmasq0)/
      mix(t, pi0, sigmasq0, sigmasq1)
  }

  thetas <- matrix(nrow = maxiter, ncol = 3)
  thetas[1,] <- c(pi0, sigmasq0, sigmasq1)
  diff = 100
  i <- 2
  while( i < maxiter ) {
    P_y1_given_t <- prob_y1_given_t(test_statistics, pi0, sigmasq0, sigmasq1)
    P_y0_given_t <- prob_y0_given_t(test_statistics, pi0, sigmasq0, sigmasq1)

    if (!sigmasq0_fixed) {
      sigmasq0 <- sum(P_y0_given_t*test_statistics^2, na.rm = TRUE)/
        sum(P_y0_given_t, na.rm = TRUE)
    }

    if (!sigmasq1_fixed) {
      sigmasq1 <- sum(P_y1_given_t*test_statistics^2, na.rm = TRUE)/
        (3*sum(P_y1_given_t, na.rm = TRUE))
    }

    if (!pi0_fixed) {
      pi0 <- mean(P_y0_given_t, na.rm = TRUE)
    }

    thetas[i,] <- c(pi0, sigmasq0, sigmasq1)
    diff = sum((thetas[i,] - thetas[i-1,])^2)
    if (is.na(diff)) { break } else if (diff < tol) { break }

    i <- i + 1
  }
  thetas <- thetas[!is.na(thetas[,1]),]
  colnames(thetas) = c('pi0', 'sigmasq0', 'sigmasq1')

  return(list(
    'parameters' = list(sigmasq0 = sigmasq0, sigmasq1 = sigmasq1, pi0 = pi0),
    'thetas' = thetas,
    'iters' = i
  ))
}

fit_synthetic_generator_asym <- function(
  test_statistics,
  sigmasq0,
  sigmasq1n,
  sigmasq1p,
  pi0,
  pi1n,
  maxiter = maxiter,
  tol = tol
) {
  prob_y1_given_t <- function(t, pi0, sigmasq0, sigmasq1n, sigmasq1p, pi1n) {
    (1-pi0)*alt_asymmetric(t, sigmasq1n, sigmasq1p, pi1n)/
      mix_asymmetric(t, pi0, sigmasq0, sigmasq1n, sigmasq1p, pi1n)
  }
  prob_y0_given_t <- function(t, pi0, sigmasq0, sigmasq1n, sigmasq1p, pi1n) {
    pi0*null(t, sigmasq0)/
      mix_asymmetric(t, pi0, sigmasq0, sigmasq1n, sigmasq1p, pi1n)
  }

  thetas <- matrix(nrow = maxiter, ncol = 5)
  thetas[1,] <- c(pi0, sigmasq0, sigmasq1n, sigmasq1p, pi1n)
  diff = 100
  i <- 2
  while( i < maxiter ) {
    P_y1_given_t <- prob_y1_given_t(test_statistics, pi0, sigmasq0, sigmasq1n, sigmasq1p, pi1n)
    P_y0_given_t <- prob_y0_given_t(test_statistics, pi0, sigmasq0, sigmasq1n, sigmasq1p, pi1n)

    sigmasq0 <- sum(P_y0_given_t*test_statistics^2, na.rm = TRUE)/
      sum(P_y0_given_t, na.rm = TRUE)

    sigmasq1n <- sum(P_y1_given_t*test_statistics^2*as.numeric(test_statistics < 0), na.rm = TRUE)/
      (3*sum(P_y1_given_t*as.numeric(test_statistics < 0), na.rm = TRUE))

    sigmasq1p <- sum(P_y1_given_t*test_statistics^2*as.numeric(test_statistics > 0), na.rm = TRUE)/
      (3*sum(P_y1_given_t*as.numeric(test_statistics > 0), na.rm = TRUE))

    pi0 <- mean(P_y0_given_t, na.rm = TRUE)
    pi1n <- mean((test_statistics < 0)*(1-P_y0_given_t), na.rm = TRUE) / mean(1-P_y0_given_t, na.rm = TRUE)

    thetas[i,] <- c(pi0, sigmasq0, sigmasq1n, sigmasq1p, pi1n)
    diff = sum((thetas[i,] - thetas[i-1,])^2)
    if (is.na(diff)) { break } else if (diff < tol) { break }

    i <- i + 1
  }
  thetas <- thetas[!is.na(thetas[,1]),]
  colnames(thetas) = c('pi0', 'sigmasq0', 'sigmasq1n', 'sigmasq1p', 'pi1n')

  return(list(
    'parameters' = list(sigmasq0 = sigmasq0, pi0 = pi0, sigmasq1n = sigmasq1n, sigmasq1p = sigmasq1p, pi1n = pi1n),
    'thetas' = thetas,
    'iters' = i
  ))
}

#' @title Simulate from synthetic generator
#'
#' @param n integer, sample size
#' @param synthetic_generator list, result of fit_synthetic_generator()
#' @param to_pval_function function, converts test statistics vector to a
#' p-value vector.
#'
#' @importFrom stats quantile
#' @return data.frame, n simulated t-statistic and truth pairs
#' @export
simulate_from_synthetic_generator <- function(n, synthetic_generator, to_pval_function = p_from_t) {

  if ("pi1n" %in% names(synthetic_generator$parameters)) {
    type = 'asymmetric'
  } else {
    type = 'symmetric'
  }

  n0 <- max(round(synthetic_generator$parameters$pi0*n), 0)
  n1 <- max(n-n0, 0)

  if (type == 'symmetric') {
    test_statistics = c(
      sample_null(
        n = n0,
        sigmasq0 = synthetic_generator$parameters$sigmasq0
      ),
      sample_alternative(
        N = n1,
        sigmasq1 = synthetic_generator$parameters$sigmasq1,
      )
    )
    true_fdr = synthetic_generator$parameters$pi0*null(
      t = test_statistics,
      sigmasq0 = synthetic_generator$parameters$sigmasq0
    ) / mix(
      t = test_statistics,
      pi0 = synthetic_generator$parameters$pi0,
      sigmasq0 = synthetic_generator$parameters$sigmasq0,
      sigmasq1 = synthetic_generator$parameters$sigmasq1
    )
    mixture = function(t) {
      mix(
        t,
        pi0 = synthetic_generator$parameters$pi0,
        sigmasq0 = synthetic_generator$parameters$sigmasq0,
        sigmasq1 = synthetic_generator$parameters$sigmasq1
      )
    }
  } else if (type == 'asymmetric') {
    test_statistics = c(
      sample_null(
        n = n0,
        sigmasq0 = synthetic_generator$parameters$sigmasq0
      ),
      sample_alternative_asymmetric(
        N = n1,
        sigmasq1n = synthetic_generator$parameters$sigmasq1n,
        sigmasq1p = synthetic_generator$parameters$sigmasq1p,
        pi1n = synthetic_generator$parameters$pi1n,
      )
    )
    true_fdr = synthetic_generator$parameters$pi0*null(
      t = test_statistics,
      sigmasq0 = synthetic_generator$parameters$sigmasq0
    ) / mix_asymmetric(
      t = test_statistics,
      pi0 = synthetic_generator$parameters$pi0,
      sigmasq0 = synthetic_generator$parameters$sigmasq0,
      sigmasq1n = synthetic_generator$parameters$sigmasq1n,
      sigmasq1p = synthetic_generator$parameters$sigmasq1p,
      pi1n = synthetic_generator$parameters$pi1n
    )
    mixture = function(t) {
      mix_asymmetric(
        t,
        pi0 = synthetic_generator$parameters$pi0,
        sigmasq0 = synthetic_generator$parameters$sigmasq0,
        sigmasq1n = synthetic_generator$parameters$sigmasq1n,
        sigmasq1p = synthetic_generator$parameters$sigmasq1p,
        pi1n = synthetic_generator$parameters$pi1n
      )
    }
  }

  this_dat <- list(
    t = test_statistics,
    true_fdr = true_fdr,
    truth = c(
      rep(0, n0),
      rep(1, n1)
    )
  )
  this_dat$p = to_pval_function(
    test_statistics = this_dat$t
  )

  this_dat$true_Fdr = get_true_Fdr(
    test_statistics = this_dat$t,
    truth = this_dat$truth
  )

  this_dat$mixture = mixture

  return(this_dat)
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
  dnorm(t, 0, sqrt(sigmasq0))
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

alt_asymmetric <- function(t, sigmasq1n, sigmasq1p, pi0L) {
  pi0L * 2 * as.numeric(t < 0) * (t^2/sigmasq1n) * (2*pi*sigmasq1n)^(-1/2) * exp(-t^2/(2*sigmasq1n)) +
  (1 - pi0L) * 2 * as.numeric(t > 0) * (t^2/sigmasq1p) * (2*pi*sigmasq1p)^(-1/2) * exp(-t^2/(2*sigmasq1p))
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

sample_alternative_asymmetric <- function(
  N,
  sigmasq1n, sigmasq1p, pi1n,
  burn_in = 1000,
  t_init = 1
) {
  t_vec = vector(length = N)
  t = t_init
  for (i in seq_len(N + burn_in)) {
    t = metropolis_hastings_step_alternative(t, sigmasq1n, sigmasq1p, pi1n)

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

metropolis_hastings_step_alternative <- function(t, sigmasq1n, sigmasq1p, pi1n) {
  # sample t_next|t ~ q
  t_next <- stats::rnorm(1, mean = t, sd = 4)

  # compute acceptance probability
  acceptance <- min(
    1,
    alt_asymmetric(t_next, sigmasq1n, sigmasq1p, pi1n)*proposal_pdf(of = t_next, given = t, sd = 4)/
      (alt_asymmetric(t, sigmasq1n, sigmasq1p, pi1n)*proposal_pdf(of = t, given = t_next, sd = 4))
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

mix_asymmetric <- function(t, pi0, sigmasq0, sigmasq1n, sigmasq1p, pi1n) {
  pi0*null(t, sigmasq0) + (1-pi0)*alt_asymmetric(t, sigmasq1n, sigmasq1p, pi1n)
}

#' @title Mixture pdf with t as null distribution
#'
#' @param t double, test statistic
#' @param pi0 double, proportion null
#' @param sigmasq1 double, variance of alternative distribution
#'
#' @return value of mixture density pdf at value t
#' @noRd
mix_t <- function(t, pi0, df, sigmasq1) {
  pi0*null_t(t, df) + (1-pi0)*alt(t, sigmasq1)
}

