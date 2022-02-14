#' @title Fit for Simulation
#' @description Fit a null and alternative distribution in order to
#' simulate labeled data for the grid search.
#'
#' @param test_statistics vector of test statistics
#' @param type one of c('symmetric','asymmetric')
#'
#' @return vector of named parameters for the fit densities and type
fit_sim <- function(
  test_statistics,
  type = 'symmetric'
) {

  neg_log_likelihood <- get_neg_log_likelihood(
    test_statistics = test_statistics,
    type = 'symmetric'
  )

  # reasonable initial guesses
  bounds <- quantile(test_statistics, c(0.25, 0.75))
  inner_ts <- test_statistics[
    test_statistics > bounds[1] &
    test_statistics < bounds[2]
  ]
  sigma_0_init <- sd(inner_ts)
  sigma_1_init <- sd(test_statistics)
  pi_0_init <- 0.9

  parameters <- optim(
    c(sigma_0_init, sigma_1_init, pi_0_init),
    neg_log_likelihood,
    method = "L-BFGS-B",
    lower = c(0.0001, 0.1, 0),
    upper = c(Inf, Inf, 0.999)
  )$par

  parameters_list = list(
    sigma_0 = parameters[1],
    sigma_1 = parameters[2],
    pi0 = parameters[3]
  )

  if (type == 'asymmetric') {
    neg_log_likelihood <- get_neg_log_likelihood(
      test_statistics = test_statistics,
      type = 'asymmetric'
    )

    # inital guesses based on symmetric fit
    sigma_0_init <- parameters_list$sigma_0
    sigma_1L_init <- parameters_list$sigma_1
    sigma_1R_init <- parameters_list$sigma_1
    pi_0_init <- parameters_list$pi0
    pi_1L_init <- mean(test_statistics < 0)

    parameters <- optim(
      c(
        sigma_0_init, sigma_1L_init, sigma_1R_init,
        pi_1L_init, pi_0_init
      ),
      neg_log_likelihood,
      method = "L-BFGS-B",
      lower = c(0.0001, 0.1, 0.1, 0, 0),
      upper = c(Inf, Inf, Inf, 1, 0.999)
    )$par

    parameters_list = list(
      sigma_0 = parameters[1],
      sigma_1l = parameters[2],
      sigma_1r = parameters[3],
      pi1l = parameters[4],
      pi0 = parameters[5]
    )
  }

  return(list(
    'parameters' = parameters_list,
    'type' = type
  ))
}


#' @title Null pdf
#'
#' @param t test statistic
#' @param sigma_0 variance of null density
#'
#' @return value of null density pdf at value t
null <- function(t, sigma_0) {
  dnorm(t, mean = 0, sd = sigma_0)
}

#' @title Sample from null density
#'
#' @param n sample size
#' @param sigma_0 variance of null distribution
#'
#' @return vector of test statistics sampled from null density N(0, sigma_0)
sample_null <- function(n, sigma_0) {
  rnorm(n, mean = 0, sd = sigma_0)
}

#' @title Alternative pdf
#'
#' @param t test statistic
#' @param sigma_1 variance of alternative distribution if symmetric
#' @param sigma_1L variance of negative alternative distribution if asymmetric
#' @param sigma_1R variance of positive alternative distribution if asymmetric
#' @param pi_1L proportion of alternative distribution > 0 if asymmetric
#' @param type one of c('symmetric','asymmetric')
#'
#' @return value of alternative density pdf at value t
alternative <- function(
  t, sigma_1 = NULL, sigma_1L = NULL,
  sigma_1R = NULL, pi_1L = 0.5,
  type = 'symmetric'
) {
  if (type == 'symmetric') {

    ( (t^2/sigma_1^2) * (2*pi*sigma_1^2)^(-1/2) )*exp( -t^2/(2*sigma_1^2) )

  } else if (type == 'asymmetric') {

    pi_1L*2*as.numeric(t <= 0)*(
      (t^2/sigma_1L^2)*(2*pi*sigma_1L^2)^(-1/2)
    )*exp( -t^2/(2*sigma_1L^2)) +
    (1-pi_1L)*2*as.numeric(t > 0)*(
      (t^2/sigma_1R^2)*(2*pi*sigma_1R^2)^(-1/2)
    )*exp( -t^2/(2*sigma_1R^2))

  }
}

#' Acceptance sample from a pdf
#'
#' @param n sample size
#' @param pdf distribution to sample from
#' @param c accept/reject constant
#' @param compare_dist comparison/proposal distribution, one of c('norm','uniform')
#' @param compare_params parameters for compare_dist, either c(mean, sd) or c(min, max)
#'
#' @return vector of sampled values
acceptance_sample <- function(
  n, pdf,
  c = 10,
  compare_dist = 'norm',
  compare_params = c(0, 1)
) {
  U = runif(n = n, min = 0, max = 1)
  if (compare_dist == 'norm') {
    X = rnorm(n = n, mean = compare_params[1], sd = compare_params[2])
    accept <- (U <= pdf(X)/(
      c * dnorm(x = X, mean  = compare_params[1], sd = compare_params[2])
    ))
  } else if (compare_dist == 'unif') {
    X = runif(n = n, min = compare_params[1], max = compare_params[2])
    accept <- (U <= pdf(X)/(
      c * dunif(x = X, min  = compare_params[1], max = compare_params[2])
    ))
  }
  if (sum(accept) == n) {
    return(X)
  } else {
    return(c(X[accept], acceptance_sample(
      n - sum(accept), pdf, c = c,
      compare_dist = compare_dist,
      compare_params = compare_params
    )))
  }
}

#' @title Sample from alternative density
#' @description Using acceptance sampling
#'
#' @param n sample size
#' @param sigma_1 variance of alternative distribution if symmetric
#' @param sigma_1L variance of negative alternative distribution if asymmetric
#' @param sigma_1R variance of positive alternative distribution if asymmetric
#' @param pi_1L proportion of alternative distribution > 0 if asymmetric
#' @param type one of c('symmetric','asymmetric')
#'
#' @return vector of test statistics sampled from alternative density
sample_alternative <- function(
  n,
  sigma_1 = NULL,
  sigma_1L = NULL,
  sigma_1R = NULL,
  pi_1L = 0.5,
  type = 'symmetric'
) {

  pdf <- function(z) {
    alternative(z, sigma_1, sigma_1L, sigma_1R, pi_1L, type = type)
  }
  return(acceptance_sample(n, pdf, compare_params = c(0, 4)))
}

#' @title Mixture pdf
#'
#' @param t test statistic
#' @param sigma_0 variance of null distribution
#' @param sigma_1 variance of alternative distribution if symmetric
#' @param sigma_1L variance of negative alternative distribution if asymmetric
#' @param sigma_1R variance of positive alternative distribution if asymmetric
#' @param pi_1L proportion of alternative distribution > 0 if asymmetric
#' @param type one of c('symmetric','asymmetric')
#'
#' @return value of mixture density pdf at value t
mixture <- function(
  z, sigma_0, sigma_1, sigma_1L, sigma_1R,
  pi_0, pi_1L, type = 'symmetric'
) {
  pi_0*null(z, sigma_0) +
    (1-pi_0)*alternative(z, sigma_1, sigma_1L, sigma_1R, pi_1L, type = type)
}

#' @title Negative log likelihood of mixture density
#'
#' @param test_statistics vector of test statistics
#' @param type one of c('symmetric','asymmetric')
#'
#' @return function of parameters vector that computes NLL for the data, t
get_neg_log_likelihood <- function(
  test_statistics,
  type = 'symmetric'
) {
  if (type == 'symmetric') {
    return(
      # param = c(sigma_0, sigma_1, pi_0)
      function(param) {
        -1*sum(sapply(test_statistics, function(t) {log(mixture(
          t,
          sigma_0 = param[1],
          sigma_1 = param[2],
          pi_0 = param[3],
          type = type
        ))}))
      }
    )
  } else if (type == 'asymmetric') {
    # param = c(sigma_0, sigma_1L, sigma_1R, pi_1L, pi_0)
    function(param) {
      -1*sum(sapply(test_statistics, function(t) {log(mixture(
        t,
        sigma_0 = param[1],
        sigma_1L = param[2],
        sigma_1R = param[3],
        pi_1L = param[4],
        pi_0 = param[5],
        type = type
      ))}))
    }
  }
}


#' @title Simulate from fit
#' @description simulate a dataset of size n from the mixture and record truth label
#' of each value (null vs alternative)
#'
#' @param n sample size
#' @param fit result of fit_sim()
#'
#' @return data frame of n simulated t-statistic truth pairs
simulate_from_fit <- function(n, fit) {

  n0 <- round(fit$parameters$pi0*n)

  if (fit$type == 'symmetric'){
    return(data.frame(
      t = c(
        sample_null(
          n = n0,
          sigma_0 = fit$parameters$sigma_0
        ),
        sample_alternative(
          n = n-n0,
          sigma_1 = fit$parameters$sigma_1,
          type = 'symmetric'
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
          sigma_0 = fit$parameters$sigma_0
        ),
        sample_alternative(
          n = n-n0,
          sigma_1L = fit$parameters$sigma_1l,
          sigma_1R = fit$parameters$sigma_1r,
          type = 'asymmetric'
        )
      ),
      truth = c(
        rep(0, n0),
        rep(1, n-n0)
      )
    ))
  }
}
