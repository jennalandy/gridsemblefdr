#' @title Fit for Simulation
#' @description Fit a null and alternative distribution in order to
#' simulate labeled data for the grid search.
#'
#' @param t vector of test statistics
#' @param type one of c('symmetric','asymmetric')
#'
#' @return vector of parameters for the fit densities
fit_sim <- function(t, type = 'symmetric') {

  neg_log_likelihood <- get_neg_log_likelihood(t, type = 'symmetric')

  # reasonable initial guesses
  bounds <- quantile(t, c(0.25, 0.75))
  inner_ts <- t[t > bounds[1] & t < bounds[2]]
  sigma_0_init <- sd(inner_ts)
  sigma_1_init <- sd(t)
  pi_0_init <- 0.9 # runif(1, min = 0.5, max = 0.999)

  parameters <- optim(
    c(sigma_0_init, sigma_1_init, pi_0_init),
    neg_log_likelihood,
    method = "L-BFGS-B",
    lower = c(0.0001, 0.1, 0),
    upper = c(Inf, Inf, 0.999)
  )$par

  if (type == 'asymmetric') {
    neg_log_likelihood <- get_neg_log_likelihood(
      t, type = 'asymmetric'
    )

    sigma_0_init <- parameters[1]
    sigma_1L_init <- parameters[2]
    sigma_1R_init <- parameters[2]
    pi_0_init <- parameters[3]
    pi_1L_init <- mean(t < 0)

    parameters <- optim(
      c(sigma_0_init, sigma_1L_init, sigma_1R_init, pi_1L_init, pi_0_init),
      neg_log_likelihood,
      method = "L-BFGS-B",
      lower = c(0.0001, 0.1, 0.1, 0, 0),
      upper = c(Inf, Inf, Inf, 1, 0.999)
    )$par
  }

  return(parameters)
}


#' @title Null density
#'
#' @param t test statistic
#' @param sigma_0 variance of null distribution
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

#' @title Alternative density
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
    )*exp( -z^2/(2*sigma_1R^2))

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
  n, sigma_1 = NULL, sigma_1L = NULL, sigma_1R = NULL,
  pi_1L = 0.5, type = 'symmetric'
) {

  pdf <- function(z) {
    alternative(z, sigma_1, sigma_1L, sigma_1R, pi_1L, type = type)
  }
  # density_sample(n, pdf, range = c(-5*sigma_1, 5*sigma_1))
  return(acceptance_sample(n, pdf, compare.params = c(0, 4)))
}

#' @title Mixture density
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
#' @param t vector of test statistics
#' @param type one of c('symmetric','asymmetric')
#'
#' @return function of parameters vector that computes NLL for the data, t
get_neg_log_likelihood <- function(Z, type = 'symmetric') {
  if (type == 'symmetric') {
    return(
      # param = c(sigma_0, sigma_1, pi_0)
      function(param) {
        -1*sum(sapply(Z, function(z) {log(mixture(
          z,
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
      -1*sum(sapply(Z, function(z) {log(mixture(
        z,
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
