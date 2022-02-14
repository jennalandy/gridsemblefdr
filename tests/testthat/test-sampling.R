t = c(rnorm(900,0,3), runif(50, -10,-3), runif(50, 3, 10))

test_that("acceptance sampling works", {
  # use acceptance sampling for U(1,2) distribution
  acceptance_sampled <- acceptance_sample(
    n = 1000,
    pdf = function(x) {ifelse(1<x & x<2, 1, 0)},
    compare_dist = 'unif',
    compare_params = c(0,3)
  )
  real_sampled <- runif(1000, 1, 2)
  qqplot(acceptance_sampled, real_sampled)
  abline(a = 0, b = 1, col = 'red')

  expect_equal(length(acceptance_sampled), 1000)
})

test_that("symmetric fit works", {
  fit = fit_sim(t, type = 'symmetric')

  expect_equal(names(fit), c('parameters','type'))
  expect_equal(names(fit$parameters), c('sigma_0','sigma_1','pi0'))
})

test_that("asymmetric fit works", {
  fit = fit_sim(t, type = 'asymmetric')

  expect_equal(names(fit), c('parameters','type'))
  expect_equal(names(fit$parameters), c('sigma_0','sigma_1l','sigma_1r','pi1l','pi0'))
})

test_that("simulate from symmetric fit works", {
  fit = fit_sim(t, type = 'symmetric')
  sim = simulate_from_fit(100, fit)

  expect_equal(colnames(sim), c('t','truth'))
  expect_equal(nrow(sim), 100)
})

test_that("simulate from asymmetric fit works", {
  fit = fit_sim(t, type = 'asymmetric')
  sim = simulate_from_fit(100, fit)

  expect_equal(colnames(sim), c('t','truth'))
  expect_equal(nrow(sim), 100)
})
