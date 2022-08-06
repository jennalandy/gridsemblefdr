t = c(rnorm(900,0,3), runif(50, -10,-3), runif(50, 3, 10))

test_that("symmetric fit works", {
  fit = fit_sim(t, type = 'symmetric')

  expect_equal(names(fit), c('parameters','thetas','type','iters'))
  expect_equal(names(fit$parameters), c('sigmasq0','sigmasq1','pi0'))
})

test_that("simulate from symmetric fit works", {
  fit = fit_sim(t, type = 'symmetric')
  sim = simulate_from_fit(100, fit)

  expect_equal(colnames(sim), c('t','true_fdr','truth'))
  expect_equal(nrow(sim), 100)
})
