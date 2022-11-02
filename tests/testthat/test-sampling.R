t = c(rnorm(900,0,3), runif(50, -10,-3), runif(50, 3, 10))

test_that("symmetric fit works", {
  generating_model = fit_generating_model(t)

  expect_equal(names(generating_model), c('parameters','thetas','type','iters'))
  expect_equal(names(generating_model$parameters), c('sigmasq0','sigmasq1','pi0'))
})

test_that("simulate from symmetric fit works", {
  generating_model = fit_generating_model(t)
  sim = simulate_from_generating_model(100, generating_model)

  expect_equal(names(sim), c('t','true_fdr','truth','topq','p','true_Fdr','mixture'))
  expect_equal(length(sim$t), 100)
})
