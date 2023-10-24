t = c(rnorm(900,0,3), runif(50, -10,-3), runif(50, 3, 10))

test_that("symmetric fit works", {
  working_model = fit_working_model(t)

  expect_equal(names(working_model), c('parameters','thetas','type','iters'))
  expect_equal(names(working_model$parameters), c('sigmasq0','sigmasq1','pi0'))
})

test_that("simulate from symmetric fit works", {
  working_model = fit_working_model(t)
  sim = simulate_from_working_model(100, working_model)

  expect_equal(names(sim), c('t','true_fdr','truth','p','true_Fdr','mixture'))
  expect_equal(length(sim$t), 100)
})
