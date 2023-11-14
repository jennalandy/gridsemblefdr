t = c(rnorm(900,0,3), runif(50, -10,-3), runif(50, 3, 10))

test_that("Normal fit works", {
  working_model = fit_working_model(t)

  expect_equal(names(working_model), c('parameters','thetas','type','iters'))
  expect_equal(names(working_model$parameters), c('sigmasq0','sigmasq1','pi0'))
})

test_that("simulate from Normal fit works", {
  working_model = fit_working_model(t)
  sim = simulate_from_working_model(100, working_model)

  expect_equal(names(sim), c('t','true_fdr','truth','p','true_Fdr','mixture'))
  expect_equal(length(sim$t), 100)
})


test_that("t fit works", {
  working_model = fit_working_model(t, df = 10, type = "t")

  expect_equal(names(working_model), c('parameters','thetas','type','iters'))
  expect_equal(names(working_model$parameters), c('sigmasq1','pi0'))
})

test_that("simulate from t fit works", {
  working_model = fit_working_model(t, df = 10, type = "t")
  sim = simulate_from_working_model(100, working_model, df = df)

  expect_equal(names(sim), c('t','true_fdr','truth','p','true_Fdr','mixture'))
  expect_equal(length(sim$t), 100)
})

