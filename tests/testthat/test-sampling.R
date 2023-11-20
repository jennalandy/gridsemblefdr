t = c(rnorm(900,0,1), runif(50, -5,-2), runif(50, 2, 5))

test_that("Normal fit works", {
  working_model = fit_working_model(t)

  expect_equal(names(working_model), c('parameters','thetas','type','iters'))
  expect_equal(names(working_model$parameters), c('sigmasq0','sigmasq1','pi0'))
})

test_that("simulate from Normal fit works", {
  working_model = fit_working_model(t)
  sim = simulate_from_working_model(100, working_model)

  expect_equal(sd(sim$t), sd(t), tolerance = 0.5)
  expect_equal(mean(sim$t), mean(t), tolerance = 0.5)
  expect_equal(names(sim), c('t','true_fdr','truth','p','true_Fdr','mixture'))
  expect_equal(length(sim$t), 100)
})


test_that("t fit works", {
  working_model = fit_working_model(t, df = 10, type = "t")

  expect_equal(names(working_model), c('parameters','thetas','type','iters', 'scale'))
  expect_equal(names(working_model$parameters), c('sigmasq1','pi0'))
})

test_that("simulate from t fit works", {
  working_model = fit_working_model(t, df = 10, type = "t")
  sim = simulate_from_working_model(100, working_model, df = 10)

  expect_equal(sd(sim$t), sd(t), tolerance = 0.5)
  expect_equal(mean(sim$t), mean(t), tolerance = 0.5)
  expect_equal(names(sim), c('t','true_fdr','truth','p','true_Fdr','mixture'))
  expect_equal(length(sim$t), 100)
})


test_that("t scaled fit works", {
  working_model = fit_working_model(t, df = 10, type = "t", standardize = TRUE)

  expect_equal(names(working_model), c('parameters','thetas','type','iters', 'scale'))
  expect_equal(names(working_model$parameters), c('sigmasq1','pi0'))
})

test_that("simulate from t fit works", {
  working_model = fit_working_model(t, df = 10, type = "t", standardize = TRUE)
  sim = simulate_from_working_model(100, working_model, df = 10)

  expect_equal(sd(sim$t), sd(t), tolerance = 0.5)
  expect_equal(mean(sim$t), mean(t), tolerance = 0.5)
  expect_equal(names(sim), c('t','true_fdr','truth','p','true_Fdr','mixture'))
  expect_equal(length(sim$t), 100)
})
