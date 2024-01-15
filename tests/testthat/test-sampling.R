t = c(rnorm(900,0,1), runif(50, -5,-2), runif(50, 2, 5))

test_that("Normal fit works", {
  working_model = fit_working_model(t)

  expect_equal(names(working_model), c('parameters','thetas','iters'))
  expect_equal(names(working_model$parameters), c('sigmasq0','sigmasq1','pi0'))
})

test_that("simulate from Normal fit works", {
  df = 30

  working_model = fit_working_model(t)
  sim = simulate_from_working_model(
    100, working_model,
    to_pval_function = function(test_statistics) {p_from_t(test_statistics, df = df)}
  )

  expect_equal(sd(sim$t), sd(t), tolerance = 0.5)
  expect_equal(mean(sim$t), mean(t), tolerance = 0.5)
  expect_equal(names(sim), c('t','true_fdr','truth','p','true_Fdr','mixture'))
  expect_equal(length(sim$t), 100)
})
