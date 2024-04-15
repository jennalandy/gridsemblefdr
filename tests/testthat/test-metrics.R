test_statistics <- c(
  rnorm(1800, 0, 3),
  runif(100, -10, -3),
  runif(100, 3, 10)
)

truth <- c(
  rep(0, 1800),
  rep(1, 200)
)

test_that("get_true_Fdr works", {
  true_Fdr = get_true_Fdr(
    t = test_statistics, truth = truth
  )
  expect_equal(length(true_Fdr), length(test_statistics))
  expect_equal(true_Fdr[10], mean(1 - truth[
    abs(test_statistics) >= abs(test_statistics[10])
  ]))
})

test_that("get_MSE works", {
  true_Fdr = runif(n = length(test_statistics))
  my_Fdr = rep(0, length(test_statistics))

  Fdr_error = get_MSE(my_Fdr, true_Fdr)
  expect_equal(Fdr_error, mean((true_Fdr - my_Fdr)**2))
  expect_type(Fdr_error, "double")
})

test_that("metrics works", {
  my_fdr = runif(n = length(test_statistics))
  my_true_fdr = runif(n = length(test_statistics))
  my_Fdr = c(0,0,0)
  my_true_Fdr = c(0,0,0)
  my_metrics = metrics(my_fdr, my_true_fdr, my_Fdr, my_true_Fdr)

  expect_type(my_metrics, 'list')
})




