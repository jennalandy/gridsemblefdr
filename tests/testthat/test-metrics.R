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
  expect_equal(true_Fdr[10], mean(truth[test_statistics <= test_statistics[10]]))
})

test_that("get_Fdr_error works", {
  true_Fdr = runif(n = length(test_statistics))
  my_Fdr = rep(0, length(test_statistics))

  Fdr_error = get_Fdr_error(my_Fdr, true_Fdr)
  expect_equal(Fdr_error, mean((true_Fdr - my_Fdr)**2))
  expect_type(Fdr_error, "double")
})

test_that("get_brier works", {
  my_fdr = runif(n = length(test_statistics))
  prob_1 = 1-my_fdr
  brier = get_brier(my_fdr, truth)
  expect_equal(brier, mean((prob_1-truth)**2))
  expect_type(brier, "double")
})

test_that("get_roc works", {
  my_fdr = runif(n = length(test_statistics))
  roc = get_roc(my_fdr, truth)
  expect_type(roc, "double")
})

test_that("get_prauc works", {
  my_fdr = runif(n = length(test_statistics))
  auc = get_prauc(my_fdr, truth)
  expect_type(auc, "double")
})

test_that("metrics works", {
  topq = test_statistics >= quantile(abs(test_statistics), 0.75)
  my_fdr = runif(n = length(test_statistics))
  my_Fdr = rep(0, length(test_statistics))
  true_Fdr = true_Fdr = runif(n = length(test_statistics))
  my_metrics = metrics(my_fdr, my_Fdr, truth, true_Fdr, topq)

  expect_type(my_metrics, 'list')
})




