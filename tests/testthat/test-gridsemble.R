test_statistics <- c(
  rnorm(1800, 0, 3),
  runif(100, -10, -3),
  runif(100, 3, 10)
)

test_that("gridsemble works", {
  out = gridsemble(test_statistics, verbose = TRUE)

  expect_equal(length(out$fdr), length(test_statistics))
  expect_equal(length(out$Fdr), length(test_statistics))
  expect_type(out$pi0, 'double')
  expect_equal(length(out$default_locfdr$fdr), length(test_statistics))
  expect_equal(length(out$default_fdrtool$lfdr), length(test_statistics))
  expect_equal(length(out$default_qvalue$lfdr), length(test_statistics))
  # tested top_grid and all_grids in test-grid_search
})