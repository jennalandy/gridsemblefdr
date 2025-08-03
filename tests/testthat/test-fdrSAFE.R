test_statistics <- c(
  rnorm(800, 0, 3),
  runif(100, -10, -3),
  runif(100, 3, 10)
)

test_that("fdrSAFE works", {
  out = fdrSAFE(test_statistics, verbose = FALSE)

  expect_equal(length(out$fdr), length(test_statistics))
  expect_equal(length(out$Fdr), length(test_statistics))
  expect_type(out$pi0, 'double')
  expect_equal(length(out$default_locfdr$fdr), length(test_statistics))
  expect_equal(length(out$default_fdrtool$lfdr), length(test_statistics))
  expect_equal(length(out$default_qvalue$lfdr), length(test_statistics))
  # tested top_grid and all_grids in test-grid_search
})

test_that("fdrSAFE works not parallel", {
  out = fdrSAFE(test_statistics, verbose = FALSE, parallel = FALSE)

  expect_equal(length(out$fdr), length(test_statistics))
  expect_equal(length(out$Fdr), length(test_statistics))
  expect_type(out$pi0, 'double')
  expect_equal(length(out$default_locfdr$fdr), length(test_statistics))
  expect_equal(length(out$default_fdrtool$lfdr), length(test_statistics))
  expect_equal(length(out$default_qvalue$lfdr), length(test_statistics))
  # tested top_grid and all_grids in test-grid_search
})

test_that("fdrSAFE works n_workers = 1", {
  out = fdrSAFE(test_statistics, verbose = FALSE, n_workers = 1)

  expect_equal(length(out$fdr), length(test_statistics))
  expect_equal(length(out$Fdr), length(test_statistics))
  expect_type(out$pi0, 'double')
  expect_equal(length(out$default_locfdr$fdr), length(test_statistics))
  expect_equal(length(out$default_fdrtool$lfdr), length(test_statistics))
  expect_equal(length(out$default_qvalue$lfdr), length(test_statistics))
  # tested top_grid and all_grids in test-grid_search
})

test_that("fdrSAFE works with custom to_pval_function", {
  my_to_pval_function = function(test_statistics) {
    one_sided <- unlist(lapply(test_statistics, function(z) {
      stats::pnorm(-1*abs(z), mean = 0, sd = 3)
    }))
    2*one_sided
  }

  out = fdrSAFE(test_statistics, to_pval_function = my_to_pval_function,
                   verbose = FALSE, n_workers = 1)

  expect_equal(length(out$fdr), length(test_statistics))
  expect_equal(length(out$Fdr), length(test_statistics))
  expect_type(out$pi0, 'double')
  expect_equal(length(out$default_locfdr$fdr), length(test_statistics))
  expect_equal(length(out$default_fdrtool$lfdr), length(test_statistics))
  expect_equal(length(out$default_qvalue$lfdr), length(test_statistics))
})
