test_statistics <- c(
  rnorm(1800, 0, 3),
  runif(100, -10, -3),
  runif(100, 3, 10)
)

test_that("run_locfdr_row works", {
  locfdr_out <- run_locfdr_row(
    t = test_statistics,
    locfdr_grid = data.frame(
      pct = c(0),
      pct0 = c(1/3),
      nulltype = c(1),
      type = c(0)
    ),
    row = 1
  )

  expect_equal(
    names(locfdr_out),
    c('fdr','Fdr','pi0')
  )
})

test_that("run_locfdr_row works", {
  fdrtool_out <- run_fdrtool_row(
    t = test_statistics,
    fdrtool_grid = data.frame(
      cutoff.method = c('pct0'),
      pct0 = c(1/3)
    ),
    row = 1
  )

  expect_equal(
    names(fdrtool_out),
    c('fdr','Fdr','pi0')
  )
})
