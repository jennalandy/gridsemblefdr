test_that("run_locfdr_row works", {
  test_statistics <- c(
    rnorm(1800, 0, 3),
    runif(100, -10, -3),
    runif(100, 3, 10)
  )

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
