test_that("Fdr_from_fdr works", {
  test_statistics <- c(
    rnorm(1800, 0, 3),
    runif(100, -10, -3),
    runif(100, 3, 10)
  )
  locfdr_run <- locfdr(test_statistics)
  my_Fdr <- Fdr_from_fdr(
    fdr = locfdr_run$fdr,
    t = test_statistics
  )

  expect_equal(length(my_Fdr), length(test_statistics))
})
