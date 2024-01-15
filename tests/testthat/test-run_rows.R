test_statistics <- c(
  rnorm(1800, 0, 3),
  runif(100, -10, -3),
  runif(100, 3, 10)
)

test_that("run_locfdr_row works", {
  locfdr_out <- run_locfdr_row(
    test_statistics = test_statistics,
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

test_that("run_qvalue_row works", {
  df = 30

  qvalue_out <- run_qvalue_row(
    test_statistics = test_statistics,
    to_pval_function = function(test_statistics) {p_from_t(test_statistics, df = df)},
    qvalue_grid = data.frame(
      transf = c('probit'),
      adj = c(0.8),
      pi0.method = c('bootstrap'),
      smooth.log.pi0 = c(TRUE)
    ),
    row = 1
  )

  if (is.null(qvalue_out)) {
    qvalue_out <- run_qvalue_row(
      test_statistics = test_statistics,
      to_pval_function = function(test_statistics) {p_from_t(test_statistics, df = df)},
      qvalue_grid = data.frame(
        transf = c('probit'),
        adj = c(0.8),
        pi0.method = c('bootstrap'),
        smooth.log.pi0 = c(TRUE)
      ),
      row = 1,
      lambda0 = TRUE
    )
  }

  expect_equal(
    names(qvalue_out),
    c('fdr','Fdr','pi0')
  )
})
