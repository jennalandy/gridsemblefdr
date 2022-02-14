test_that("nrow_null0 works", {
  n1 = nrow_null0(data.frame(a = c(1,2,3), b = c(2,3,4)))
  n2 = nrow_null0(NULL)

  expect_equal(n1, 3)
  expect_equal(n2, 0)
})

test_statistics <- c(
  rnorm(1800, 0, 3),
  runif(100, -10, -3),
  runif(100, 3, 10)
)

truth <- c(
  rep(0, 1800),
  rep(1, 200)
)

test_that("grid_search works", {
  default_locfdr_grid <- build_locfdr_grid(
    test_statistics = test_statistics
  )
  default_fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics
  )
  default_qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics
  )

  fit <- fit_sim(test_statistics = test_statistics, type = 'symmetric')

  gs1 <- grid_search(
    n = length(test_statistics),
    nsim = 10,
    topn = 1,
    fit = fit,
    df = NULL,
    locfdr_grid = default_locfdr_grid,
    fdrtool_grid = default_fdrtool_grid,
    qvalue_grid = default_qvalue_grid,
    focus_metric = 'pr',
    large_abs_metric = TRUE,
    params_type = 'symmetric',
    verbose = F
  )

  expect_equal(gs1$fit, fit)
  expect_equal(nrow(gs1$top_grid), 10)
  expect_equal(nrow(gs1$all_grids),
               10*(nrow_null0(default_locfdr_grid) +
               nrow_null0(default_fdrtool_grid) +
               nrow_null0(default_qvalue_grid)))

  gs2 <- grid_search(
    n = length(test_statistics),
    nsim = 2,
    topn = 3,
    fit = fit,
    df = NULL,
    locfdr_grid = default_locfdr_grid,
    fdrtool_grid = default_fdrtool_grid,
    qvalue_grid = default_qvalue_grid,
    focus_metric = 'pr',
    large_abs_metric = TRUE,
    params_type = 'symmetric',
    verbose = F
  )

  expect_equal(gs2$fit, fit)
  expect_equal(nrow(gs2$top_grid), 6)

  fit_asym <- fit_sim(test_statistics, type = 'symmetric')

  gs3 <- grid_search(
    n = length(test_statistics),
    nsim = 1,
    topn = 1,
    fit = fit_asym,
    df = NULL,
    locfdr_grid = default_locfdr_grid,
    fdrtool_grid = default_fdrtool_grid,
    qvalue_grid = default_qvalue_grid,
    focus_metric = 'pr',
    large_abs_metric = TRUE,
    params_type = 'asymmetric',
    verbose = F
  )

  expect_equal(gs3$fit, fit_asym)
  expect_equal(nrow(gs3$top_grid), 1)
})
