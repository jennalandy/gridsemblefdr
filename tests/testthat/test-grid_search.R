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

default_locfdr_grid <- build_locfdr_grid(
  test_statistics = test_statistics
)

default_fdrtool_grid <- build_fdrtool_grid(
  test_statistics = test_statistics
)

default_qvalue_grid <- build_qvalue_grid(
  test_statistics = test_statistics
)


test_that("grid_search works", {
  synthetic_generator <- fit_synthetic_generator(
    test_statistics = test_statistics,
    verbose = FALSE
  )

  gs1 <- grid_search(
    synthetic_generator = synthetic_generator,
    nsim = 10,
    synthetic_size = length(test_statistics),
    ensemble_size = 1,
    df = NULL,
    locfdr_grid = default_locfdr_grid,
    fdrtool_grid = default_fdrtool_grid,
    qvalue_grid = default_qvalue_grid,
    verbose = FALSE
  )

  expect_equal(nrow(gs1$top_grid), 1)
  expect_equal(nrow(gs1$all_grids),
               10*(nrow_null0(default_locfdr_grid) +
               nrow_null0(default_fdrtool_grid) +
               nrow_null0(default_qvalue_grid)))
  expect_true(sum(unlist(gs1$all_grids$fdrerror), na.rm = TRUE) > 0)

  gs2 <- grid_search(
    synthetic_generator = synthetic_generator,
    nsim = 2,
    synthetic_size = length(test_statistics),
    ensemble_size = 3,
    df = NULL,
    locfdr_grid = default_locfdr_grid,
    fdrtool_grid = default_fdrtool_grid,
    qvalue_grid = default_qvalue_grid,
    verbose = FALSE
  )

  expect_equal(nrow(gs2$top_grid), 3)


  gs3 <- grid_search(
    synthetic_generator = synthetic_generator,
    nsim = 1,
    synthetic_size = length(test_statistics),
    ensemble_size = 1,
    df = NULL,
    locfdr_grid = default_locfdr_grid,
    fdrtool_grid = default_fdrtool_grid,
    qvalue_grid = default_qvalue_grid,
    verbose = FALSE
  )

  expect_equal(nrow(gs3$top_grid), 1)
})
