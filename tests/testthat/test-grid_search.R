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

method_list = c()
row_list = c()

default_locfdr_grid <- build_locfdr_grid(
  test_statistics = test_statistics
)
method_list = c(method_list, rep('locfdr', nrow(default_locfdr_grid)))
row_list = c(row_list, 1:nrow(default_locfdr_grid))

default_fdrtool_grid <- build_fdrtool_grid(
  test_statistics = test_statistics
)
method_list = c(method_list, rep('fdrtool', nrow(default_fdrtool_grid)))
row_list = c(row_list, 1:nrow(default_fdrtool_grid))

default_qvalue_grid <- build_qvalue_grid(
  test_statistics = test_statistics
)
method_list = c(method_list, rep('qvalue', nrow(default_qvalue_grid)))
row_list = c(row_list, 1:nrow(default_qvalue_grid))


test_that("grid_search works", {
  generating_model <- fit_generating_model(test_statistics = test_statistics)

  gs1 <- grid_search(
    generating_model = generating_model,
    nsim = 10,
    sim_size = length(test_statistics),
    ensemble_size = 1,
    df = NULL,
    locfdr_grid = default_locfdr_grid,
    fdrtool_grid = default_fdrtool_grid,
    qvalue_grid = default_qvalue_grid,
    verbose = F
  )

  expect_equal(gs1$generating_model, generating_model)
  expect_equal(nrow(gs1$top_grid), 1)
  expect_equal(nrow(gs1$all_grids),
               10*(nrow_null0(default_locfdr_grid) +
               nrow_null0(default_fdrtool_grid) +
               nrow_null0(default_qvalue_grid)))
  expect_equal(nrow(gs1$avg_grid),
               (nrow_null0(default_locfdr_grid) +
                     nrow_null0(default_fdrtool_grid) +
                     nrow_null0(default_qvalue_grid)))
  expect_true(sum(unlist(gs1$all_grids$fdrerror), na.rm = TRUE) > 0)

  gs2 <- grid_search(
    generating_model = generating_model,
    nsim = 2,
    sim_size = length(test_statistics),
    ensemble_size = 3,
    df = NULL,
    locfdr_grid = default_locfdr_grid,
    fdrtool_grid = default_fdrtool_grid,
    qvalue_grid = default_qvalue_grid,
    verbose = F
  )

  expect_equal(gs2$generating_model, generating_model)
  expect_equal(nrow(gs2$top_grid), 3)


  gs3 <- grid_search(
    generating_model = generating_model,
    nsim = 1,
    sim_size = length(test_statistics),
    ensemble_size = 1,
    df = NULL,
    locfdr_grid = default_locfdr_grid,
    fdrtool_grid = default_fdrtool_grid,
    qvalue_grid = default_qvalue_grid,
    verbose = F
  )

  expect_equal(gs3$generating_model, generating_model)
  expect_equal(nrow(gs3$top_grid), 1)
})

test_that('grid_search on random grids works', {
  fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    cutoff.method = c('fndr','pct0','locfdr'),
    pct0_range = c(0, 1/2),
    method = 'random',
    seed = 123
  )
  locfdr_grid <- build_locfdr_grid(
    test_statistics = test_statistics,
    pct_range = c(0, 0.1),
    pct0 =c(0, 1/3),
    nulltype = 1:3,
    type = 0:1,
    method = 'random',
    seed = 123
  )
  qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    transf = c('probit','logit'),
    adj = runif(n = 5, min = 0.5, max = 1.5),
    pi0.method = c('bootstrap','smoother'),
    smooth.log.pi0 = c(TRUE, FALSE)
  )

  generating_model <- fit_generating_model(test_statistics = test_statistics)

  gs <- grid_search(
    generating_model = generating_model,
    nsim = 1,
    sim_size = length(test_statistics),
    ensemble_size = 10,
    df = NULL,
    locfdr_grid = default_locfdr_grid,
    fdrtool_grid = default_fdrtool_grid,
    qvalue_grid = default_qvalue_grid,
    verbose = F
  )

  expect_equal(gs$generating_model, generating_model)
  expect_equal(nrow(gs$top_grid), 10)
  expect_equal(nrow(gs$all_grids),
               (nrow_null0(default_locfdr_grid) +
                     nrow_null0(default_fdrtool_grid) +
                     nrow_null0(default_qvalue_grid)))
  expect_true(sum(unlist(gs$all_grids$fdrerror), na.rm = TRUE) > 0)

})
