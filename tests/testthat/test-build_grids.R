test_statistics <- c(
  rnorm(1800, 0, 2),
  runif(100, -10, -3),
  runif(100, 3, 10)
)

test_that("build_locfdr_grid works", {
  default_locfdr_grid <- build_locfdr_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    seed = 111
  )
  expect_equal(
    names(default_locfdr_grid),
    c('pct','pct0','nulltype','type')
  )
  expect_gt(
    nrow(default_locfdr_grid),
    20
  )

  default_locfdr_grid <- build_locfdr_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    method = 'grid'
  )
  expect_equal(
    names(default_locfdr_grid),
    c('pct','pct0','nulltype','type')
  )
  expect_gt(
    nrow(default_locfdr_grid),
    20
  )
})

test_that("build_fdrtool_grid works", {
  default_fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    lower_pi0 = 0,
    seed = 123
  )
  expect_equal(
    names(default_fdrtool_grid),
    c('cutoff.method','pct0')
  )
  expect_gt(
    nrow(default_fdrtool_grid),
    30
  )

  default_fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    method = 'grid',
    lower_pi0 = 0
  )
  expect_equal(
    names(default_fdrtool_grid),
    c('cutoff.method','pct0')
  )
  expect_gt(
    nrow(default_fdrtool_grid),
    30
  )

  default_fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    method = 'grid'
  )
  expect_equal(
    names(default_fdrtool_grid),
    c('cutoff.method','pct0')
  )
  expect_gt(
    nrow(default_fdrtool_grid),
    10
  )
})

test_that("build_qvalue_grid works", {
  default_qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    lower_pi0 = 0
  )
  expect_equal(
    names(default_qvalue_grid),
    c('transf','adj','pi0.method','smooth.log.pi0')
  )
  expect_gt(
    nrow(default_qvalue_grid),
    20
  )

  default_qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    method = 'grid',
    lower_pi0 = 0
  )
  expect_equal(
    names(default_qvalue_grid),
    c('transf','adj','pi0.method','smooth.log.pi0')
  )
  expect_gt(
    nrow(default_qvalue_grid),
    20
  )
})
