# Building grids takes in a grid depth for numeric variables,
# ranges of numeric variables, and specific values of categorical variables
# (reduced to avoid errors on given test_statistics)
# A parallel parameter can be specified to test each set of parameters
# on given test_statistics across multiple cores.

test_statistics <- c(
  rnorm(1800, 0, 2),
  runif(100, -10, -3),
  runif(100, 3, 10)
)

test_statistics2 <- c(
  rnorm(1800, 0, 2),
  runif(100, -16, -8),
  runif(100, 8, 16)
)

test_that("grid locfdr grid works", {

  locfdr_grid <- build_locfdr_grid(
    test_statistics = test_statistics,
    grid_depth = 5
  )

  expect_equal(
    names(locfdr_grid),
    c('pct','pct0','nulltype','type')
  )

  expect_gt(
    nrow(locfdr_grid),
    20
  )

})

test_that("grid locfdr grid works no nulltype 2, 3", {

  locfdr_grid <- build_locfdr_grid(
    test_statistics = test_statistics,
    pct_range = c(0, 0.1),
    pct0_range = c(0, 1/3),
    nulltype = 1,
    type = 0:1,
    grid_depth = 5
  )

  expect_equal(
    names(locfdr_grid),
    c('pct','pct0','nulltype','type')
  )

  expect_gt(
    nrow(locfdr_grid),
    20
  )

})

test_that("grid fdrtool grid works", {

  fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    grid_depth = 5
  )
  expect_equal(
    names(fdrtool_grid),
    c('cutoff.method','pct0')
  )
  expect_gt(
    nrow(fdrtool_grid),
    5
  )

  fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    cutoff.method = c('fndr','locfdr'),
    pct0_range = c(0, 0.5),
    grid_depth = 5
  )

  expect_equal(
    names(fdrtool_grid),
    c('cutoff.method','pct0')
  )
  expect_equal(
    nrow(fdrtool_grid),
    2
  )
})

test_that("grid qvalue grid works", {

  qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    grid_depth = 5
  )

  expect_equal(
    names(qvalue_grid),
    c('transf','adj','pi0.method','smooth.log.pi0')
  )

  expect_gt(
    nrow(qvalue_grid),
    20
  )
})

test_that("grid qvalue grid works with df", {

  qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    grid_depth = 5,
    df = 67
  )

  expect_equal(
    names(qvalue_grid),
    c('transf','adj','pi0.method','smooth.log.pi0')
  )

  expect_gt(
    nrow(qvalue_grid),
    20
  )
})

test_that("grid qvalue grid works no smoother", {
  qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    transf = c('probit', 'logit'),
    adj_range = c(0.5, 1.5),
    pi0.method = c('bootstrap'),
    smooth.log.pi0 = c(TRUE, FALSE),
    grid_depth = 5
  )

  expect_equal(
    names(qvalue_grid),
    c('transf','adj','pi0.method','smooth.log.pi0')
  )
  expect_gt(
    nrow(qvalue_grid),
    5 # will be 10 if no models fail
  )
})
