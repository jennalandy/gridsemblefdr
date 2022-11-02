# Building grids takes in a max grid size (reduced to avoid errors on given test_statistics)
# and a method ('random' for a random search and 'grid' for evenly spaced options).
# If method = 'random', a seed can be specified.
# For either method, a lower cutoff for pi0 can be specified. This will further reduce the grid.
# A parallel parameter can be specified to test each set of parameters on given test_statistics
# across multiple cores.

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

test_that("random locfdr grid works", {

  locfdr_grid <- build_locfdr_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    method = 'random',
    seed = 111
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

test_that("parallel random locfdr grid works and is faster", {

  t1 = system.time({ locfdr_grid <- build_locfdr_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    method = 'random',
    seed = 111
  )})

  n_workers = parallel::detectCores() - 2
  parallel_param <- BiocParallel::MulticoreParam(
    workers = n_workers,
    tasks = n_workers
  )

  t2 = system.time({ locfdr_grid <- build_locfdr_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    method = 'random',
    seed = 111,
    parallel_param = parallel_param
  ) })

  expect_gt(t1['elapsed'] - t2['elapsed'], -0.05)

  expect_equal(
    names(locfdr_grid),
    c('pct','pct0','nulltype','type')
  )

  expect_gt(
    nrow(locfdr_grid),
    20
  )

})

test_that("grid locfdr grid works", {

  locfdr_grid <- build_locfdr_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    method = 'grid'
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
    grid_size = 40,
    method = 'grid'
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

test_that("random fdrtool grid works", {

  fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    lower_pi0 = 0,
    method = 'random',
    seed = 123
  )

  expect_equal(
    names(fdrtool_grid),
    c('cutoff.method','pct0')
  )
  expect_gt(
    nrow(fdrtool_grid),
    30
  )
})

test_that("random fdrtool grid works, no pct0", {

  fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    cutoff.method = c('fndr','locfdr'),
    pct0_range = c(0, 0.5),
    grid_size = 40,
    lower_pi0 = 0,
    method = 'random',
    seed = 123
  )

  expect_equal(
    names(fdrtool_grid),
    c('cutoff.method','pct0')
  )
  expect_gt(
    nrow(fdrtool_grid),
    0
  )
})


test_that("parallel random fdrtool grid works and is faster", {

  t1 = system.time({ fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    lower_pi0 = 0,
    method = 'random',
    seed = 123
  )})

  n_workers = parallel::detectCores() - 2
  parallel_param <- BiocParallel::MulticoreParam(
    workers = n_workers,
    tasks = n_workers
  )

  t2 = system.time({ fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    lower_pi0 = 0,
    method = 'random',
    seed = 123,
    parallel_param = parallel_param
  ) })

  expect_gt(t1['elapsed'] - t2['elapsed'], -0.05)

  expect_equal(
    names(fdrtool_grid),
    c('cutoff.method','pct0')
  )
  expect_gt(
    nrow(fdrtool_grid),
    30
  )

})

test_that("grid fdrtool grid works", {

  fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    method = 'grid',
    lower_pi0 = 0
  )
  expect_equal(
    names(fdrtool_grid),
    c('cutoff.method','pct0')
  )
  expect_gt(
    nrow(fdrtool_grid),
    30
  )
})

test_that("lower pi0 works for fdrtool", {

  fdrtool_grid <- build_fdrtool_grid(
    test_statistics = test_statistics,
    grid_size = 40,
    method = 'grid'
  )

  expect_equal(
    names(fdrtool_grid),
    c('cutoff.method','pct0')
  )
  expect_gt(
    nrow(fdrtool_grid),
    10
  )

})

test_that("random qvalue grid works", {
  qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    method = 'random',
    seed = 123,
    lower_pi0 = 0
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

test_that("random qvalue grid works no smoother", {
  qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    transf = c('probit', 'logit'),
    adj_range = c(0.5, 1.5),
    pi0.method = c('bootstrap'),
    smooth.log.pi0 = c(TRUE, FALSE),
    method = 'random',
    seed = 123,
    lower_pi0 = 0
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

test_that("parallel random qvalue grid works and is not slower", {

  t1 = system.time({ qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    method = 'random',
    seed = 123,
    lower_pi0 = 0
  )})

  n_workers = parallel::detectCores() - 2
  parallel_param <- BiocParallel::MulticoreParam(
    workers = n_workers,
    tasks = n_workers
  )

  t2 = system.time({ qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    method = 'random',
    seed = 123,
    lower_pi0 = 0
  ) })

  expect_gt(t1['elapsed'] - t2['elapsed'], -0.05)

  expect_equal(
    names(qvalue_grid),
    c('transf','adj','pi0.method','smooth.log.pi0')
  )
  expect_gt(
    nrow(qvalue_grid),
    20
  )

})

test_that("grid qvalue grid works", {

  qvalue_grid <- build_qvalue_grid(
    test_statistics = test_statistics,
    method = 'grid',
    lower_pi0 = 0
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
    method = 'grid',
    df = 67,
    lower_pi0 = 0
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
