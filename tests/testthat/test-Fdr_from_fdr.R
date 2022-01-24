test_that("Fdr_from_fdr works", {
  test_statistics <- c(
    rnorm(1800, 0, 3),
    runif(100, -10, -3),
    runif(100, 3, 10)
  )
  locfdr_run <- locfdr(test_statistics)
  my_Fdr <- Fdr_from_fdr(
    fdr = locfdr_run$fdr,
    t = test_statistics,
    direction = 'left'
  )

  plot(
    test_statistics[order(test_statistics)],
    my_Fdr[order(test_statistics)],
    type = 'l',
    ylim = c(0, 1),
    xlab = 'test statistic',
    ylab = 'Fdr'
  )
  lines(
    locfdr_run$mat[,'x'],
    locfdr_run$mat[,'Fdrleft'],
    col = 'red',
    ylim = c(0,1)
  )
  legend(
    x = -1, y = 0.2,
    legend = c('locfdr','calculated'),
    col = c('red','black'),
    lty = 1,
    cex = 0.4
  )
})
