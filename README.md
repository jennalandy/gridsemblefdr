
# gridsemblefdr

<!-- badges: start -->
<!-- badges: end -->

The goal of gridsemblefdr is to ...

## Installation

You can install the current version of gridsemblefdr from github:

``` r
library(devtools)
devtools::install_github("jennalandy/gridsemblefdr")
```

## Example

This is a basic example calculating local and tail-end false disocvery rates from a vector of test statistics:

``` r
test_statistics <- c(
  rnorm(900, 0, 3),
  runif(50, -10, -3),
  runif(50, 3, 10)
)

library(gridsemblefdr)
gridsemblefdr(test_statistics)
```

