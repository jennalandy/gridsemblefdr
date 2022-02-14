
# gridsemblefdr

<!-- badges: start -->
<!-- badges: end -->

## Overview

`gridsemblefdr` computes local and tail-end false discovery rates using ensemble methods.

False discovery rate methodologies are popular in multiple hypothesis testing, for example, with high-dimensional genomic data. Though most implementations use hyperparameters, it is an unsupervised problem, so standard hyperparameter optimization techniques cannot be used. In practice, this means hyperparameters are chosen differently by each researcher, and often left as the default values in their package of choice. To address this issue, we introduce a novel approach to tail end (Fdr) and local (fdr) false discovery rate estimation that builds upon existing methods in two important ways. First, we utilize ensemble methods to combine information from models fit with different sets of hyperparameters. Simulation studies show that an ensemble outperforms existing methods with default hyperparameters. Second, we determine informed sets of hyperparameters to ensemble over through a grid search on simulated, labeled test statistics based on the true, unlabeled data. Simulation studies show that in an ensemble, these informed sets of hyperparameters outperform random sets of hyperparameters. 

For more details on the methodology, see the paper associated with this package [[1]()].

## Citing this package

The methodology implemented in this package was introduced in the following publication.

**todo**

## Installation
```{r setup, eval = FALSE}
# install package from GitHub:
# install.packages("devtools")
library(devtools)
devtools::install_github("jennalandy/gridsemblefdr")

library(gridsemblefdr)
```

## Quick start guide

Given a vector of test statistics, the `gridsemblefdr` function can be used to compute local and tail-end false discovery rates.

```{r eval = FALSE}
library(gridsemblefdr)

test_statistics = c(
  rnorm(900, 0, 3),   # 90% of values near 0
  runif(50, -10, -3), # add 5% negative extreme values
  runif(50, 3, 10)    # add 5% positivive extreme values
)

gridsemble <- gridsemble(test_statistics)
```

Local and tail-end false discovery rates, as well as the proportion of test statistics from the null distribution, pi0, can be accessed from the `gridsemble` object.

```{r eval = FALSE}
fdr <- gridsemble$fdr
Fdr <- gridsemble$Fdr
pi0 <- gridsemble$pi0
```

#### Options

By default, `gridsemble` will build `nsim = 10` simulated datasets, run a grid search on each, and choose the top `topn = 1` set of hyperparameters from each to ensemble over. Further, the hyperparameter sets considered are set to the default grid values in the functions `build_locfdr_grid`, `build_fdrtool_grid`, and `build_qvalue_grid`. If a user wants a more personalized set of hyperparameters considered, they can build their own grids of options with these functions, or even exclude one of the three package from the ensemble entirely if desired.

```{r eval = FALSE}
my_locfdr_grid <- build_locfdr_grid(
  test_statistics,
  pct = c(0.001, 0.003),
  pct0 = c(1/3,1/4),
  nulltype = c(1,2),
  type = c(0)
)

my_fdrtool_grid <- build_fdrtool_grid(
  test_statistics,
  cutoff.method = c('fndr','pct0'),
  pct0 = c(1/3,1/4)
)

my_qvalue_grid <- build_qvalue_grid(
  test_statistics,
  transf = c('probit','logit'),
  adj = c(1),
  pi0.method = c('bootstrap'),
  smooth.log.pi0 = c('FALSE')
)

gridsemble <- gridsemble(
  test_statistics,
  locfdr_grid = my_locfdr_grid,
  fdrtool_grid = my_fdrtool_grid,
  qvalue_grid = my_qvalue_grid
)
```
