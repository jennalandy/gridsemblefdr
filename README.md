# gridsemblefdr

<!-- badges: start -->

<!-- badges: end -->

## Overview

`gridsemblefdr` is a data-driven selective ensembling algorithm for estimating local (fdr) and tail-end (Fdr) false discovery rates in large-scale multiple hypothesis testing. This method is introduced in the paper [*Gridsemble: Selective Ensembling for False Discovery Rates*](https://arxiv.org/abs/2401.12865). For the code to replicate all results reported in the paper, see the [jennalandy/gridsemble_PAPER](https://github.com/jennalandy/gridsemble_PAPER) repository.

## Installation

```{r setup, eval = FALSE}
# install package from GitHub:
# install.packages("devtools")
library(devtools)
devtools::install_github("jennalandy/gridsemblefdr")

library(gridsemblefdr)
```

## Quick Start in [vignette](vignettes/gridsemblefdr.Rmd)
