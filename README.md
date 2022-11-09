
# gridsemblefdr

<!-- badges: start -->
<!-- badges: end -->

## Overview

`gridsemblefdr` computes local and tail-end false discovery rates using ensemble methods. This method is introduced in the in-progress paper "Ensembling for Unsupervised Learning with Application to False Discovery Rates".

Unsupervised models present a unique challenge for hyperparameter optimization because measures of accuracy used in standard supervised techniques cannot be computed. In practice, this means hyperparameters are chosen differently by each researcher and are often left as the default values in their package of choice. We introduce a novel ensemble framework to address this issue for unsupervised problems with latent labels. This framework selects models to ensemble by their approximate performances, which are estimated using simulated labeled data informed by domain knowledge of the latent label structure. We implement our framework to improve existing false discovery rate methodology, viewing multiple hypothesis testing as an unsupervised classification problem with binary latent labels. Our simulation studies show that an ensemble outperforms three popular methods with their default hyperparameters and that, within an ensemble, combining models chosen based on their approximate performances outperforms an ensemble over a random subset of models. The R package for the false discovery rate implementation of this framework, `gridsemblefdr`, can be installed here.

## Installation
```{r setup, eval = FALSE}
# install package from GitHub:
# install.packages("devtools")
library(devtools)
devtools::install_github("jennalandy/gridsemblefdr")

library(gridsemblefdr)
```

## Quick Start in [vignette](vignettes/gridsemblefdr.Rmd)
