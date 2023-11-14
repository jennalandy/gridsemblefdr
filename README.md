
# gridsemblefdr

<!-- badges: start -->
<!-- badges: end -->

## Overview

`gridsemblefdr` computes local and tail-end false discovery rates using ensemble methods. This method is introduced in the in-progress paper "Gridsemble: Model Selection and Ensembling for False Discovery Rates".

Testing many hypotheses simultaneously requires Type I error correction, for example using local false discovery rates (fdr). There are many computational methods to estimate fdr, each with their own hyperparameters, yet in practice, there is no consensus and limited practical guidance on how models should be chosen. Because fdr is never observed, model performances cannot be computed, making model selection a challenge. In this paper, we introduce a novel model selection framework for fdr that estimates model performances with synthetic datasets. Simulation studies show that our method identifies models that outperform three popular R software packages with their default hyperparameter values. Ensembling multiple models chosen in this way yields further improvement.

## Installation
```{r setup, eval = FALSE}
# install package from GitHub:
# install.packages("devtools")
library(devtools)
devtools::install_github("jennalandy/gridsemblefdr")

library(gridsemblefdr)
```

## Quick Start in [vignette](vignettes/gridsemblefdr.Rmd)
