% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{Fdr_from_fdr}
\alias{Fdr_from_fdr}
\title{Fdr from fdr}
\usage{
Fdr_from_fdr(fdr, test_statistics)
}
\arguments{
\item{fdr}{vector, local false discovery rate estimates}

\item{test_statistics}{vector, test statistics
#'}
}
\value{
vector, tail-end false discovery rates
}
\description{
Calculate tail-end false discovery rate (Fdr)
from local false discovery rate (fdr)
}
