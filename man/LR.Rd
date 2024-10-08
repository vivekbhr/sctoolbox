% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic_regression.R
\name{LR}
\alias{LR}
\title{fit a linear model and do a likelihood ratio test}
\usage{
LR(counts, cell_group, cutoff = 0.9 * nrow(counts))
}
\arguments{
\item{counts}{vector of counts}

\item{cell_group}{vector with cell groups}

\item{cutoff}{count cutoff}
}
\value{
pvalue of LR test
}
\description{
fit a linear model and do a likelihood ratio test
}
