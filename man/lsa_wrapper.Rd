% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TF-IDF_functions.R
\name{lsa_wrapper}
\alias{lsa_wrapper}
\title{Wrapper for LSA}
\usage{
lsa_wrapper(counts, nPC, nK, umap_spread = 5, umap_mindist = 0.1, outPdf = NA)
}
\arguments{
\item{counts}{a count matrix of genes/regions (rows) * cells(columns), or a tibble of "gene", "cell" and "count"}

\item{nPC}{No. of Dims to reduce the data to}

\item{nK}{No. of neighbours louvain clustering and UMAP}

\item{outPdf}{File name to save the output UMAP}
}
\value{
list object with pca and umap output for both cells and regions.
}
\description{
Wrapper for LSA
}
