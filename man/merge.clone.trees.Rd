% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clonevol.r
\name{merge.clone.trees}
\alias{merge.clone.trees}
\title{Merge clonnal evolution trees from multiple samples into a single tree}
\usage{
\method{merge}{clone.trees}(trees, samples = NULL, sample.groups = NULL,
  merge.similar.samples = FALSE)
}
\arguments{
\item{trees:}{a list of clonal evolution trees' data frames}

\item{samples:}{name of samples that will be used in node labels}

\item{sample.groups:}{named vector of sample grouping}

\item{merge.similar.samples:}{drop a sample if there is already
another sample with the same tree}
}
\description{
Merge a list of clonal evolution trees (given by the clonal
evolution tree data frames) in multiple samples into a single clonal
evolution tree, and label the leaf nodes (identified in individual tree)
with the corresponding samples in the merged tree.
}
