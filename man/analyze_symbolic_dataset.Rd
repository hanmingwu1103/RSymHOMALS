\name{analyze_symbolic_dataset}
\alias{analyze_symbolic_dataset}
\title{Compare Symbolic HOMALS Variants on One Dataset}
\description{
Fits multiple Symbolic HOMALS variants to the same multi-valued categorical
dataset and returns a comparison table of objective values and runtime.
}
\usage{
analyze_symbolic_dataset(
  responses,
  n_categories = NULL,
  response_weights = NULL,
  ndim = 2L,
  methods = c("euclidean", "wasserstein", "hausdorff"),
  max_iter = 80L,
  tol = 1e-06,
  seed = 123,
  verbose = FALSE,
  include_mca = FALSE,
  category_labels = NULL,
  variable_names = NULL,
  object_labels = NULL
)
}
\arguments{
  \item{responses}{A list of variables containing multi-valued responses.}
  \item{n_categories}{Optional number of categories for each variable.}
  \item{response_weights}{Optional list of within-set weights, parallel to
  \code{responses}, or a character string identifying a built-in weighting
  scheme.}
  \item{ndim}{Target embedding dimension.}
  \item{methods}{Methods to compare.}
  \item{max_iter}{Maximum number of iterations.}
  \item{tol}{Relative convergence tolerance.}
  \item{seed}{Random seed used for initialization.}
  \item{verbose}{Logical flag for iteration messages.}
  \item{include_mca}{Logical flag indicating whether a naive MCA benchmark should be added.}
  \item{category_labels}{Optional list of category labels for each variable.}
  \item{variable_names}{Optional character vector of variable names.}
  \item{object_labels}{Optional character vector of object labels.}
}
\value{
A list containing fitted objects and a comparison data frame.
}
