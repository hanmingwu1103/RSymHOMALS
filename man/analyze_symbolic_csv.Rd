\name{analyze_symbolic_csv}
\alias{analyze_symbolic_csv}
\title{Analyze a Set-Valued CSV Dataset}
\description{
Loads a CSV file containing set-valued symbolic variables and runs the
Symbolic HOMALS comparison workflow, optionally including a naive MCA
benchmark.
}
\usage{
analyze_symbolic_csv(
  path,
  response_weights = NULL,
  ndim = 2L,
  methods = c("euclidean", "wasserstein", "hausdorff"),
  max_iter = 80L,
  tol = 1e-06,
  seed = 123,
  verbose = FALSE,
  include_mca = TRUE
)
}
\arguments{
  \item{path}{Path to a CSV file whose symbolic columns are named with the
  pattern \code{Y1_}, \code{Y2_}, and so on.}
  \item{response_weights}{Optional list of within-set weights, parallel to
  \code{responses}, or a character string identifying a built-in weighting
  scheme.}
  \item{ndim}{Target embedding dimension.}
  \item{methods}{Methods to compare.}
  \item{max_iter}{Maximum number of iterations.}
  \item{tol}{Relative convergence tolerance.}
  \item{seed}{Random seed used for initialization.}
  \item{verbose}{Logical flag for iteration messages.}
  \item{include_mca}{Logical flag indicating whether a naive MCA benchmark
  should be included.}
}
\value{
A list containing the loaded dataset, fitted methods, comparison summaries,
and dataset-level descriptive summaries.
}
