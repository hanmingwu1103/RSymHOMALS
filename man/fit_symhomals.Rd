\name{fit_symhomals}
\alias{fit_symhomals}
\title{Fit a Symbolic HOMALS Model}
\description{
Fits a Symbolic HOMALS model for multi-valued categorical data using the
Euclidean, Wasserstein, or Hausdorff objective.
}
\usage{
fit_symhomals(
  responses,
  n_categories = NULL,
  ndim = 2L,
  method = c("euclidean", "wasserstein", "hausdorff"),
  max_iter = 80L,
  tol = 1e-06,
  seed = 123,
  verbose = FALSE
)
}
\arguments{
  \item{responses}{A list of variables, where each variable is a list of integer
  vectors giving the selected categories for each object.}
  \item{n_categories}{Optional number of categories for each variable.}
  \item{ndim}{Target embedding dimension.}
  \item{method}{Loss function used by Symbolic HOMALS.}
  \item{max_iter}{Maximum number of iterations.}
  \item{tol}{Relative convergence tolerance.}
  \item{seed}{Random seed used for initialization.}
  \item{verbose}{Logical flag for iteration messages.}
}
\value{
An object of class \code{symhomals_fit}.
}

