\name{simulate_symhomals_study}
\alias{simulate_symhomals_study}
\title{Run the Symbolic HOMALS Simulation Study}
\description{
Runs a multi-scenario simulation study comparing Euclidean, Wasserstein, and
Hausdorff Symbolic HOMALS.
}
\usage{
simulate_symhomals_study(
  n_rep = 12L,
  scenarios = c("geometry", "heterogeneous", "sparse"),
  methods = c("euclidean", "wasserstein", "hausdorff"),
  ndim = 2L,
  max_iter = 80L,
  tol = 1e-06,
  seed = 20260322,
  verbose = TRUE
)
}
\arguments{
  \item{n_rep}{Number of Monte Carlo replications per scenario.}
  \item{scenarios}{Simulation scenarios to evaluate.}
  \item{methods}{Methods to compare.}
  \item{ndim}{Target embedding dimension.}
  \item{max_iter}{Maximum number of iterations for each fit.}
  \item{tol}{Relative convergence tolerance.}
  \item{seed}{Base random seed.}
  \item{verbose}{Logical flag for progress messages.}
}
\value{
A list with raw replication-level results and an aggregated summary table.
}

