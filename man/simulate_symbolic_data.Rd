\name{simulate_symbolic_data}
\alias{simulate_symbolic_data}
\title{Simulate Multi-Valued Symbolic Categorical Data}
\description{
Generates latent low-dimensional object scores, category quantifications, and
multi-valued responses for Symbolic HOMALS studies.
}
\usage{
simulate_symbolic_data(
  n = NULL,
  m = NULL,
  k = NULL,
  ndim = 2L,
  scenario = c("geometry", "heterogeneous", "sparse"),
  seed = NULL
)
}
\arguments{
  \item{n}{Number of objects.}
  \item{m}{Number of variables.}
  \item{k}{Number of categories per variable.}
  \item{ndim}{Latent geometric dimension.}
  \item{scenario}{Simulation scenario.}
  \item{seed}{Optional random seed.}
}
\value{
A list containing the simulated responses, true object scores, true category
quantifications, set sizes, and category counts.
}

