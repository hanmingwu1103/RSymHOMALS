\name{load_symbolic_csv}
\alias{load_symbolic_csv}
\title{Load a Set-Valued Symbolic CSV File}
\description{
Parses a CSV file into the response structure used by \code{RSymHOMALS}.
Each symbolic variable is expected to appear as one CSV column whose name
begins with \code{Y1_}, \code{Y2_}, and so on.
}
\usage{
load_symbolic_csv(path)
}
\arguments{
  \item{path}{Path to a CSV file containing set-valued symbolic variables.}
}
\value{
A list containing the raw data frame, parsed set-valued responses, category
labels, variable names, object labels, and category counts.
}
