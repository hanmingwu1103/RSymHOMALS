\name{sym_example_path}
\alias{sym_example_path}
\title{Locate Packaged Example Files}
\description{
Returns the path to the package's example files stored under
\code{inst/extdata}.
}
\usage{
sym_example_path(file = NULL)
}
\arguments{
  \item{file}{Optional file name under \code{inst/extdata}. If omitted, the
  directory path is returned.}
}
\value{
A character string giving the requested path.
}
