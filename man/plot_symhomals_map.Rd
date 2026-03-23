\name{plot_symhomals_map}
\alias{plot_symhomals_map}
\title{Plot a Symbolic HOMALS Map}
\description{
Creates a simple two-dimensional object-category map for a fitted Symbolic
HOMALS model.
}
\usage{
plot_symhomals_map(
  fit,
  variable = 1L,
  main = NULL,
  show_labels = FALSE,
  object_col = "grey40",
  category_col = "firebrick",
  object_pch = 16,
  category_pch = 17,
  ...
)
}
\arguments{
  \item{fit}{A fitted \code{symhomals_fit} object.}
  \item{variable}{Which variable's category points should be plotted.}
  \item{main}{Optional plot title.}
  \item{show_labels}{Logical flag for category labels.}
  \item{object_col}{Color used for object points.}
  \item{category_col}{Color used for category points.}
  \item{object_pch}{Plot symbol for objects.}
  \item{category_pch}{Plot symbol for categories.}
  \item{\dots}{Additional plotting arguments.}
}
\value{
No return value. The function is called for its side effect.
}

