\name{plot_symbolic_embedding}
\alias{plot_symbolic_embedding}
\title{Plot a Multi-Variable Symbolic Embedding}
\description{
Draws a two-dimensional Symbolic HOMALS or naive MCA map showing object scores
and category points across all variables.
}
\usage{
plot_symbolic_embedding(
  fit,
  main = NULL,
  show_category_labels = TRUE,
  category_label_mode = c("selected", "all", "none"),
  label_top_n = 2L,
  show_object_labels = FALSE,
  object_cex = 0.6,
  category_cex = 1,
  legend_cex = 0.8,
  ...
)
}
\arguments{
  \item{fit}{A fitted object returned by \code{fit_symhomals()} or the naive
  MCA comparison routine.}
  \item{main}{Optional plot title.}
  \item{show_category_labels}{Logical flag indicating whether category labels
  should be drawn.}
  \item{category_label_mode}{Whether to draw no category labels, all category
  labels, or only a selected subset. The default \code{"selected"} labels the
  most extreme category points within each variable.}
  \item{label_top_n}{Number of category labels to show per variable when
  \code{category_label_mode = "selected"}.}
  \item{show_object_labels}{Logical flag indicating whether object labels
  should be drawn.}
  \item{object_cex}{Point size for object scores.}
  \item{category_cex}{Point size for category scores.}
  \item{legend_cex}{Text size for the legend.}
  \item{\dots}{Additional graphical arguments passed to \code{plot()}.}
}
\value{
No return value. The function is called for its side effect.
}
