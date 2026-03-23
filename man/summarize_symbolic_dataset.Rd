\name{summarize_symbolic_dataset}
\alias{summarize_symbolic_dataset}
\title{Summarize a Set-Valued Symbolic Dataset}
\description{
Computes variable-level set-size summaries, label frequencies, and the most
high-cardinality objects for a parsed symbolic dataset.
}
\usage{
summarize_symbolic_dataset(dataset, top_n = 5L)
}
\arguments{
  \item{dataset}{A dataset object returned by \code{load_symbolic_csv()}.}
  \item{top_n}{Number of highest-cardinality objects to report.}
}
\value{
A list containing variable summaries, label-frequency tables, and a table of
the richest objects in the dataset.
}
