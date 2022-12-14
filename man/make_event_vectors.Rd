% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_functions.R
\name{make_event_vectors}
\alias{make_event_vectors}
\title{Make Event Vectors}
\usage{
make_event_vectors(
  event_names = LETTERS[1:10],
  vector_length = 10,
  type = "field",
  sparsity = 0.1
)
}
\arguments{
\item{event_names}{character vector, names for each event}

\item{vector_length}{integer, number of 1s for each vector}

\item{type}{string, 'field' sets all values in each event vector to 1, which is used for the concatenation approach, 'riv' creates sparse random index vectors used in the superposition approach.}

\item{sparsity}{numeric, proportion of non-zero elements (1s, or -1s) in a random index vector}
}
\value{
matrix
}
\description{
Make Event Vectors
}
\examples{
make_event_vectors()
}
