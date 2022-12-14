% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_functions.R
\name{make_event_timeline}
\alias{make_event_timeline}
\title{Make timeline of events}
\usage{
make_event_timeline(
  num_trials = 5,
  events = c("A", "B", "C"),
  onsets = list(1, round(rnorm(5, 50, 10)), 1),
  durations = list(25, 25, 199)
)
}
\arguments{
\item{num_trials}{integer, the number of trials}

\item{events}{character vector, names of each event}

\item{onsets}{list, containing an integer, or vector, defining onsets for each event}

\item{durations}{list, containing an integer, or vector, defining durations for each event}
}
\value{
dataframe
}
\description{
Make timeline of events
}
\examples{
make_event_timeline(
 num_trials = 5,
 events = c('A', 'B', 'C'),
 onsets = list(1,round(rnorm(5,50,10)),1),
 durations = list(25,25,199)
 )
}
