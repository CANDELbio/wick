% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/queries_dataset.R
\name{get_all_timepoints}
\alias{get_all_timepoints}
\title{Return all timepoint IDs for a given dataset}
\usage{
get_all_timepoints(dataset.name, ...)
}
\arguments{
\item{...}{Additional arguments to be passed to \code{\link{do_query}} (typically
\code{db.name}, unless it has been set globally)}
}
\value{
Returns a \code{data.frame} containing information about all the
  timepoints in a dataset, including regimen name and relative order when
  applicable.
}
\description{
Return all timepoint IDs for a given dataset
}
