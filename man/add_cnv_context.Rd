% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/results.R
\name{add_cnv_context}
\alias{add_cnv_context}
\title{Add CNV context}
\usage{
add_cnv_context(meas.table, ...)
}
\arguments{
\item{meas.table}{A \code{data.frame} of measurements as returned by
\code{link{get_all_measurements}}}

\item{...}{Additional parameters passed to \code{\link{get_all_cnv}}
e.g. \code{db.name}}
}
\value{
Returns the input data with additional information about the CNVs
}
\description{
This function adds additional copy number variant information to a table of measurements
}
