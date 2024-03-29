% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/queries_dataset.R
\name{get_dataset_summary}
\alias{get_dataset_summary}
\title{Return dataset summary}
\usage{
get_dataset_summary(dataset.name, ...)
}
\arguments{
\item{...}{Additional arguments to be passed to \code{\link{do_query}} (typically
\code{db.name}, unless it has been set globally)}
}
\value{
Returns a \code{data.frame} containing assay names, assay
  descriptions, and a nested \code{data.frame} of associated measurement-sets
  for each assay. Useful for looking up assay names in a dataset for further
  selection and query, as with \code{wick::get_all_measurements()}.
}
\description{
Return dataset summary
}
