% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/results.R
\name{select_targets}
\alias{select_targets}
\title{Select specific targets in a measurement table}
\usage{
select_targets(tab, targets.to.include = NULL, targets.to.exclude = NULL)
}
\arguments{
\item{tab}{A matrix of measurements, with columns representing measurement
targets (e.g. genes, cell-populations etc.)}

\item{targets.to.include}{Vector of strings representing specific targets
(column names in \code{tab}) to include. If \code{NULL} (default), includes
all targets.}

\item{targets.to.exclude}{Vector of strings representing specific targets
(column names in \code{tab}) to exclude. If \code{NULL} (default), includes
all targets (does not exclude any).}
}
\value{
Matrix of measurements with specified columns included and excluded.
}
\description{
Select specific targets in a measurement table
}
