% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wick.R
\name{get_query_server_url}
\alias{get_query_server_url}
\title{Return the query server URL associated with a specific database name}
\usage{
get_query_server_url(db.name)
}
\arguments{
\item{db.name}{If using the pret provision workflow, this is the database
name, ie "laceys-database". If using the pret request-db workflow, this is
the query url sent to you in your email notification, ie
"http://ec2-54-80-237-178.compute-1.amazonaws.com/query/cdel"}
}
\value{
Returns a server.url suitable to be passed to \code{datalogr::do_query}
}
\description{
Return the query server URL associated with a specific database name
}
