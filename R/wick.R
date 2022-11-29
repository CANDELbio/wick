#' wick
#'
#' @name wick
#'
#'
#' @importFrom datalogr d query find not not_join or pull where with args
NULL

pkg.env <- new.env()
pkg.env$query.server <- "http://localhost:8988"


#' Globally set the database name
#'
#' This function sets a database name globally for the entire package (so that you
#'   don't have to pass around a \code{db.name} argument everywhere)
#'
#' @param db.name If using the pret provision workflow, this is the database
#'   name, ie "laceys-database". If using the pret request-db workflow, this is
#'   the query url sent to you in your email notification, ie
#'   "http://ec2-12-12-123-678.compute-1.amazonaws.com/query/cdel"
#' @return Returns the previous \code{db.name} value, invisibly
#' @export
set_dbname <- function(db.name) {
    old <- pkg.env$db.name
    pkg.env$db.name <- db.name
    check_schema_compatibility(get_db_schema_version(db.name))
    return(invisible(old))
}

#' Globally set the auth token
#'
#' This function globally sets the authorization token. Users should use
#' \code{\link{login}}, instead of calling this directly
#'
#' @param token Character. The token
#'
#' @return Returns the previous token, invisibly
#'
#' @export
#'
set_token <- function(token) {
    old <- pkg.env$auth.token
    pkg.env$auth.token <- token
    return(invisible(old))
}

#' Sets the query server
#'
#' This function is for development purposes only. Sets the query server for this
#' session.
#'
#' @param query.server The query server address
#'
#' @return Returns the previous \code{query.server} value, invisibly
#'
#' @export
set_query_server <- function(query.server) {
    old <- pkg.env$query.server
    pkg.env$query.server <- query.server
    return(invisible(old))
}


#' Get the oauth url from candelabra
#'
#' @return Returns a list with two elements
#'   \itemize{
#'     \item oauth.url The oauth url
#'     \item redirect.url The redirect url
#'   }
get_candelabra_oauth_info <- function() {
    redirect.url <- httr::GET(url = paste(pkg.env$query.server, "cdel-oauth-viewer-url", sep = "/"))
    redirect.url <- httr::content(redirect.url, simplifyVector = TRUE)$"cdel-oauth-viewer-url"

    oauth.url <- httr::GET(url = paste(pkg.env$query.server, "client-oauth-url", sep = "/"))
    oauth.url <- httr::content(oauth.url, simplifyVector = TRUE)$"client-oauth-url"
    return(
        list(
            oauth.url = oauth.url,
            redirect.url = redirect.url
        )
    )
}

#' Login into the CANDEL system
#'
#' This function is used to login into the CANDEL system. It only needs
#' to be invoked once per session (unless the token expires, in which
#' case it needs to be called again, see below). Upon calling this function,
#' if a cached authorization token is not already available,
#' a browser window will open asking the user to authenticate into his/her
#' account. When this is done, a page will display an authorization code
#' that needs to be entered when prompted by the login function, thus
#' completing the login process
#'
#' @param email DEPRECATED. The user email
#' @param force By default this function reads a token cached in the \code{.candel}
#'   directory in the user home. If this token has expired, you can pass the
#'   \code{force} option to force a refresh
#'
#' @export
login <- function(email = NULL, force = FALSE) {
    pkg.env$auth.token <- ""
    return()
}

#' Return the query server URL associated with a specific database name
#'
#' @inheritParams set_dbname
#' @return Returns a server.url suitable to be passed to \code{datalogr::do_query}
#'
get_query_server_url <- function(db.name) {
    url <- NULL
    if (substr(db.name, 1, 7) == 'http://')
        url <- db.name
    else
        url <- paste(pkg.env$query.server, "query", db.name, sep = "/")
    return(url)
}

#' Ensure db.name has been set correctly
#'
#' @param db.name The name of the database, possibly \code{NULL}
#'
#' @return Returns \code{db.name} if not \code{NULL}. Otherwise returns the global
#'   \code{pkg.env$db.name} value
#'
#' @export
ensure_dbname <- function(db.name) {
    if(is.null(db.name)) {
        if(is.null(pkg.env$db.name))
            stop("Either provide a db.name or globally set the db name using set_dbname")
        db.name <- pkg.env$db.name
    } else {
        check_schema_compatibility(get_db_schema_version(db.name))
    }
    return(db.name)
}

#' Runs a query against the database
#'
#' @param db.name Optional. The name of the database against which the
#'   query is run. This can be set globally using \code{\link{set_dbname}}
#' @param ... Additional arguments. Please refer to the documentation of
#'   \code{datalogr::do_query}
#' @return Returns the query results. Refer to the documentation of
#'   \code{datalogr::do_query} for more details
#'
#' @export
do_query <- function(query, db.name = NULL, ...) {
    if(grepl("candelabra", pkg.env$query.server) && is.null(pkg.env$auth.token))
        stop("Please authenticate your session with the login function before issuing queries")
    db.name <- ensure_dbname(db.name)
    server.url <- get_query_server_url(db.name)
    ret <- datalogr::do_query(query = query, server.url = server.url, auth.token = pkg.env$auth.token, ...)

    if(length(ret) > 0 && datalogr::is_pull_query(query)) {
        # Remove synthetic id's
        cols.to.exclude <- c("measurement.id", "clinical.observation.id",
                             "drug.regimen.id", "therapy.id")
        ret <- ret[, !(names(ret) %in% cols.to.exclude), drop = FALSE]
    }
    return(ret)

}

#'
#'#' @return Returns the named list of schema versions that are compatible
#'           with each other.
get_schema_compatibility_overrides <- function() {
    # literal for now, we can make a separate file if this gets unwieldy
    return(list(`1.1.0` = c("1.2.0", "1.2.1"),
                `1.1.1` = c("1.2.0", "1.2.1"),
                `1.2.0` = c("1.1.0", "1.1.1"),
                `1.2.1` = c("1.1.0", "1.1.1")))
}

#' Check schema compatibility
#'
#' This function checks compatibility between the version of \code{wick} and the version
#' of the database schema, and issues a warning if they are not compatible
#'
#' @param schema.ver The schema version (as returned by \link{get_db_schema_version})
check_schema_compatibility <- function(schema.ver) {
    wick.ver <- packageVersion("wick")

    # if explicitly marked as compatible, return without further checks.
    overrides <- get_schema_compatibility_overrides()
    # paste as a shortcut to non-recursive unlist
    if (schema.ver %in% overrides[[paste(wick.ver)]])
        return(invisible(NULL))

    # default logic is based on differing major or minor version.
    wick.major.ver <- paste(strsplit(as.character(wick.ver), "\\.")[[1]][1:2], collapse = ".")
    schema.major.ver <- paste(strsplit(schema.ver, "\\.")[[1]][1:2], collapse = ".")

    if(schema.major.ver != wick.major.ver)
        warning(sprintf("The schema of the database (%s) and your installed wick package (%s) are incompatible
  some functions will likely not work. Consider installing a different version of wick
  or connecting to a different database", schema.ver, wick.ver))
    return(invisible(NULL))
}

#' Return the database schema version
#'
#' @param db.name The name of the database
#' @return Returns the schema version as a string
#'
#' @export
get_db_schema_version <- function(db.name) {
    if(grepl("candelabra", pkg.env$query.server)) {
        if(is.null(pkg.env$auth.token))
            stop("Please authenticate your session with the login function before issuing queries")
        else {
            url <- paste(pkg.env$query.server, "schema-version", db.name, sep = "/")
            x <- httr::RETRY("GET", url = url, encode = "raw", config = httr::add_headers(Authorization = paste("Token", pkg.env$auth.token)),
                             quiet = TRUE)
            content <- httr::content(x)

            if(httr::status_code(x) == 200)
                return(content$version)
            else if(httr::status_code(x) == 404) {
                # Support for old query api, remove eventually
                get_db_schema_version_via_query(db.name)
                # stop(sprintf("URL not found %s", url)) We should use this eventually
            } else if(httr::status_code(x) == 401) {
                stop("Authentication token is expired, please login again with force = TRUE")
            }
            else
                stop(content)

        }
    } else { # Support for non-candelabra stack, should be removed eventually
        get_db_schema_version_via_query(db.name)

    }

}


#' Return a matrix file as a data.frame.
#'
#' @param matrix.file The matrix file to retrieve.
#' @param db.name Optional: the name of the database
#' @return Returns the schema version as a string
#'
#' @export
get_matrix_file <- function(matrix_file, db.name = NULL) {
    if(grepl("candelabra", pkg.env$query.server) && is.null(pkg.env$auth.token))
        stop("Please authenticate your session with the login function before retrieving matrices.")
    db.name <- ensure_dbname(db.name)
    url <- paste(pkg.env$query.server, "matrix-url", db.name,
                 matrix_file, sep = "/")
    x <- httr::RETRY("GET", url = url, encode = "raw",
                     config = httr::add_headers(Authorization = paste("Token", pkg.env$auth.token)),
                     quiet = TRUE)
    content <- httr::content(x)

    if(httr::status_code(x) == 200) {
        s3_url <- content$`matrix-url`
        matrix.text <- httr::RETRY("GET", url = s3_url,
                                   times = 1000, quiet = T)
        result.tibble <- httr::content(matrix.text, as="parsed")
        return(as.data.frame(result.tibble))
        } else if(httr::status_code(x) == 404) {
        stop(sprintf("Specified matrix file not available in database: %s",
                     db.name))
    } else if(httr::status_code(x) == 401) {
        stop("Authentication token is expired, please login again with force = TRUE")
    }
    else
        stop(content)
}



#' Return the database schema version via a query
#'
#' @param db.name The name of the database
#' @return Returns the schema version as a string
#'
#' @export
get_db_schema_version_via_query <- function(db.name) {
    q <- query(
        find(?v),
        where(d(candel/schema, candel.schema/version, ?v))
    )

    # Needs to skip the ensure_dbname call flow to avoid infinite recursion
    server.url <- get_query_server_url(db.name)


    return(datalogr::do_query(query = q,
                              server.url = server.url,
                              auth.token = pkg.env$auth.token))
}




