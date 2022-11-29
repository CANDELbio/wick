
#' Group data by assay and measurement set
#'
#' This function groups a measurements query result by
#' assay and measurement set
#'
#' @param tab The query result

#' @param measure.var The name of the column in \code{tab} representing
#'   the measurement value (e.g. \code{"count"})
#' @param target.var Optional. The name(s) of the column(s) in \code{tab} representing
#'   the target of the measurement (e.g. \code{"gene.hugo"} or \code{c("cell.pop", "epitope.id")})
#'   If not specified, the result will only contain the data in \code{measure.var}. Must be supplied
#'   if \code{cast.output} is specified
#' @param cast.output Whether the resulting data should be \code{cast} using \code{reshape2}
#'   or kept in \code{melt} format.
#' @param fun.aggregate Only used if \code{cast.output == TRUE}. Function to aggregate
#'   across multiple measurements with the same sample and target.
#'   Defaults to \code{NULL} and will emit a warning and
#'   use \code{length} if an aggregation function is needed but not supplied.
#'
#' @return Returns a named list of datasets, one for each assay and
#'   measurement set combination, or a single dataset if there
#'   is only one such combinations. Each dataset is either a \code{matrix},
#'   if \code{cast.output == FALSE} or a \code{data.frame} otherwise
group_by_assay_meas_set <- function(tab,
                                    measure.var,
                                    target.var = NULL,
                                    cast.output = TRUE,
                                    fun.aggregate = NULL) {
    cast.formula <- sprintf("sample.id ~ %s", paste(target.var, collapse = "+"))

    tab.list <- plyr::dlply(tab, ~assay.name * measurement.set.name, function(x) {
        if(is.null(target.var)) {
            ret <- as.matrix(x[, c(measure.var)])
            row.names(ret) <- x$sample.id
        } else if(cast.output)
            ret <- reshape2::acast(data = x,
                                   formula = as.formula(cast.formula),
                                   value.var = measure.var,
                                   fun.aggregate = fun.aggregate)
        else {
            ret <- x[, setdiff(colnames(x), c("assay.name", "set.name"))]
            ret <- data.frame(ret, check.names = FALSE, stringsAsFactors = FALSE)
        }
        return(ret)
    })

    # Remove plyr-specific information
    tab.list <- lapply(tab.list, function(x) {
        if(is.matrix(x))
            as.matrix(x)
        else
            as.data.frame(x, check.names = FALSE)
    })

    if(length(tab.list) == 1)
        tab.list <- tab.list[[1]]
    return(tab.list)

}


#' Select single assay and measurement set from a list of measurement query results
#'
#' @param tab A \code{data.frame} returned from a "get_all_X_measurements" query function
#' @param assay.name The name of the assay. Required if there are multiple assays for this data type in this dataset
#' @param measurement.set.name The name of the measurement set. Required if there are multiple measurement sets
#'   and / or multiple assays in this dataset
#'
#'
#' @return Returns a \code{data.frame} of data corresponding to a specific
#'   assay and measurement set combination, or \code{tab} if they are
#'   both \code{NULL}
#' @export
select_assay_exp <- function(tab,
                             assay.name = NULL,
                             measurement.set.name = NULL){

    if(is.list(tab)) {
        if(is.null(assay.name) || is.null(measurement.set.name))
            stop("Multiple assays / measurement sets found. Please specify assay.name and measurement.set.name")

        if (assay.name == "all" && measurement.set.name == "all"){
            message("Combining across all assays and measurement sets.
                    Note that when duplicate measurements appear, median will be reported.")
            tab <- do.call(rbind, lapply(FUN = reshape2::melt,
                                         X = tab))
            tab <- reshape2::acast(data = tab,
                                   formula = Var1 ~ Var2,
                                   fun.aggregate = function(x) median(x, na.rm = TRUE))
        } else {
            list.name <- paste(assay.name, measurement.set.name, sep = ".")
            if (!(list.name %in% names(tab)))
                stop(sprintf("Assay and measurement set not found.
                             Possible options for assay.meas-set are: %s",
                             names(tab)))
            tab <- tab[[list.name]]
        }
    } else if (!(is.null(assay.name) & is.null(measurement.set.name)))
        message("Warning: Did not detect multiple assays or measurement sets.
                Ignoring inputs assay.name and measurement.set.name")

    return(tab)
}



#' Select specific targets in a measurement table
#'
#' @param tab A matrix of measurements, with columns representing measurement
#'   targets (e.g. genes, cell-populations etc.)
#' @param targets.to.include Vector of strings representing specific targets
#'   (column names in \code{tab}) to include. If \code{NULL} (default), includes
#'   all targets.
#' @param targets.to.exclude Vector of strings representing specific targets
#'   (column names in \code{tab}) to exclude. If \code{NULL} (default), includes
#'   all targets (does not exclude any).
#'
#' @return Matrix of measurements with specified columns included and excluded.
#' @export
select_targets <- function(tab,
                           targets.to.include = NULL,
                           targets.to.exclude = NULL){

    if (!is.null(targets.to.include)){
        tab <- tab[, dimnames(tab)[[2]] %in% targets.to.include]
    }

    if (!is.null(targets.to.exclude)){
        tab <- tab[, !(dimnames(tab)[[2]] %in% targets.to.exclude)]
    }

    return(tab)
}


#' Ensure that data is molten
#'
#' This function ensures that data returned from e.g. \code{\link{get_all_measurements}}
#' is in molten format
#'
#' @param meas.table Input \code{data.frame} as returned for instance from
#'   \code{\link{get_all_measurements}}
#' @param variable.name Only used if the data is not molten already. This will be passed
#'   to \code{reshape2::melt} and represents the name of variable used to store measured variable names
#'
#' @return Returns the input data in molten format (or unchanged if it was molten already)
ensure_molten <- function(meas.table, variable.name = "variable") {
    if(is.matrix(meas.table) || is.null(meas.table$sample.id)) {#The data is not in molten format
        meas.table <- data.frame(meas.table, stringsAsFactors = FALSE, check.names = FALSE)
        meas.table$sample.id <- row.names(meas.table)
        meas.table <- reshape2::melt(meas.table, id.vars = "sample.id", variable.name = variable.name)
    }
    return(meas.table)
}

#' Add subject context to a table of data
#'
#'
#' @inheritParams add_sample_context
#' @export
#'
add_subject_context <- function(meas.table,
                                dataset.name,
                                include.response.data = TRUE, ...) {
    if(!("subject.id" %in% names(meas.table)))
        stop("subject.id column missing from meas.table, cannot add subject context")

    subjects <- get_all_subjects(dataset.name = dataset.name,
                                 include.response.data = include.response.data,
                                 ...)
    meas.table <- merge(meas.table, subjects, by = "subject.id")
    return(meas.table)
}


#' Add tcr context to a table of data with a tcr.id column
#'
#'
#' @param meas.table The table w/tcr.id that needs tcr data joined on to it.
#' @param dataset.name The name of the dataset which contains the tcrs
#' @param assay.name The name of the assay which contains the tcrs.
#' @param measurement.set.name The name of the measurement set which contains
#'                             the tcrs.
#' @param ... Additional params, such as \code{database.name}, are passed to
#'   \code{get_all_tcrs}
#' @export
#'
add_tcr_context <- function(meas.table,
                            dataset.name,
                            assay.name,
                            measurement.set.name,
                            ...) {
    if(!("tcr.id" %in% names(meas.table)))
        stop("tcr.id column missing from meas.table, cannot add tcr context")

    tcrs <- get_all_tcrs(dataset.name, assay.name, measurement.set.name, ...)
    meas.table <- merge(meas.table, tcrs, by = "tcr.id", all.x = TRUE)
    return(meas.table)
}


#' Add sample, treatment, and timepoint context to a table of measurements
#'
#' @inheritParams get_all_samples
#' @param meas.table A \code{data.frame}, or a \code{matrix} which will be
#'   converted to a \code{data.frame}, containing measurements,
#'   as returned by functions like \code{get_all_measurements}
#' @param include.subject.context Whether to also include available information
#'   about the subjects, beyond the subject.id
#' @param ... Additional params, such as \code{database.name}, are passed to
#'   \code{get_all_samples}
#' @return Returns a \code{data.frame} containing the measurements in melted
#'   format along with the sample, timepoint, and treatment regimen contextual
#'   information.
#' @export
add_sample_context <- function(meas.table,
                               dataset.name,
                               include.subject.context = TRUE,
                               include.response.data = TRUE,
                               ...) {

    if ('subject.id' %in% names(meas.table) | sum(grepl('sample', setdiff(names(meas.table), 'sample.id'))) > 0)
        stop('Sample info found in meas.table. Has sample context already been added?')

    samples <- get_all_samples(dataset.name, ...)
    samples$timepoint.treatment.regimen.name <- samples$treatment.regimen.name
    samples$treatment.regimen.name <- NULL

    if (!is.data.frame(meas.table) && !is.matrix(meas.table))
        stop('Measurement table must not be grouped by assay and measurement.set
              to add sample.context.
              Re-run with group.by.assay.meas.set = FALSE.')


    meas.table <- data.frame(meas.table,
                             check.names = FALSE,
                             stringsAsFactors = FALSE)

    meas.table <- ensure_molten(meas.table)

    meas.table <- merge(meas.table, samples, by = "sample.id")

    if(include.subject.context) {
        meas.table <- add_subject_context(meas.table = meas.table, dataset.name = dataset.name, ...)
    }

    meas.table$timepoint.id <- gsub("^(.)*/", "", meas.table$timepoint.id)

    return(meas.table)
}

#' Add variants context
#'
#' This function adds additional variant information to a table of measurements
#'
#' @param meas.table A \code{data.frame} of measurements as returned by
#'   \code{link{get_all_measurements}}
#' @param ... Additional parameters passed to \code{\link{get_all_variants}}
#'   e.g. \code{db.name}
#' @return Returns the input data with additional information about the variants
#'
#' @export
#'
add_variants_context <- function(meas.table, ...) {
    variants <- get_all_variants(...)
    meas.table <- ensure_molten(meas.table)


    if(!("variant.id" %in% colnames(meas.table))){
        stop("Column 'variant.id' not found. Ensure 'cast.output = FALSE' in get_all_measurements()")
    }

    ret <- merge(meas.table, variants, by = "variant.id", all.x = T)
    return(ret)
}

#' Add CNV context
#'
#' This function adds additional copy number variant information to a table of measurements
#'
#' @param meas.table A \code{data.frame} of measurements as returned by
#'   \code{link{get_all_measurements}}
#' @param ... Additional parameters passed to \code{\link{get_all_cnv}}
#'   e.g. \code{db.name}
#' @return Returns the input data with additional information about the CNVs
#'
#' @export
#'
add_cnv_context <- function(meas.table, ...) {
    cnv <- get_all_cnv(...)

    meas.table <- ensure_molten(meas.table)

    if(!("cnv.id" %in% colnames(meas.table))){
        stop("Column 'cnv.id' not found. Ensure 'cast.output = FALSE' in get_all_measurements()")
    }

    ret <- merge(meas.table, cnv, by = "cnv.id", all.x = T)
    return(ret)
}

#' Maps gene symbols
#'
#' This function maps gene symbols on the HGNC symbols used within the database.
#' Note that if you have multiple symbols to match, it is quicker to call this function once
#' with all the symbols, as opposed to calling it multiple times with the individual symbols
#'
#' @param all.genes Information about all the genes as returned by \link{get_all_genes}
#' @param gene.symbol Character vector of symbols to map
#' @param suppress.warnings Boolean suppressing warning messages indicating the number of successfully mapped genes
#'
#' @return Returns a named vector where the names are the original symbols and the corresponding values
#'   the ones they have been mapped to
#' @export
map_gene_symbols <- function(all.genes, gene.symbol, suppress.warnings = TRUE) {
    names.map <- as.list(all.genes$gene.hgnc.symbol)
    names(names.map) <- toupper(all.genes$gene.hgnc.symbol)

    for(col.name in c("gene.previous.hgnc.symbols", "gene.alias.hgnc.symbols")) {
        x <- all.genes[, col.name]
        names(x) <- all.genes[, "gene.hgnc.symbol"]
        x <- x[!sapply(x, is.null)]

        n.map <- sapply(names(x), function(n) {
            cbind(n, x[[n]])
        })

        n.map <- do.call(rbind, n.map)
        n.map.list <- as.list(n.map[, 1])
        names(n.map.list) <- toupper(n.map[, 2])

        names.map <- c(names.map, n.map.list)
    }

    ret <- names.map[toupper(gene.symbol)]
    names(ret) <- gene.symbol
    ret[sapply(ret, is.null)] <- NA
    ret <- unlist(ret)
    na.counts <- sum(is.na(ret))
    if (!suppress.warnings & na.counts > 0){
        warning(sprintf("%s of %s genes failed to map.", na.counts, length(ret)))
    }

    return(ret)

}

#' Flatten therapies information
#'
#' This function flattens and simplifies therapy information for an input table. This information
#' is typically in a nested \code{data.frame} format. This function concatenates that information
#' into a plain string
#'
#' @param tab Input \code{data.frame}. It must contain a column called \code{subject.therapies}
#' @return Returns the input data with the \code{subject.therapies} column simplified to
#'   a string concatenated with \code{;}
#'
#' @export
#'
#'
flatten_subject_therapies <- function(tab) {
    arms <- sapply(tab$subject.therapies, function(x) {
        if (!is.null(x)){
            x <- x[order(x$":therapy/order"), ]
            return(paste(x[, ":therapy/treatment-regimen"][, ":treatment-regimen/name"], collapse = ";"))
        } else
            return(NA)
    })
    return(unlist(arms))
}

#' Deduplicate otu taxonomic levels
#'
#' This function processes otu information to make sure that names at each taxonomic
#' level are unique. If that's the case, names are left unchanged, otherwise
#' the name is transformed to be a concatenation of the names of upstream taxonomic levels,
#' to guarantee uniqueness. This function is used internally by \code{aggregate_otu_measurements},
#' and it is unlikely that you would have to call it directly as an end-user
#'
#' @inheritParams aggregate_otu_measurements
#' @param sep Character separator used to concatenate names if necessary
#'
#' @return Returns the data in \code{all.otus} with the names at each taxonomic level processed
#'   as necessary
#'
#' @export
deduplicate_otu_taxonomic_levels <- function(all.otus, na.value, sep = "_") {
    level <- c("otu.kingdom", "otu.phylum", "otu.class", "otu.order", "otu.family", "otu.genus", "otu.species")

    ret <- all.otus

    for(i in 1:length(level)) {
        cur.level <- level[i]
        antecendent.levels <- level[1:i]
        w.na <- is.na(ret[, cur.level])
        ret[w.na, cur.level] <- na.value
        all.otus[w.na, cur.level] <- na.value
        unique.taxons <- unique(all.otus[, antecendent.levels, drop = FALSE])
        duplicated.at.this.level <- unique.taxons[duplicated(unique.taxons[, cur.level]), cur.level]
        w.duplicated <- ret[, cur.level] %in% duplicated.at.this.level
        if(any(w.duplicated))
            ret[w.duplicated, cur.level] <- apply(all.otus[w.duplicated, antecendent.levels, drop = FALSE], 1, paste, collapse = sep)

    }
    return(ret)
}


#' Aggregate otu measurements at different taxonomic levels
#'
#' This function aggregates OTU counts at different taxonomic level
#' to produce a table of per-taxa counts
#'
#' @param tab OTU measurements as returned by \link{get_all_measurements}
#' @param all.otus A table with OTU information as returned by \link{get_all_otus}
#' @param cast.ouput If specified, the aggregated data is cast into a matrix, otherwise it will
#'   be returned molten
#' @param normalize Whether OTU read counts should be normalized to the total read count
#'   in each sample
#' @param na.value  Name to use for OTU's that have an \code{NA} at a specific taxonomic level
#'
#' @return Returns a named list, with elements corresponding to the data aggregated
#'   at the corresponding taxonomic level. Note that in the aggregated data, the names
#'   of some taxa will appear as the concatenation of the names of upstream taxonomic levels.
#'   This is necessary because, in some cases, the name at each level is not enough to uniquely identify
#'   a taxonomic group (e.g. there could be multiple families with the same name in
#'   different orders)
#'
#' @export
aggregate_otu_measurements <- function(tab, all.otus, cast.output = FALSE, normalize = TRUE,
                                       na.value = "Unclassified") {
    all.otus <- deduplicate_otu_taxonomic_levels(all.otus, na.value = na.value)
    level <- c("otu.kingdom", "otu.phylum", "otu.class", "otu.order", "otu.family", "otu.genus", "otu.species")
    tab <- ensure_molten(tab, variable.name = "otu.id")
    tab <- tab[!is.na(tab$value), ] # This is possible if the data was previously cast
    tab <- merge(tab, all.otus, by = "otu.id")
    all.otu.cols <- grep("otu\\.", names(tab), value = TRUE)
    ret <- list()

    for(i in 1:length(level)) {
        cur.level <- level[i]
        otu.levels <- level[1:i]
        m <- tab[, c("sample.id", "value", intersect(otu.levels, all.otu.cols))]

        m <- plyr::ddply(m, c("sample.id", otu.levels), function(x) {
            sum(x$value)
        })

        names(m) <- gsub("V1", "value", names(m))

        if(normalize)
            m <- plyr::ddply(m, ~ sample.id, function(x) {
                x$value <- x$value / sum(x$value)
                return(x)}
            )

        if(cast.output) {
            m <- m[, c("sample.id", cur.level, "value")]
            m <- reshape2::dcast(m, as.formula(sprintf("sample.id ~ %s", cur.level)))
            row.names(m) <- m$sample.id
            m$sample.id <- NULL
            m <- as.matrix(m)
            m[is.na(m)] <- 0
        }

        ret[[cur.level]] <- m

    }
    return(ret)
}

#' Resolve db-idents in input data
#'
#' This function resolves columns containing db id's pointing to db-idents,
#' by substituting the values with the names of the idents that the id's point to.
#' This is useful for instance if you have a column of numeric id's in your data that represent
#' the id's of database enums, and you want to substitute them with the actual enums themselves
#'
#' @param tab The input data
#' @param col.names A character of vector of column names in \code{tab} that must be processed
#'
#' @return Returns \code{tab} with the specified columns processed
#'
#' @export
resolve_db_idents <- function(tab, col.names, ...) {
    all.idents <- get_all_db_idents(...)

    for(x in col.names) {
        v <- tab[, x]
        m <- match(v, all.idents$db.id)
        returned.nas <- is.na(m) != is.na(v)
        if(any(returned.nas)){
            unmatched <- v[returned.nas]
            stop(sprintf("%s unmatched idents in %s. First unmatched ident is %s",  sum(returned.nas),x, unmatched[1]))
        }
        tab[, x] <- all.idents[m, "db.ident"]
    }
    return(tab)
}


