#' Query to return all measurements for a given measurement type
#'
#' @param dataset.name The name of the dataset
#' @param technology DEPRECATED - this argument is unused and will be
#' removed in a future release.
#' @param measurement.type The type of measurement, without the
#'   \code{:measurement} namespace portion, ie \code{"nanostring-count"} or
#'   \code{"percent-of-parent"}
#' @param include.targets Whether to include the targets of the measurement in
#'   the result. The types of targets vary according to the measurement type,
#'   for instance targets are genes for \code{measurement.type =
#'   "nanostring-count"}, cell populations for \code{measurement.type =
#'   "percent-of-leukocytes"} etc.
#' @param assay.name Optional, if provided restricts the results to measurements
#'   associated to the assay with the given name
#' @param measurement.set.name Optional, if provided restricts the results to
#'   measurements associated with the measurement set with the given name
#' @param sel.cell.pops Optional. Restrict the results to the specified cell population names
#' @param sel.epitopes Optional. Restrict the results to the specific epitope names
#' @return Returns a query suitable to get all the measurements for a given
#'   measurement type
#'
#' @export
all_measurements <- function(dataset.name,
                             measurement.type,
                             technology = NULL,
                             assay.name = NULL,
                             measurement.set.name = NULL,
                             sel.epitopes = NULL,
                             sel.cell.pops = NULL,
                             include.targets = TRUE) {

    if (!is.null(technology))
        message("Warning: the technology argument to get_all_measurements is deprecated.
                You should remove it from your function calls as
                it will be removed in a future release.")

    measurement.type <- paste(":measurement", measurement.type, sep = "/")

    gene.clause <- query(
        find(?gene-hgnc-symbol),
        where(
            d(?m, measurement/gene-product, ?gp),
            d(?gp, gene-product/gene, ?g),
            d(?g, gene/hgnc-symbol, ?gene-hgnc-symbol)
        )
    )

    cell.pop.clause <- query(
        find(?cell-population-name),
        where(
            d(?m, measurement/cell-population, ?c),
            d(?c, cell-population/name, ?cell-population-name)
        )
    )

    if(!is.null(sel.cell.pops)) {
        cell.pop.clause <- c_query(
            cell.pop.clause,
            query(
                where(
                    generate_or(!!sel.cell.pops, ?c, cell-population/name)
                )
            )
        )
    }

    variant.clause <- query(
        find(?variant-id),
        where(
            d(?m, measurement/variant, ?v),
            d(?v, variant/id, ?variant-id)
        )
    )

    cnv.clause <- query(
        find(?cnv-id),
        where(
            d(?m, measurement/cnv, ?c),
            d(?c, cnv/id, ?cnv-id)
        )
    )

    epitope.clause <- query(
        find(?epitope-id),
        where(
            d(?m, measurement/epitope, ?p),
            d(?p, epitope/id, ?epitope-id)
        )
    )

    if(!is.null(sel.epitopes)) {
        epitope.clause <- c_query(
            epitope.clause,
            query(
                where(
                    generate_or(!!sel.epitopes, ?p, epitope/id)
                )
            )
        )
    }

    tcr.clause <- query(
        find(?tcr-id),
        where(
            d(?m, measurement/tcr, ?t),
            d(?t, tcr/id, ?tcr-id)
        )
    )

    otu.clause <- query(
        find(?otu-id),
        where(
            d(?m, measurement/otu, ?o),
            d(?o, otu/id, ?otu-id)
        )
    )

    extra.clauses <- switch(measurement.type,
                            ':measurement/vaf' = variant.clause,
                            ':measurement/t-ref-count' = variant.clause,
                            ':measurement/t-alt-count' = variant.clause,
                            ':measurement/t-depth' = variant.clause,
                            ':measurement/n-ref-count' = variant.clause,
                            ':measurement/n-alt-count' = variant.clause,
                            ':measurement/n-depth' = variant.clause,
                            ':measurement/nanostring-count' = gene.clause,
                            ':measurement/cell-count' = cell.pop.clause,
                            ':measurement/percent-of-parent' = cell.pop.clause,
                            ':measurement/percent-of-nuclei' = cell.pop.clause,
                            ':measurement/percent-of-lymphocytes' = cell.pop.clause,
                            ':measurement/percent-of-leukocytes' = cell.pop.clause,
                            ':measurement/percent-of-live' = cell.pop.clause,
                            ':measurement/percent-of-singlets' = cell.pop.clause,
                            ':measurement/percent-of-total-cells' = cell.pop.clause,
                            ':measurement/tpm' = gene.clause,
                            ':measurement/fpkm' = gene.clause,
                            ':measurement/rpkm' = gene.clause,
                            ':measurement/array-log-intensity' = gene.clause,
                            ':measurement/array-log-ratio' = gene.clause,
                            ':measurement/protein-array-log-intensity' = epitope.clause,
                            ':measurement/rsem-raw-count' = gene.clause,
                            ':measurement/rsem-normalized-count' = gene.clause,
                            ':measurement/rsem-scaled-estimate' = gene.clause,
                            ':measurement/read-count' = gene.clause,
                            ':measurement/read-count-otu' = otu.clause,
                            ':measurement/read-count-otu-rarefied' = otu.clause,
                            ':measurement/segment-mean-lrr' = cnv.clause,
                            ':measurement/a-allele-cn' = cnv.clause,
                            ':measurement/b-allele-cn' = cnv.clause,
                            ':measurement/absolute-cn' = cnv.clause,
                            ':measurement/baf' = cnv.clause,
                            ':measurement/baf-n' = cnv.clause,
                            ':measurement/loh' = cnv.clause,
                            ':measurement/ng-mL' = epitope.clause,
                            ':measurement/pg-mL' = epitope.clause,
                            ':measurement/luminex-mfi' = epitope.clause,
                            ':measurement/olink-npx' = epitope.clause,
                            ':measurement/median-channel-value' = list(cell.pop.clause,
                                                                       epitope.clause),
                            ':measurement/tcr-frequency' = tcr.clause,
                            ':measurement/tcr-count' = tcr.clause)

    q <- query(
        find(?assay-name, ?measurement-set-name, ?sample-id, ?value),
        where(
            d(?d, dataset/name, ?dataset-name),
            d(?d, dataset/assays, ?a),
            d(?a, assay/name, ?assay-name),
            d(?a, assay/measurement-sets, ?e),
            d(?e, measurement-set/name, ?measurement-set-name),
            d(?e, measurement-set/measurements, ?m),
            d(?m, !!measurement.type, ?value)
        ),

        args(?dataset-name <- !!dataset.name)
    )

    if(include.targets) {
        if(is.null(names(extra.clauses))){ # Only way to distinguish whether it's a list of list
            for(i in 1:length(extra.clauses))
                q <- datalogr::c_query(q, extra.clauses[[i]])
        } else
            q <- datalogr::c_query(q, extra.clauses)
    }


    if(!is.null(sel.epitopes) ||
       !is.null(sel.cell.pops))
        q$query$where <- rev(q$query$where)


    q <- datalogr::c_query(q, query(
        where(
            d(?m, measurement/sample, ?s),
            d(?s, sample/id, ?sample-id)
        )
    ))


    if(!is.null(assay.name))
        q <- datalogr::c_query(q, query(
            args(?assay-name <- !!assay.name)
        ))

    if(!is.null(measurement.set.name))
        q <- datalogr::c_query(q, query(
            args(?measurement-set-name <- !!measurement.set.name)
        ))





    return(q)

}

#' Query to find a measurement matrix backing file. The backing file is
#' the key that needs to be passed to candelabra to retrieve a full
#' measurement matrix file.
#'
#' @param dataset.name The dataset that this measurement matrix is part of.
#' @param assay.name The assay this measurement matrix corresponds to.
#' @param measurement.set.name The measurement set that includes this
#'   measurement matrix.
#' @return Returns a query suitable to find a measurement matrix backing file.
#'
#' @export
measurement_matrices <- function(dataset.name,
                               assay.name,
                               measurement.set.name) {
    q <- query(
        find(?matrix-name, ?backing-file),
        where(
            d(?d, dataset/name, ?dataset-name),
            d(?d, dataset/assays, ?a),
            d(?a, assay/name, ?assay-name),
            d(?a, assay/measurement-sets, ?ms),
            d(?ms, measurement-set/name, ?ms-name),
            d(?ms, measurement-set/measurement-matrices, ?mat),
            d(?mat, measurement-matrix/name, ?matrix-name),
            d(?mat, measurement-matrix/backing-file, ?backing-file)
        ),
        args(?dataset-name <- !!dataset.name,
             ?assay-name <- !!assay.name,
             ?ms-name <- !!measurement.set.name)
    )
    return(q)
}

#' Query to find all measurement matrices and backing files.
#'
#' @return Returns a query suitable to find all measurement matrix backing
#' files. Joins each matrix to its assay, dataset, measurement set context.
#'
#' @export
all_measurement_matrices <- function() {
    q <- query(
        find(?matrix-name,
             ?backing-file,
             ?measurement-set-name,
             ?assay-name,
             ?dataset-name),
        where(
            d(?m, measurement-matrix/name, ?matrix-name),
            d(?m, measurement-matrix/backing-file, ?backing-file),
            d(?ms, measurement-set/measurement-matrices, ?m),
            d(?ms, measurement-set/name, ?measurement-set-name),
            d(?a, assay/measurement-sets, ?ms),
            d(?a, assay/name, ?assay-name),
            d(?d, dataset/assays, ?a),
            d(?d, dataset/name, ?dataset-name)
        )
    )
    return(q)
}


#' Query to find all measurement matrices and backing files.
#'
#' @inheritParams measurement_matrices
#' @param measurement.matrix.name The name of the measurement matrix to find
#' the backing file for.
#' @return Returns a query suitable to find all measurement matrix backing
#' files. Joins each matrix to its assay, dataset, measurement set context.
#'
#' @export
matrix_backing_file <- function(dataset.name,
                                assay.name,
                                measurement.set.name,
                                measurement.matrix.name,
                                ...) {
    q <- query(
        find(?backing-file),
        where(
            d(?m, measurement-matrix/name, ?matrix-name),
            d(?m, measurement-matrix/backing-file, ?backing-file),
            d(?ms, measurement-set/measurement-matrices, ?m),
            d(?ms, measurement-set/name, ?measurement-set-name),
            d(?a, assay/measurement-sets, ?ms),
            d(?a, assay/name, ?assay-name),
            d(?d, dataset/assays, ?a),
            d(?d, dataset/name, ?dataset-name)
        ),
        args(?dataset-name <- !!dataset.name,
             ?assay-name <- !!assay.name,
             ?ms-name <- !!measurement.set.name,
             ?matrix-name <- !!measurement.matrix.name)
    )
}

#' Get all the measurement matrix metadata for a given dataset, assay,
#' measurement set, and measurement matrix (name) context.
#'
#' @inheritParams measurement_matrices
#' @param ... Additional parameters passed to \code{do_query}
#' @return All measurement matrices in a database.
#' @export
get_measurement_matrices <- function(dataset.name,
                                     assay.name,
                                     measurement.set.name,
                                     ...) {
    q <- measurement_matrices(dataset.name, assay.name,
                              measurement.set.name)
    res <- do_query(q, ...)
    return(res)
}


#' Get the full measurement matrix for a given dataset, assay, measurement set,
#' and measurement matrix (name) context.
#'
#' @inheritParams matrix_backing_file
#' @param ... Additional parameters passed to \code{do_query}
#' @return All measurement matrices in a database.
#' @export
get_measurement_matrix <- function(dataset.name,
                                   assay.name,
                                   measurement.set.name,
                                   measurement.matrix.name,
                                   ...) {
    q <- matrix_backing_file(dataset.name,
                             assay.name,
                             measurement.set.name,
                             measurement.matrix.name)
    backing.file <- do_query(q, ...)
    if (is.null(backing.file))
        stop("No matrix file found for specified measurement matrix.")
    mx <- get_matrix_file(backing.file, ...)
    name.mapping <- list(':measurement-matrix/single-cells' = 'single.cell.id',
                         ':measurement-matrix/samples' = 'sample.id',
                         ':measurement-matrix/cell-populations' = 'cell.population.name',
                         ':measurement-matrix/epitopes' = 'epitope.id',
                         ':measurement-matrix/gene-products' = 'gene.hgnc.symbol')
    for (i in names(name.mapping)){
        names(mx)[names(mx) == i] <- name.mapping[i]
    }
    return(mx)
}

#' Retrieve a table of measurement matrices in a database, with information on
#' dataset, assay, and measurement-set context for retrieving the full matrix
#' with \code{get_measurement_matrix}.
#'
#' @param ... Additional parameters passed to \code{do_query}
#' @return All measurement matrices in a database.
#' @export
get_all_measurement_matrices <- function(...) {
    q <- all_measurement_matrices()
    res <- do_query(q, ...)
    return(res)
}


#' Get all measurements for a given measurement type
#'
#' @inheritParams all_measurements
#' @inheritParams group_by_assay_meas_set
#' @param include.sample.context Whether to include the sample context in the results
#'   using \code{link{get_sample_context}}. Note that if this option is specified
#'   the output will be molten, irrespective of the value of \code{cast.output}
#' @param ... Additional parameters passed to \code{do_query}
#' @inherit group_by_assay_meas_set return
#' @export
get_all_measurements <- function(dataset.name,
                                 measurement.type,
                                 assay.name = NULL,
                                 technology = NULL,
                                 measurement.set.name = NULL,
                                 include.targets = TRUE,
                                 cast.output = TRUE,
                                 fun.aggregate = NULL,
                                 include.sample.context = FALSE,
                                 group.by.assay.meas.set = TRUE,
                                 sel.epitopes = NULL,
                                 sel.cell.pops = NULL,
                                 ...) {

    if (!is.null(technology))
        message("Warning: the technology argument to get_all_measurements is deprecated.
                You should remove it from your function calls as
                it will be removed in a future release.")

    ds <- get_dataset_summary(dataset.name, ...)
    measurement.sets <- (Reduce(rbind, ds$assay.measurement.sets))$`:measurement-set/name`

    if (!is.null(measurement.set.name) &&
        !(measurement.set.name %in% measurement.sets)) {
       stop(sprintf("The measurement set '%s' does not exist!", measurement.set.name))
    }

    if (!is.null(assay.name) &&
        !(assay.name %in% ds$assay.name)) {
        stop(sprintf("The assay name '%s' does not exist!", assay.name))
    }

    ret <- do_query(all_measurements(dataset.name = dataset.name,
                                     measurement.type = measurement.type,
                                     assay.name = assay.name,
                                     measurement.set.name = measurement.set.name,
                                     include.targets = include.targets,
                                     sel.epitopes = sel.epitopes,
                                     sel.cell.pops = sel.cell.pops),
                    optimize = FALSE, ...)

    target <- NULL
    if(include.targets) {
        s <- setdiff(colnames(ret), c("assay.name", "measurement.set.name", "sample.id", "value"))
        if(length(s) != 0)
            target <- s
    }

    if (group.by.assay.meas.set)
        ret <- group_by_assay_meas_set(tab = ret,
                                       measure.var = "value",
                                       target.var = target,
                                       cast.output = cast.output,
                                       fun.aggregate = fun.aggregate)
    else if (cast.output)
        warning('Warning: Returning measurements in melted format.
                You cannot have cast.output = TRUE and group.by.assay.meas.set = FALSE.')

    if(include.sample.context) {
        if(cast.output)
            warning("Including sample context. Output will be molten even though cast.output is TRUE")
        if(inherits(ret, "list"))
            stop("If you want to include the sample context please select a single assay and measurment set
                 or set group.by.assay.meas.set = FALSE.")
        ret <- add_sample_context(ret, dataset.name, ...)

    }

    return(ret)
}




#' Query for all subjects
#'
#' @param dataset.name The name of the dataset
#' @return Returns a query suitable for getting all subjects
#'   for a given dataset
#' @export
all_subjects <- function(dataset.name) {
    q <- query(
        find(pull(?s, c(.,
                        {subject/sex = c(db/ident)},
                        {subject/race = c(db/ident)},
                        {subject/ethnicity = c(db/ident)},
                        {subject/meddra-disease = c(meddra-disease/preferred-name)},
                        {subject/disease-stage = c(db/ident)},
                        {subject/smoker = c(db/ident)},
                        {subject/cause-of-death = c(db/ident)},
                        {subject/therapies = c(therapy/order,
                                               {therapy/treatment-regimen = c(treatment-regimen/name)})}))),
        where(
            d(?d, dataset/name, ?dataset-name),
            d(?d, dataset/subjects, ?s)),
        args(?dataset-name <- !!dataset.name)
    )

    return(q)
}

#' Return all subject information for a given dataset
#'
#' This function gets information on all subjects in a dataset, including
#' demographics, treatment regimen, and (optionally) response information. If
#' response data is included, this function will error if multiple entities
#' exist for clinical-observations where we expect a single value, ie multiple
#' pfs or os values.
#'
#' @inheritParams get_all_gene_symbols
#' @param include.response.data Boolean. If TRUE (default), also pulls response
#'   data for all subjects, including PFS, OS, and Best Overall Response (bor).
#'   Will error if it finds more than one PFS or OS clinical-observation entity
#'   per subject. Current calculates BOR rather than pulling it from the
#'   database.
#'
#' @importFrom dplyr %>%
#' @export
get_all_subjects <- function(dataset.name,
                             include.response.data = TRUE,
                             ...) {

    subjects <- do_query(all_subjects(dataset.name), ...)

    # clean up subject-level info
    subjects$treatment.regimen.name <- flatten_subject_therapies(subjects)
    if ('subject.sex' %in% names(subjects))
        subjects$subject.sex <- sub(':subject.sex/', '', subjects$subject.sex)
    if ('subject.race' %in% names(subjects))
        subjects$subject.race <- sub(':subject.race/', '', subjects$subject.race)
    if ('subject.ethnicity' %in% names(subjects))
        subjects$subject.ethnicity <- sub(':subject.ethnicity/', '', subjects$subject.ethnicity)
    if ('subject.smoker' %in% names(subjects))
        subjects$subject.smoker <- sub(':subject.smoker/', '', subjects$subject.smoker)
    if ('subject.cause.of.death' %in% names(subjects))
        subjects$subject.cause.of.death <- sub(':subject.cause-of-death/', '', subjects$subject.cause.of.death)

    if (include.response.data){
        bor.tab <- get_all_clinical_observations(dataset.name = dataset.name,
                                                    obs.type = 'bor',
                                                    ...)

        if (nrow(bor.tab) == 0){
            recist.tab <- get_all_clinical_observations(dataset.name = dataset.name,
                                                        obs.type = 'recist',
                                                        ...)
            if (nrow(recist.tab)){
                message('BOR not found in the database, but RECIST found. Calculating BOR.')

                recist.tab$recist <- gsub(':clinical.observation.recist/', '',
                                          recist.tab$recist)


                bor.tab <- recist.tab %>%
                    dplyr::group_by(subject.id) %>%
                    dplyr::summarise(bor = dplyr::case_when(sum(recist == 'CR') > 0 ~ 'CR',
                                                            sum(recist == 'PR') > 0 ~ 'PR',
                                                            sum(recist == 'SD') > 0 ~ 'SD',
                                                            sum(recist == 'PD') > 0 ~ 'PD',
                                                            TRUE ~ 'Unknown')) %>%
                    data.frame(.,
                               stringsAsFactors = FALSE,
                               check.names = FALSE)
            }
        }

        bor.tab$bor <- gsub(':clinical.observation.recist/', '', bor.tab$bor)
        bor.tab$timepoint.id <- NULL
        if (length(unique(bor.tab$subject.id)) < nrow(bor.tab))
            stop('Error - there is more than one BOR value per subject in the database.')

        # get PFS and make sure there's only one measurement per patient
        pfs.tab <- get_all_clinical_observations(dataset.name = dataset.name,
                                                 obs.type = 'pfs',
                                                 ...)
        if (length(unique(pfs.tab$subject.id)) < nrow(pfs.tab))
            stop('Error - there is more than one PFS value per subject in the database.')
        pfs.tab$timepoint.id <- NULL
        pfs.tab$pfs <- as.numeric(pfs.tab$pfs)

        # get OS and make sure there's only one measurement per patient
        os.tab <- get_all_clinical_observations(dataset.name = dataset.name,
                                                obs.type = 'os',
                                                ...)
        if (length(unique(os.tab$subject.id)) < nrow(os.tab))
            stop('Error - there is more than one OS value per subject in the database.')
        os.tab$timepoint.id <- NULL
        os.tab$os <- as.numeric(os.tab$os)

        if (length(subjects) > 0)
            subjects <- merge(subjects,
                                  bor.tab, by = 'subject.id',
                                  all.x = TRUE) %>%
                merge(.,
                      pfs.tab, by = 'subject.id',
                      all.x = TRUE) %>%
                merge(.,
                      os.tab, by = 'subject.id',
                      all.x = TRUE)
    }

    return(subjects)
}


#' Query to get all cell populations
#'
#' @param dataset.name The name of the dataset
#' @param assay.name The assay name
#' @param measurement.set.name The name of the measurement set
#'
#' @return Returns a query suitable for getting all the cell populations
#' @export
all_cell_populations <- function(dataset.name, assay.name, measurement.set.name) {
    query(
        find(pull(?p, c(.,
                        {cell-population/cell-type = c(cell-type/co-name)},
                        {cell-population/negative-markers = c(epitope/id)},
                        {cell-population/positive-markers = c(epitope/id)}))),
        where(
            d(?d, dataset/name, ?dataset-name),
            d(?d, dataset/assays, ?a),
            d(?a, assay/name, ?assay-name),
            d(?a, assay/measurement-sets, ?e),
            d(?e, measurement-set/name, ?measurement-set-name),
            d(?e, measurement-set/cell-populations, ?p)
        ),
        args(?dataset-name <- !!dataset.name,
             ?assay-name <- !!assay.name,
             ?measurement-set-name <- !!measurement.set.name
        )
    )
}

#' Get all cell populations
#' @inheritParams all_cell_populations
#' @param ... Additional parameters passed to \code{do_query}
#' @return Returns a character vector with all the cell populations names
#' @export
get_all_cell_populations <- function(dataset.name, assay.name,
                                     measurement.set.name, ...) {
    return(do_query(all_cell_populations(dataset.name, assay.name, measurement.set.name), ...))
}

#' Query to get all tcrs
#'
#' @inheritParams all_cell_populations
#'
#' @return Returns a query suitable for getting all the cell populations
#' @export
all_tcrs <- function(dataset.name, assay.name, measurement.set.name) {
    query(
        find(pull(?tcr, c(.))),
        where(
            d(?d, dataset/name, ?dataset-name),
            d(?d, dataset/assays, ?a),
            d(?a, assay/name, ?assay-name),
            d(?a, assay/measurement-sets, ?e),
            d(?e, measurement-set/name, ?measurement-set-name),
            d(?e, measurement-set/tcrs, ?tcr)
        ),
        args(?dataset-name <- !!dataset.name,
             ?assay-name <- !!assay.name,
             ?measurement-set-name <- !!measurement.set.name
        )
    )
}


#' Get all TCR data.
#' @inheritParams all_cell_populations
#' @param ... Additional parameters passed to \code{do_query}
#' @return Returns a character vector with all the tcr names
#' @export
get_all_tcrs <- function(dataset.name, assay.name, measurement.set.name, ...) {
    return(do_query(all_tcrs(dataset.name,
                             assay.name,
                             measurement.set.name),
                    ...))
}


#' Query to get all OTU's
#'
#' @inheritParams all_cell_populations
#'
#' @return Returns a query suitable for getting all OTU's for a specific measurement set
#' @export
all_otus <- function(dataset.name, assay.name, measurement.set.name) {
    query(
        find(pull(?o, c(.))),
        where(
            d(?d, dataset/name, ?dataset-name),
            d(?d, dataset/assays, ?a),
            d(?a, assay/name, ?assay-name),
            d(?a, assay/measurement-sets, ?e),
            d(?e, measurement-set/name, ?measurement-set-name),
            d(?e, measurement-set/otus, ?o)
        ),
        args(?dataset-name <- !!dataset.name,
             ?assay-name <- !!assay.name,
             ?measurement-set-name <- !!measurement.set.name
        )
    )
}

#' Get all OTU's
#'
#' Get all OTU's for a specific measurement set
#'
#' @inheritParams get_all_cell_populations
#' @export
get_all_otus <- function(dataset.name, assay.name,
                         measurement.set.name, ...) {
    return(do_query(all_otus(dataset.name = dataset.name,
                             assay.name = assay.name,
                             measurement.set.name = measurement.set.name), ...))
}


#' Query for all samples
#' @inheritParams all_subjects
#' @return Returns a query suitable for getting all samples
#'   for a given dataset
#' @export
all_samples <- function(dataset.name) {

    q <- query(
        find(pull(?s, c(.,
                        {sample/specimen = c(db/ident)},
                        {sample/timepoint = c(timepoint/id,
                                              {timepoint/treatment-regimen = c(treatment-regimen/name)})},
                        {sample/study-day = c(study-day/id)},
                        {sample/subject = c(subject/id)},
                        {sample/container = c(db/ident)},
                        {sample/gdc-anatomic-site = c(gdc-anatomic-site/name)}))),
        where(d(?d, dataset/name, ?dataset-name),
              d(?d, dataset/samples, ?s)),
        args(?dataset-name <- !!dataset.name)
    )

    return(q)
}

#' Return all samples information for a given dataset
#'
#' @inheritParams get_all_gene_symbols
#'
#' @return Returns a \code{data.frame} containing information about all the samples
#'   in the dataset, joined with information about the subjects they originate
#'   from and treatment regimens associated with those subjects
#' @export
get_all_samples <- function(dataset.name, ...) {
    return(do_query(all_samples(dataset.name), ...))
}



#' For a given datatset, return all measurement sets and the samples and
#' timepoints they measure for a given dataset
#'
#' @inheritParams get_all_gene_symbols
#'
#' @return Returns a \code{data.frame} containing information about all the
#'   measurement sets in the dataset, their asssays, and the samples they
#'   measure, joined with information about the subjects they originate from and
#'   treatment regimens associated with those subjects
#' @export
get_all_meas_set_samples <- function(dataset.name, ...) {

    result <- try(do_query(wick::all_meas_set_samples(dataset.name = dataset.name), ...),
                  silent = TRUE)


    if (class(result) == 'try-error'){
        message('Dataset seems too large to profile in a single query.
                Iterating over samples, this may take awhile!')

        samples <- wick::get_all_samples(dataset.name = dataset.name, ...)

        result <- list()

        if (length(samples) > 0){

            for (s in samples$sample.id){
                result[[s]] <- do_query(all_meas_set_single_sample(dataset.name = dataset.name,
                                                            sample.id = s), ...)
            }

            result <- do.call(rbind, result)
        }
    }

    return(result)
}

#' Query for all measurement sets and assays for one specific sample
#' @inheritParams all_subjects
#' @return Returns a query suitable for getting all the above data for a given
#'   dataset and sample ID
#' @export
all_meas_set_single_sample <- function(dataset.name, sample.id) {

    q <- query(find(?subject-id,
                    ?assay-name,
                    ?measurement-set-name,
                    ?timepoint-id),
               where(
                   d(?sm, sample/id, ?sample-id),
                   d(?d, dataset/samples, ?sm),
                   d(?d, dataset/name, ?dataset-name),
                   d(?m, measurement/sample, ?sm),
                   d(?ms, measurement-set/measurements, ?m),
                   d(?ms, measurement-set/name, ?measurement-set-name),
                   d(?a, assay/measurement-sets, ?ms),
                   d(?a, assay/name, ?assay-name),
                   d(?sm, sample/subject, ?s),
                   d(?s, subject/id, ?subject-id),
                   d(?sm, sample/timepoint, ?tp),
                   d(?tp, timepoint/id, ?timepoint-id)),
               args(?dataset-name <- !!dataset.name,
                    ?sample-id <- !!sample.id))

    return(q)
}



#' Query for all measurement sets and assays for all samples in a dataset
#' @inheritParams all_subjects
#' @return Returns a query suitable for getting all the above data for a given
#'   dataset
#' @export
all_meas_set_samples <- function(dataset.name) {

    q <- query(find(?subject-id,
                    ?assay-name,
                    ?measurement-set-name,
                    ?timepoint-id),
               where(
                   d(?d, dataset/name, ?dataset-name),
                   d(?d, dataset/subjects, ?s),
                   d(?s, subject/id, ?subject-id),
                   d(?sm, sample/subject, ?s),
                   d(?m, measurement/sample, ?sm),
                   d(?ms, measurement-set/measurements, ?m),
                   d(?ms, measurement-set/name, ?measurement-set-name),
                   d(?a, assay/measurement-sets, ?ms),
                   d(?a, assay/name, ?assay-name),
                   d(?sm, sample/timepoint, ?tp),
                   d(?tp, timepoint/id, ?timepoint-id)),
               args(?dataset-name <- !!dataset.name))

    return(q)
}

#' For a given measurement-set, return all measurement types
#' of measurements in that set
#'
#' @inheritParams all_cell_populations
#' @param simple.strings reformats the datomic entity types to strings
#' containing just the name of the measurement type.
#'
#' @return Returns a \code{vector} containing all the measurement types
#' of measurements in the given measurement-set.
#' @export
get_all_meas_set_meas_types <- function(dataset.name,
                                        measurement.set.name,
                                        simple.strings = TRUE,
                                        ...){
    # target.types are filtered out from get_all_meas_set_measurement_attributes
    target.types <- c(":measurement/sample",
                      ":measurement/atac-peak",
                      ":measurement/cell-population",
                      ":measurement/cnv",
                      ":measurement/epitope",
                      ":measurement/gene-product",
                      ":measurement/nanostring-signature",
                      ":measurement/tcr",
                      ":measurement/variant")
    meas_attrs <- get_all_meas_set_measurement_attributes(measurement.set.name = measurement.set.name,
                                                          dataset.name = dataset.name,
                                                          ...)
    meas_types <- meas_attrs[!(meas_attrs %in% target.types)]
    if(simple.strings){
        meas_types <- unlist(lapply(meas_types,function(x){stringr::str_remove(x, ":measurement/")}))}
    return(meas_types)
}

#' For a given measurement-set, return all attributes
#' of measurements in that set
#'
#' @inheritParams all_cell_populations
#'
#' @return Returns a \code{vector} containing all the attributes
#' of measurements in the given measurement-set.
#' @export
get_all_meas_set_measurement_attributes <- function(dataset.name,
                                                    measurement.set.name,
                                                    ...){
    result <- do_query(wick:::all_meas_set_meas_types(dataset.name = dataset.name,
                                                      measurement.set.name = measurement.set.name),
                       ...)
    return(result)
}

#' Query for all measurement types in a measurement set
#' @inheritParams all_cell_populations
#' @return Returns a query suitable for getting all the above data for a given
#'   dataset
#' @export
all_meas_set_meas_types <- function(dataset.name, measurement.set.name){
    query(
        find(?attr),
        where(
            d(?d, dataset/name, ?dataset-name),
            d(?ms, measurement-set/name, ?measurement-set-name),
            d(?ms, measurement-set/measurements, ?m),
            d(?m, ?a),
            d(?a, db/ident, ?attr)
        ),
        args(?measurement-set-name <- !!measurement.set.name),
        args(?dataset-name <- !!dataset.name)
        )
}

#' Query for all timepoints in a dataset
#'
#' Query uses pull to include isolated timepoints with no regimen, and including
#' regimen name and relative order when applicable.
#'
#' @inheritParams all_subjects
#' @return Returns a query suitable for getting all timepoints in a given dataset
#' @export
all_timepoints <- function(dataset.name) {
    q <- query(
        find(pull(?t, c(timepoint/id,
                        {timepoint/type = c(db/ident)},
                        {timepoint/treatment-regimen = c(treatment-regimen/name)},
                        timepoint/relative-order,
                        timepoint/cycle,
                        timepoint/day,
                        timepoint/offset))),
        where (
            d(?d, dataset/name, ?dataset-name),
            d(?d, dataset/timepoints, ?t)),

        args(?dataset-name <- !!dataset.name)
    )
    return(q)
}

#' Return all timepoint IDs for a given dataset
#'
#' @inheritParams get_all_gene_symbols
#'
#' @return Returns a \code{data.frame} containing information about all the
#'   timepoints in a dataset, including regimen name and relative order when
#'   applicable.
#' @export
#'
get_all_timepoints <- function(dataset.name, ...) {
    return(do_query(all_timepoints(dataset.name), ...))
}





#' Query for all clinical observations
#'
#' @inheritParams all_subjects
#' @param clinical.observation.set.name The name of the clinical-observation-set. All the observations in that
#'   set will be pulled
#' @param obs.type A string. The type of the observation to return (e.g. \code{tumor-sum-diam}).
#'   Only used if clinical.observation.set.name is \code{NULL}, in which case all the observation
#'   of the given type will be returned.. Note that you should not include the
#'   namespace (i.e. \code{:clinical-observation} in the string)
#' @return Returns a query suitable for getting all clinical observations
#'   of a given type for a specified dataset
#'
#' @export
all_clinical_observations <- function(dataset.name, clinical.observation.set.name = NULL, obs.type = NULL) {
    q <- NULL

    if(!is.null(clinical.observation.set.name)) {
        q <- query(find(pull(?o, c(.,
                                   {clinical-observation/subject = c(subject/id)},
                                   {clinical-observation/timepoint = c(timepoint/id)},
                                   {clinical-observation/recist = c(db/ident)},
                                   {clinical-observation/pfs-reason = c(db/ident)},
                                   {clinical-observation/os-reason = c(db/ident)},
                                   {clinical-observation/dfi-reason = c(db/ident)},
                                   {clinical-observation/ttf-reason = c(db/ident)},
                                   {clinical-observation/disease-stage = c(db/ident)},
                                   {clinical-observation/metastasis-gdc-anatomic-sites = c(gdc-anatomic-site/name)},
                                   {clinical-observation/study-day = c(study-day/id)}))),
                   where(
                       d(?d, dataset/name, ?dataset.name),
                       d(?d, dataset/clinical-observation-sets, ?cos),
                       d(?cos, clinical-observation-set/name, ?clinical.observation.set.name),
                       d(?cos, clinical-observation-set/clinical-observations, ?o)
                   ),
                   args(?dataset.name <- !!dataset.name, ?clinical.observation.set.name <- !!clinical.observation.set.name)
        )
    } else if(!is.null(obs.type)) {
        obs.type <- paste(":clinical-observation", obs.type, sep = "/")
        q <- query(find(?subject-id, ?timepoint-id, ?value),
                   where(
                       d(?d, dataset/name, ?dataset.name),
                       d(?d, dataset/clinical-observation-sets, ?cos),
                       d(?cos, clinical-observation-set/clinical-observations, ?o),
                       d(?o, !!obs.type, ?value),
                       d(?o, clinical-observation/subject, ?s),
                       d(?o, clinical-observation/timepoint, ?t),
                       d(?s, subject/id, ?subject-id),
                       d(?t, timepoint/id, ?timepoint-id)
                   ),
                   args(?dataset.name <- !!dataset.name)
        )

    } else {
        stop("Please provide either clinical.observation.set.name or obs.type")
    }
    return(q)
}




#' Return all clinical observations for a given dataset
#'
#' @inheritParams all_clinical_observations
#'
#' @return Returns a \code{data.frame} containing all the clinical observation
#'   and their associated subjects and timepoints for a given dataset
#'
#' @export
get_all_clinical_observations <- function(dataset.name, clinical.observation.set.name = NULL, obs.type = NULL, ...) {
    ret <- do_query(all_clinical_observations(dataset.name = dataset.name,
                                              clinical.observation.set.name = clinical.observation.set.name,
                                              obs.type = obs.type),
                    ...)

    if(is.null(clinical.observation.set.name)) {
        if(obs.type %in% c("recist", "bor", "ir-recist", "dfi-reason", "pfs-reason",
                            "os-reason", "dfi-event", "pfs-event", "os-event"))
            ret <- resolve_db_idents(tab = ret, col.names = "value", ...)
        names(ret) <- make.names(gsub("value", obs.type, names(ret)))
    }
    return(ret)
}


#' Query for all adverse events
#'
#' @inheritParams all_subjects
#' @param clinical.observation.set.name The name of the clinical-observation-set. All the AEs in that
#'   set will be pulled
#' @return Returns a query suitable for getting all AEs
#'   for a specified dataset
#'
#' @export
all_adverse_events <- function(dataset.name, clinical.observation.set.name = NULL) {
    q <- NULL

    if(!is.null(clinical.observation.set.name)) {
        q <- query(find(pull(?o, c(.,
                                   {adverse-event/subject = c(subject/id)},
                                   {adverse-event/timepoint = c(timepoint/id)},
                                   {adverse-event/meddra-adverse-event = c(meddra-disease/preferred-name)},
                                   {adverse-event/ctcae-grade = c(db/ident)},
                                   {adverse-event/ae-causality = c(db/ident)},
                                   {adverse-event/study-day = c(study-day/id)}))),
                   where(
                       d(?d, dataset/name, ?dataset.name),
                       d(?d, dataset/clinical-observation-sets, ?cos),
                       d(?cos, clinical-observation-set/name, ?clinical.observation.set.name),
                       d(?cos, clinical-observation-set/adverse-events, ?o)
                   ),
                   args(?dataset.name <- !!dataset.name, ?clinical.observation.set.name <- !!clinical.observation.set.name)
        )
    } else {
        stop("Please provide clinical.observation.set.name")
    }
    return(q)
}

#' Return all adverse events for a given dataset
#'
#' @inheritParams all_clinical_observations
#'
#' @return Returns a \code{data.frame} containing all the clinical observation
#'   and their associated subjects and timepoints for a given dataset
#'
#' @export
get_all_adverse_events <- function(dataset.name, clinical.observation.set.name = NULL, ...) {
    ret <- do_query(all_adverse_events(dataset.name = dataset.name,
                                              clinical.observation.set.name = clinical.observation.set.name),
                    ...)
    return(ret)
}



#' Return timepoint information
#'
#' This function returns timepoint information for timepoints that are common to all
#' the treatment regimen. The timepoints are renamed by removing
#' the portion before the \code{/} separator. This function
#' is intended to generate simple timepoint information which is valid across
#' all treatment regimens for cases where all the treatment regimens have the
#' same timepoint structure. If the timepoint structure is not identical across all
#' treatment.regimens, this function will issue a warning. Note that in this case
#' the relative order of timepoints will be inconsistent across different treatment regimens.
#' THIS FUNCTIONALITY WILL BE REMOVED IN THE FUTURE AND ONLY COMMON TIMEPOINTS WILL BE RETURNED
#'
#' @inheritParams get_all_gene_symbols
#'
#' @return Returns a \code{data.frame} with timepoint information
#' @export
get_common_timepoints <- function(dataset.name, ...) {
    timepoints <- get_all_timepoints(dataset.name, ...)

    if(length(timepoints) == 0)
        stop('There are no timepoints in the specified dataset.')

    timepoints <- timepoints[!is.na(timepoints$treatment.regimen.name), ]

    timepoints$timepoint.id <- gsub("^(.)*/", "",  timepoints$timepoint.id)
    m <- as.matrix(table(timepoints[, c("treatment.regimen.name", "timepoint.id")]))

    excluded.timepoints <- names(which(colSums(m) != nrow(m)))

    #TODO: Remove this functionality and only return common timepoints

    if(length(excluded.timepoints) > 0) {
        warning(sprintf("The timepoints in dataset %s are not identical among
                        all treatment regimens, omitting non-shared timepoints: %s", dataset.name, paste(excluded.timepoints, sep = ",")))
        timepoints <- timepoints[!(timepoints$timepoint.id %in% excluded.timepoints), ]
    }

    timepoints <- timepoints[!duplicated(timepoints$timepoint.id), ]
    timepoints <- timepoints[order(timepoints$timepoint.relative.order), ]
    timepoints$treatment.regimen.name <- NULL
    return(timepoints)
}

#' Return a unified timeline for mulitple treatment regimens
#'
#' This function attempts to construct a unified timeline for multiple treatment regimens.
#' It assumes a situation where timepoints that are conceptually similar between multiple
#' treatment regimens have the same name, but potentially each treatment regimen has additional
#' timepoints interspersed that are not shared with other regimens. The unified timeline
#' is constructed by taking the union of all the timepoint ids, and, for each timepoint id,
#' selecting the maximum among all its relative orders across the treatment regimens
#'
#' @inheritParams get_all_gene_symbols
#' @param treatment.regimen.name optional character vector. If provided, only limit
#'   the timeline to the selected treatment regimens
#' @param ... Additional arguments passed to \code{get_all_timepoints}
#' @return Returns a \code{data.frame} with timepoint information
#'
#' @export
get_unified_timeline <- function(dataset.name, treatment.regimen.name = NULL, ...) {
    timepoints <- get_all_timepoints(dataset.name, ...)

    if(!is.null(treatment.regimen.name))
        timepoints <- timepoints[timepoints$treatment.regimen.name %in% treatment.regimen.name, ]

    timepoints$timepoint.id <- gsub("^(.)*/", "",  timepoints$timepoint.id)
    ret <- plyr::ddply(timepoints, ~timepoint.id, function(x) {
        ret <- x[1,, drop = FALSE]
        ret$timepoint.relative.order <- max(x$timepoint.relative.order)
        return(ret)
    })
    ret <- ret[order(ret$timepoint.relative.order), ]
    ret$treatment.regimen.name <- NULL

    return(ret)
}


#' Return dataset summary
#'
#'
#' @inheritParams get_all_gene_symbols
#'
#' @return Returns a \code{data.frame} containing assay names, assay
#'   descriptions, and a nested \code{data.frame} of associated measurement-sets
#'   for each assay. Useful for looking up assay names in a dataset for further
#'   selection and query, as with \code{wick::get_all_measurements()}.
#' @export
get_dataset_summary <- function(dataset.name, ...){
    return(do_query(dataset_summary(dataset.name = dataset.name),
                    ...))
}

#' Query to get a summary of assays, assay description, and measurement sets in a dataset
#'
#' @inheritParams get_all_gene_symbols
#'
#' @return A query suitable for returning a summary of dataset assays
#' @export
#'
dataset_summary <- function(dataset.name, ...) {
    q <- query(find(pull(?a, c(assay/name, assay/description,
                               {assay/technology = c(db/ident)},
                               {assay/measurement-sets = c(measurement-set/name)}))),
               where(d(?d, dataset/name, ?dataset.name),
                     d(?d, dataset/assays, ?a)),
               args(?dataset.name <- !!dataset.name)
    )
    return(q)
}


database_summary <- function(...){
    q <- query(find(pull(?d, c(dataset/name,
                               {dataset/assays = c(assay/name,
                                                   {assay/measurement-sets = c(measurement-set/name)})}))),
               where(d(?d, dataset/name, ?dataset.name)))
    return(q)
}
