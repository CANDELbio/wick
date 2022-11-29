#' Get all gene symbols from the database
#' @param ... Additional arguments to be passed to \code{\link{do_query}} (typically
#'   \code{db.name}, unless it has been set globally)
#' @return Returns a character vector with all the HGNC gene symbols in the database
#' @export
get_all_gene_symbols <- function(...) {
    return(as.vector(do_query(all_gene_symbols(), ...)))
}

#' Get all gene symbols from the database
#' @return Returns a query suitable for getting all gene symbols from the database
#' @export
all_gene_symbols <- function() {
    query(
        find(?gene-hgnc-symbol),
        where(
            d(., gene/hgnc-symbol, ?gene-hgnc-symbol)
        )
    )
}



#' Get all genes information from the database
#'
#' @inheritParams get_all_gene_symbols
#' @return Returns a \code{data.frame} with complete information about all the gene
#'   entities in the database (including alternate names and symbols)
#' @export
get_all_genes <- function(...) {
    return(do_query(all_genes(), ...))
}

#' Get all genes information from the database
#' @return Returns a query suitable for getting all the genes information
#' @export
all_genes <- function() {
    query(
        find(pull(?g, c(.))),
        where(
            d(?g, gene/hgnc-symbol)
        )
    )
}


#' Get all genomic-coordinate information from the database
#'
#' @inheritParams get_all_gene_symbols
#' @return Returns a \code{data.frame} with complete information about all the
#'   genomic coordinate entities in the database
#' @export
get_all_genomic_coordinates <- function(...) {
    return(do_query(all_genomic_coordinates(), ...))
}

#' Get all genomic-coordinate information from the database
#' @return Returns a query suitable for getting all the genomic coordinate
#'      information
#' @export
all_genomic_coordinates <- function() {
    query(
        find(pull(?g, c(.))),
        where(
            d(?g, genomic-coordinate/id)
        )
    )
}



#' Get all variant information from the database
#'
#' @inheritParams get_all_gene_symbols
#' @return Returns a \code{data.frame} with complete information about all the
#'   variant entities in the database
#' @export
get_all_variants <- function(...) {
    return(do_query(all_variants(), ...))
}

#' Get all variant information from the database
#' @return Returns a query suitable for getting all the variant
#'      information
#' @export
all_variants <- function() {
    query(
        find(pull(?v, c(.,
                        {variant/classification = c(db/ident)},
                        {variant/so-consequences = c(so-sequence-feature/name)},
                        {variant/gene = c(gene/hgnc-symbol)}))),
        where(
            d(?v, variant/id)
        )
    )
}


#' Get all gene product information from the database
#'
#' @inheritParams get_all_gene_symbols
#' @return Returns a \code{data.frame} with complete information about all the gene
#'   entities in the database (including alternate names and symbols)
#' @export
get_all_gene_products <- function(...) {
    return(do_query(all_gene_products(), ...))
}

#' Get all genes information from the database
#' @return Returns a query suitable for getting all the genes information
#' @export
all_gene_products <- function() {
    query(
        find(?gene-product-id, ?gene-hgnc-symbol),
        where(
            d(?gp, gene-product/id, ?gene-product-id),
            d(?gp, gene-product/gene, ?g),
            d(?g, gene/hgnc-symbol, ?gene-hgnc-symbol)
        )
    )
}




#' Get all GDC anatomic sites
#'
#' @return Returns a query suitable for getting all GDC anatomic sites
#' @export
all_gdc_anatomic_sites <- function() {
    query(
        find(?gdc-anatomic-site-name),
        where(
            d(?g, gdc-anatomic-site/name, ?gdc-anatomic-site-name)
        )
    )
}

#' Get all GDC anatomic sites from the database
#'
#' @inheritParams get_all_gene_symbols
#' @return Returns a character vector containing all the GDC anatomic sites in the database
#' @export
get_all_gdc_anatomic_sites <- function(...) {
    return(do_query(all_gdc_anatomic_sites(), ...))
}


#' Get all cnv from the database
#' @param ... Additional arguments to be passed to \code{\link{do_query}} (typically
#'   \code{db.name}, unless it has been set globally)
#' @return Returns a data.frame with all cnv information
#' @export
get_all_cnv <- function(...) {
    return(as.vector(do_query(all_cnv(), ...)))
}

#' @export
all_cnv <- function() {
    query(
        find(pull(?c, c(.,
                        {cnv/genomic-coordinates = c(
                                                     {genomic-coordinate/assembly = c(db/ident)})},
                        {cnv/genes = c(gene/hgnc-symbol)}))),
        where(
            d(?c, cnv/id, ?cnv-id)
        )
    )
}



#' Get all genes from the database
#' @param ... Additional arguments to be passed to \code{\link{do_query}} (typically
#'   \code{db.name}, unless it has been set globally)
#' @return Returns a character vector with all the HUGO gene names in the database
#' @export
get_all_genes_with_gc <- function(...) {
    return(as.vector(do_query(all_genes_with_gc(), ...)))
}

#' @export
all_genes_with_gc <- function() {
    query(
        find(?gene-hgnc-symbol,
             ?genomic-coordinate-id,
             ?genomic-coordinate-assembly,
             ?genomic-coordinate-contig,
             ?genomic-coordinate-strand,
             ?genomic-coordinate-start,
             ?genomic-coordinate-end),
        where(
            d(?g, gene/hgnc-symbol, ?gene-hgnc-symbol),
            d(?g, gene/genomic-coordinates, ?gc),
            d(?gc, genomic-coordinate/id, ?genomic-coordinate-id),
            d(?gc, genomic-coordinate/assembly, ?a),
            d(?a, db/ident, ?genomic-coordinate-assembly),
            d(?gc, genomic-coordinate/contig, ?genomic-coordinate-contig),
            d(?gc, genomic-coordinate/strand, ?genomic-coordinate-strand),
            d(?gc, genomic-coordinate/start, ?genomic-coordinate-start),
            d(?gc, genomic-coordinate/end, ?genomic-coordinate-end)
        )
    )
}


#' Get all proteins from the database
#'
#' @inheritParams get_all_gene_symbols
#' @return Returns a character vector with all the protein names in the database
#' @export
get_all_proteins <- function(...) {
    return(as.vector(do_query(all_proteins(), ...)))
}

#' @export
all_proteins <- function() {
    query(
        find(?protein-uniprot-name),
        where(
            d(?p, protein/uniprot-name, ?protein-uniprot-name)
        )
    )
}

#' Get all epitopes from the database
#'
#' @inheritParams get_all_gene_symbols
#' @return Returns a character vector with all the epitope ids in the database
#' @export
get_all_epitopes <- function(...) {
    return(as.vector(do_query(all_epitopes(), ...)))
}

#' @export
all_epitopes <- function() {
    query(
        find(?epitope-id),
        where(
            d(?e, epitope/id, ?epitope-id)
        )
    )
}




#' Get all cell types from the database
#'
#' @inheritParams get_all_gene_symbols
#' @return Returns a character vector with all the cell types in the
#'   database
#' @export
get_all_cell_types <- function(...) {
    return(as.vector(do_query(all_cell_types(), ...)))
}

#' @export
all_cell_types <- function() {
    query(
        find(?cell-type-co-name),
        where(
            d(?c, cell-type/co-name, ?cell-type-co-name)
        )
    )
}




#' Get all diseases from the database
#' @inheritParams get_all_gene_symbols
#' @return Returns a character vector with all the disease names in the database
#' @export
get_all_meddra_diseases <- function(...) {
    return(as.vector(do_query(all_meddra_diseases(), ...)))
}

all_meddra_diseases <- function() {
    query(
        find(?meddra-disease-preferred-name),
        where(
            d(?m, meddra-disease/preferred-name, ?meddra-disease-preferred-name)
        )
    )
}




#' Get all drugs from the database
#' @inheritParams get_all_gene_symbols
#' @return Returns a character vector with all the drug names in the database
#' @export
get_all_drugs <- function(...) {
    return(as.vector(do_query(all_drugs(), ...)))
}

all_drugs <- function() {
    query(
        find(?drug-preferred-name),
        where(
            d(?d, drug/preferred-name, ?drug-preferred-name)
        )
    )
}


#' Get all datasets from the database
#'
#' @inheritParams get_all_genes
#' @return Returns a character vector with all the dataset names in the database
#' @export
get_all_datasets <- function(...) {
    return(as.vector(do_query(all_datasets(), ...)))
}

#' @export
all_datasets <- function() {
    query(
        find(?dataset-name),
        where(
            d(?d, dataset/name, ?dataset-name)
        )
    )
}

#' Query to return all db-idents
#'
#' @return Returns a query suitable for getting all the db-idents from the database
all_db_idents <- function() {
    q <- query(find(?db-id, ?db-ident),
               where(
                   d(?db-id, db/ident, ?db-ident)
               )
    )
    return(q)
}

#' Get all db-idents
#'
#'
#' @inheritParams get_all_gene_symbols
#' @return Returns a \code{data.frame} with all the db-idents in the database
#'   and their associated entity id
#'
#' @export
get_all_db_idents <- function(...) {
    return(do_query(all_db_idents(), ...))
}

