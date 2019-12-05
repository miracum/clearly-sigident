#' @title export_deg_annotations
#'
#' @description Helper function to export DEG annotations from mergedset.
#'
#' @param mergedset A large Expression Set. The output of the function
#'   `sigident.preproc::geo_merge`. Please note, that mergedset holds data,
#'   which are not yet batch corrected.
#' @param idtype A character string. The type of ID used to name the
#'   genes. One of 'entrez' or 'affy' intended to use either entrez IDs or
#'   affy IDs. Caution: when using entrez IDs, missing and duplicated IDs
#'   are being removed!
#'
#' @inheritParams plot_deg_heatmap
#'
#' @seealso sigident.preproc
#'
#' @export
export_deg_annotations <- function(mergedset,
                                   genes,
                                   idtype) {

  stopifnot(idtype %in% c("entrez", "affy"))
  ids <- mergedset@featureData@data[, "ID"]
  sym <- Biobase::fData(mergedset)["Gene Symbol"][ids, ]
  tit <- Biobase::fData(mergedset)["Gene Title"][ids, ]
  gb_acc <- Biobase::fData(mergedset)["GB_ACC"][ids, ]
  entrez <- Biobase::fData(mergedset)["ENTREZ_GENE_ID"][ids, ]
  degs_info <- data.table::data.table(cbind(ids,
                                            sym,
                                            tit,
                                            gb_acc,
                                            entrez))
  colnames(degs_info) <-
    c("probe_ID",
      "gene_symbol",
      "gene_title",
      "genebank_accession",
      "entrez_id")

  if (idtype == "entrez") {
    ret <- degs_info[get("entrez_id") %in% genes, ]
  } else if (idtype == "affy") {
    ret <- degs_info[get("probe_ID") %in% genes, ]
  }

  return(ret)
}
