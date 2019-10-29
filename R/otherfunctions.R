#' @title merge_
#'
#' @description Helper function to merge ExpressionSets
#'
#' @param esetlist A list of ExpressionSets with consistent phenoData to merge
#'
#' @export
merge_ <- function(esetlist) {
  # set1 = first eset in list, will be overwritten by combined esets
  set1 <- esetlist[[1]]

  for (i in 2:length(esetlist)) {
    # loop combines esets one by one

    set2 <- esetlist[[i]] # (selecting the next eset)
    # Expression-Matrices
    e1 <- Biobase::exprs(set1)
    e2 <- Biobase::exprs(set2)
    #% dim(e1)
    #% dim(e2)
    # get list of all overlapping genes
    # overlap is used to selectonly overlapping genes -> genes that are only
    # present in one eset will be discarded
    overlap <- sort(intersect(rownames(e1), rownames(e2)))
    # --> allows for combination of different chips (with differing feature
    # numbers), as long as the IDs are mapped
    # Merging ExpressionSets
    #
    # Generate fData
    f_data_new <- Biobase::fData(set1)[overlap, ]

    # Generate pData
    p1 <- Biobase::pData(set1)
    p2 <- Biobase::pData(set2)

    # pData is bound rowwise, therefore sample-order is unaffected
    p_data_new <- rbind(p1, p2)
    #% dim(p1)
    #% dim(p2)
    # Annotation
    #% annoNew = c(annotation(set1),annotation(set2))
    # Merge expression data

    # expression matrices are bound columnwise, therefor
    # sample-order is unaffected
    e_new <- cbind(e1[overlap, ], e2[overlap, ])

    #% dim(e_new)
    # Rebuild eset
    # new eset has to build from scratch; cant just be appended
    # (like insilicomerging tried to do)
    set1 <- methods::new("ExpressionSet", exprs = e_new)
    Biobase::pData(set1) <- p_data_new
    Biobase::fData(set1) <- f_data_new
    #% annotation(eset_new) = unique(annoNew)
  }

  # overwriting set1 with newly combined eset, then jumping back in the loop
  eset_new <- set1
  # loop is finished after combining all esets; result is returned here
  return(eset_new)
}


#' @title export_deg_annotations
#'
#' @description Helper function to export DEG annotations from mergedset.
#'
#' @inheritParams batch_correction
#' @inheritParams plot_deg_heatmap
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
