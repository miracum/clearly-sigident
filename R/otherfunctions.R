#' @title merge_
#'
#' @description Helper function to merge ExpressionSets
#'
#' @param esetlist A list of ExpressionSets with consistent phenoData to merge
#'
#' @export
merge_ <- function(esetlist){

  set1 <- esetlist[[1]] # set1 = first eset in list, will be overwritten by combined esets

  for(i in 2:length(esetlist)){ # loop combines esets one by one

    set2 <- esetlist[[i]] # (selecting the next eset)
    # Expression-Matrices
    e1 <- Biobase::exprs(set1)
    e2 <- Biobase::exprs(set2)
    #dim(e1)
    #dim(e2)
    #get list of all overlapping genes
    overlap <- sort(intersect(rownames(e1),rownames(e2))) # overlap is used to select only overlapping genes -> genes that are only present in one eset will be discarded
    # --> allows for combination of different chips (with differing feature numbers), as long as the IDs are mapped
    #Merging ExpressionSets
    #Generate fData
    fDataNew = Biobase::fData(set1)[overlap,]
    #Generate pData
    p1 <- Biobase::pData(set1)
    p2 <- Biobase::pData(set2)
    pDataNew <- rbind(p1,p2) # pData is bound rowwise, therefore sample-order is unaffected
    #dim(p1)
    #dim(p2)
    #Annotation
    #annoNew = c(annotation(set1),annotation(set2))
    #Merge expression data
    eNew <- cbind(e1[overlap,],e2[overlap,]) # expression matrices are bound columnwise, therefor sample-order is unaffected
    #dim(eNew)
    # Rebuild eset
    set1 <- methods::new("ExpressionSet", exprs=eNew) # new eset has to build from scratch; cant just be appended (like insilicomerging tried to do)
    Biobase::pData(set1) = pDataNew
    Biobase::fData(set1) = fDataNew
    #annotation(esetNew) = unique(annoNew)
  }
  esetNew <- set1 # overwriting set1 with newly combined eset, then jumping back in the loop
  return(esetNew) # loop is finished after combining all esets; result is returned here
}


#' @title exportDEGannotations_
#'
#' @description Helper function to export DEG annotations from mergedset.
#'
#' @inheritParams batchCorrection_
#' @inheritParams createDEGheatmap_
#'
#' @export
exportDEGannotations_ <- function(mergedset, genes, idtype){
  stopifnot(
    idtype %in% c("entrez", "affy")
  )
  ids <- mergedset@featureData@data[,"ID"]
  sym <- Biobase::fData(mergedset)["Gene Symbol"][ids,]
  tit <- Biobase::fData(mergedset)["Gene Title"][ids,]
  gbACC <- Biobase::fData(mergedset)["GB_ACC"][ids,]
  Entrez <- Biobase::fData(mergedset)["ENTREZ_GENE_ID"][ids,]
  DEGsInfo <- data.table::data.table(cbind(ids,
                                           sym,
                                           tit,
                                           gbACC,
                                           Entrez))
  colnames(DEGsInfo) <- c("probe_ID","gene_symbol","gene_title","genebank_accession","entrez_id")

  if (idtype == "entrez"){
    ret <- DEGsInfo[get("entrez_id") %in% genes, ]
  } else if (idtype == "affy"){
    ret <- DEGsInfo[get("probe_ID") %in% genes, ]
  }

  return(ret)
}
