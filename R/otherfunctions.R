#' @title merge_
#'
#' @description Helper function to merge expression sets
#'
#' @param esetlist A list of expression sets to merge
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
    #Merging Expressionsets
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
    #Merge Expressiondata
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

#' @title createCombat_
#'
#' @description Helper function to create a matrix with genes as rows and samples as columns
#'
#' @param mergedset A large expression set. The output of the function `merge_()`.
#'
#' @inheritParams batchCorrection_
#' @inheritParams goDiffReg_
#'
#' @export
createCombat_ <- function(mergedset, batch, design){
  # generate data frame with expression values and model matrix regardarding diagnosis
  DF <- mergedset@assayData$exprs
  edata <- sva::ComBat(DF, batch = batch, mod = design, par.prior = T)
  # mapping Entrez-IDs to expression matrix
  rownames(edata) <- as.character(mergedset@featureData@data$ENTREZ_GENE_ID)
  # remove empty characters and replicates in EntrezIDs
  edata <- edata[rownames(edata)!="",]
  edata <- edata[!duplicated(rownames(edata)),]
  return(edata)
}
