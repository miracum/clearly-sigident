#' @title createDiagnosisDesign_
#'
#' @description Helper function to create diagnosis and design
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
# create design
createDiagnosisDesign_ <- function(sampleMetadata, studyMetadata, controlname, targetname, targetcol){
  stopifnot(
    is.data.frame(sampleMetadata),
    is.data.frame(studyMetadata),
    is.character(controlname),
    is.character(targetname),
    is.character(targetcol)
  )

  discovery <- discovery_(sampleMetadata = sampleMetadata,
                          studyMetadata = studyMetadata)
  discoverydata <- sampleMetadata[which(sampleMetadata$study %in% discovery),][[targetcol]]

  diagnosis <- diagnosis_(vector = discoverydata,
                          controlname = controlname,
                          targetname = targetname)

  design <- stats::model.matrix(~diagnosis)
  return(list(diagnosis = diagnosis, design = design))
}

diagnosis_ <- function(vector, controlname, targetname){
  diag <- as.vector(vector)

  diagnosis <- gsub(controlname, "0", diag)
  diagnosis <- gsub(targetname, "1", diagnosis)

  return(as.integer(diagnosis))
}

#' @title createBatch_
#'
#' @description Helper function to create batch
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
# create batch
createBatch_ <- function(sampleMetadata, studyMetadata){

  discovery <- discovery_(sampleMetadata = sampleMetadata,
                          studyMetadata = studyMetadata)
  studylist <- list()

  for (d in discovery){
    studylist[[d]] <- sum(sampleMetadata$study == d)
  }

  x <- 1:length(studylist)

  times <- sapply(names(studylist), function(n){
    studylist[[n]]
  }, USE.NAMES = F)
  return(rep(x = x, times = times))
}


#' @title batchDetection_
#'
#' @description Helper function to detect batches
#'
#' @param batch Takes the results from \code{createBatch_()} as input.
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
batchDetection_ <- function(mergeset, batch){
  return(gPCA::gPCA.batchdetect(x = t(mergeset), batch = batch, center = FALSE))
}

#' @title batchCorrection_
#'
#' @description Helper function correcting for batch effects and mapping affy probes to Entrez IDs
#'
#' @param mergedset An ExpressionSet. The output of the function `merge_()`.
#'
#' @details This function takes a Bioconductor's ExpressionSet class (the output of the function `merge()`) and outputs a batch corrected 
#' matrix containing expression data. In order to correct for occurring batch effects and other unwanted variation in high-throughput 
#' experiments the `ComBat` function from the sva package is conducted.  
#' The affy probes are mapped to their Entrez IDs. Thereby, empty and replicated character strings are
#' removed.  
#'  
#' @inheritParams batchDetection_
#' @inheritParams goDiffReg_
#'
#' @references W.E. Johnson, C. Li, and A. Rabinovic. Adjusting batch effects in microarray data using empirical bayes methods. Biostatistics, 8(1):118–127, 2007.
#' Jeffrey T. Leek, W. Evan Johnson, Hilary S. Parker, Elana J. Fertig, Andrew E. Jaffe, John D. Storey, Yuqing Zhang and Leonardo Collado Torres (2019). sva: Surrogate Variable Analysis. R package version 3.30.1.
#'
#' @export
batchCorrection_ <- function(mergedset, batch, design){
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