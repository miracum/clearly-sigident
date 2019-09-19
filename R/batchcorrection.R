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

  diag <- as.vector(discoverydata)

  diagnosis <- gsub(controlname, "0", diag)
  diagnosis <- gsub(targetname, "1", diagnosis)

  diagnosis <- as.integer(diagnosis)
  design <- stats::model.matrix(~diagnosis)
  return(list(diagnosis = diagnosis, design = design))
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
