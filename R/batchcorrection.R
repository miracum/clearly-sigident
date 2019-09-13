# create design
createDiagnosisDesign_ <- function(samplemetadata, controlname, targetname, targetcol){
  stopifnot(
    is.data.frame(samplemetadata),
    is.character(controlname),
    is.character(targetname),
    is.character(targetcol)
  )

  diag <- as.vector(samplemetadata[[targetcol]])

  diagnosis <- gsub(controlname, "0", diag)
  diagnosis <- gsub(targetname, "1", diagnosis)

  diagnosis <- as.integer(diagnosis)
  design <- stats::model.matrix(~diagnosis)
  return(list(diagnosis = diagnosis, design = design))
}

# create batch
createBatch_ <- function(sampleMetadata, studyMetadata){

  discovery <- studyMetadata[which(studyMetadata$discovery), "study"]
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

batchCorrection_ <- function(x, batch){
  return(gPCA::gPCA.batchdetect(x = t(x), batch = batch, center = FALSE))
}
