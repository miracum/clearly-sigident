# create design
createDiagnosisDesign_ <- function(samplemetadata, controlname, targetname, targetcol){
  stopifnot(
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
createBatch_ <- function(study_length){
  # TODO make this working generically
  
  return(rep(c(1,2,3), times = c(ncol(eset1),ncol(eset2),ncol(eset3))))
}

batchCorrection_ <- function(x, batch){
  return(gPCA::gPCA.batchdetect(x = t(x), batch = batch, center = FALSE))
}

