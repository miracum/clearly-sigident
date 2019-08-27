# create design
createDiagnosisDesign_ <- function(mergeset, controlname, targetname, targetcol){
  stopifnot(
    is.character(controlname),
    is.character(targetname),
    is.character(targetcol)
  )

  diag <- as.vector(mergeset[[targetcol]])

  diagnosis <- gsub(controlname, "0", diag)
  diagnosis <- gsub(targetname, "1", diagnosis)

  diagnosis <- as.integer(diagnosis)
  design <- stats::model.matrix(~diagnosis)
  return(list(diagnosis = diagnosis, design = design))
}

# generate data frame with expression values and model matrix regardarding diagnosis
createDataMatrix_ <- function(mergeset){
  return(mergeset@assayData$exprs)
}

# create batch
createBatch_ <- function(eset1, eset2, eset3){
  # TODO make this working generically
  return(rep(c(1,2,3), times = c(ncol(eset1),ncol(eset2),ncol(eset3))))
}

batchCorrection_ <- function(x, batch){
  return(gPCA::gPCA.batchdetect(x = t(x), batch = batch, center = FALSE))
}

createCombat_ <- function(DF, batch, design){
  return(sva::ComBat(DF, batch = batch, mod = design, par.prior = T))
}
