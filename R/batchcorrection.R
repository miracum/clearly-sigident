#' @title create_diagnosisdesignbatch
#'
#' @description Helper function to create diagnosis, design and batch
#'
#' @param targetname A character string. Name of the the targets, specified
#'   in the 'target' column of `sample_metadata`.
#'
#' @inheritParams sigidentDEG
#'
#' @export
# create design
create_diagnosisdesignbatch <- function(sample_metadata,
                                        study_metadata,
                                        controlname,
                                        targetname,
                                        targetcol) {
  stopifnot(
    is.data.frame(sample_metadata),
    is.data.frame(study_metadata),
    is.character(controlname),
    is.character(targetname),
    is.character(targetcol)
  )

  discovery <- discovery_func(sample_metadata = sample_metadata,
                          study_metadata = study_metadata)
  discoverydata <-
    sample_metadata[which(sample_metadata$study %in% discovery), ][[targetcol]]

  # create diagnosis
  diagnosis <- create_diagnosis(vector = discoverydata,
                          controlname = controlname,
                          targetname = targetname)
  # create design
  design <- stats::model.matrix(~ diagnosis)

  # create batch
  batch <- create_batch(study_metadata = study_metadata,
                        sample_metadata = sample_metadata)

  return(list(
    diagnosis = diagnosis,
    design = design,
    batch = batch
  ))
}

create_diagnosis <- function(vector, controlname, targetname) {
  diag <- as.vector(vector)

  diagnosis <- gsub(controlname, "0", diag)
  diagnosis <- gsub(targetname, "1", diagnosis)

  outdat <- as.integer(diagnosis)
  return(outdat)
}

#' @title create_batch
#'
#' @description Helper function to create batch
#'
#' @inheritParams sigidentDEG
#'
#' @export
# create batch
create_batch <- function(sample_metadata,
                         study_metadata) {

  discovery <- discovery_func(sample_metadata = sample_metadata,
                          study_metadata = study_metadata)

  studylist <- list()

  for (d in discovery) {
    studylist[[d]] <- sum(sample_metadata$study == d)
  }

  x <- seq_len(length(studylist))

  times <- sapply(names(studylist), function(n) {
    studylist[[n]]
  }, USE.NAMES = F)
  outdat <- rep(x = x, times = times)
  return(outdat)
}


#' @title batch_detection
#'
#' @description Helper function to detect batches
#'
#' @param batch Takes the results from \code{create_batch()} as input.
#'
#' @inheritParams sigidentDEG
#'
#' @export
batch_detection <- function(mergeset,
                            batch) {
  outdat <-
    gPCA::gPCA.batchdetect(x = t(mergeset),
                           batch = batch,
                           center = FALSE)
  return(outdat)
}

#' @title batch_correction
#'
#' @description Helper function correcting for batch effects and mapping
#'   affy probes to Entrez IDs
#'
#' @details This function takes a Bioconductor's ExpressionSet class (the
#'   output of the function `merge()`) and outputs a batch corrected
#'   matrix containing expression data. In order to correct for occurring
#'   batch effects and other unwanted variation in high-throughput
#'   experiments the `ComBat` function from the sva package is conducted.
#'   The affy probes are mapped to their Entrez IDs. Thereby, empty and
#'   replicated character strings are removed.
#'
#' @inheritParams batch_detection
#' @inheritParams go_diff_reg
#' @inheritParams sigidentDEG
#'
#' @references W.E. Johnson, C. Li, and A. Rabinovic. Adjusting batch effects
#' in microarray data using empirical bayes methods. Biostatistics,
#' 8(1):118–127, 2007. Jeffrey T. Leek, W. Evan Johnson, Hilary S. Parker,
#' Elana J. Fertig, Andrew E. Jaffe, John D. Storey, Yuqing Zhang and
#' Leonardo Collado Torres (2019). sva: Surrogate Variable Analysis.
#' R package version 3.30.1.
#'
#' @export
batch_correction <- function(mergedset,
                             batch,
                             design,
                             idtype) {
  # generate data frame with expression values and model matrix
  # regardarding diagnosis
  #

  df <- mergedset@assayData$exprs
  edata <- sva::ComBat(
    df,
    batch = batch,
    mod = design,
    par.prior = T
  )
  edata <- id_type(expr = edata,
                   eset = mergedset,
                   idtype = idtype)
  return(edata)
}
