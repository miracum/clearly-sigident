DEG.limma_ <- function(mergeset, design){
  fit <- fitLimma_(mergeset, design)
  t <- limma::topTable(fit, coef=2,number=Inf,p.value=0.05,lfc=2)
  genes <- rownames(t)
  return(genes)
}

#' @title identifyDEGs_
#'
#' @description Helper function to identify DEGs based on the limma package
#'
#' @param design A object. The output of the function `createDiagnosisDesign_()`.
#' @param qValue A numeric value. The output of the function `qSelection_()`.
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
identifyDEGs_ <- function(mergeset, design, qValue){
  fit <- fitLimma_(mergeset, design)
  t <- limma::topTable(fit, coef=2,number=Inf,p.value=qValue, lfc=2, adjust.method = "BH")
  genes <- rownames(t)
  return(genes)
}

fitLimma_ <- function(mergeset, design){
  return(limma::eBayes(limma::lmFit(mergeset, design)))
}


#' @title limmaTopTable_
#'
#' @description Helper function to get DEG results
#'
#' @inheritParams sigidentMicroarray
#' @inheritParams identifyDEGs_
#'
#' @export
limmaTopTable_ <- function(mergeset, design, qValue){
  fit <- fitLimma_(mergeset, design)
  t <- limma::topTable(fit,
                       coef = 2,
                       number = Inf,
                       p.value = qValue,
                       lfc = 2,
                       adjust.method = "BH")
  t[,2:4] <- NULL
  t[,3] <- NULL
  t <- cbind(rownames(t),t)
  colnames(t) <- c("Probe ID", "logFC", "adj.qValue")
  return(t)
}


#' @title qSelection_
#'
#' @description Helper function to select qValues
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
qSelection_ <- function(sampleMetadata, studyMetadata, deg.q.selection = NULL){
  if (is.null(deg.q.selection)){
    discovery <- discovery_(sampleMetadata = sampleMetadata,
                            studyMetadata = studyMetadata)
    deg_q <- 1/length(sampleMetadata[sampleMetadata$study %in% discovery, "sample"])
  } else {
    deg_q <- as.numeric(deg.q.selection)
  }
  return(deg_q)
}
