#' @title compareDiagnosticModels
#'
#' @description Print an overview of all diagnostic models included in `modellist`.
#'
#' @param modellist A list object containing a list of diagnostic models. The output of the function `sigidentDiagnostic()`.
#'
#' @export
compareDiagnosticModels <- function(modellist){

  outdat <- data.table::data.table("Model" = character(0), "AUC" = numeric(0))

  for (n in names(modellist)){
    outdat <- rbind(outdat, data.table::data.table("Model" = n, "AUC" = round(modellist[[n]]$auc, 3)))
  }

  return(outdat)
}
