#' @title cleanPathName_ helper function
#'
#' @description Internal function to clean paths to have a tailing slash
#'
#' @param pathname A character string. A pathname to be cleaned (to have a tailing slash).
#'
#' @export
#'
cleanPathName_ <- function(pathname){
  return(gsub("([[:alnum:]])$", "\\1/", pathname))
}

discovery_ <- function(sampleMetadata, studyMetadata){
  return(studyMetadata[which(studyMetadata$discovery == TRUE), "study"])
}
