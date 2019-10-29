#' @title clean_path_name helper function
#'
#' @description Internal function to clean paths to have a
#'   tailing slash
#'
#' @param pathname A character string. A pathname to be
#'   cleaned (to have a tailing slash).
#'
#' @export
#'
clean_path_name <- function(pathname){
  return(gsub("([[:alnum:]])$", "\\1/", pathname))
}

discovery_func <- function(sample_metadata, study_metadata){
  return(study_metadata[which(study_metadata$discovery == TRUE), "study"])
}
