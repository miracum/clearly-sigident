discovery_func <- function(sample_metadata, study_metadata) {
  return(study_metadata[which(study_metadata$discovery == TRUE), "study"])
}
