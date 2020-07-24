#' @title gene_map_sig
#'
#' @description Helper function to map relevant input variables of a diagnostic
#'   model to corresponding IDs.
#'
#' @inheritParams plot_grid_model_plot
#' @inheritParams sigidentPrognostic
#'
#' @export
gene_map_sig <- function(mergeset, model) {
  id <- rownames(mergeset)
  # TODO warum i+1?
  index <- model[["beta"]]@i + 1
  # TODO map entrez_id on gene symbol here and include as second
  # column to ouput
  return(as.data.frame(x = cbind("ID" = id[index])))
}
