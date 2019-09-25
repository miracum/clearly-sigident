#' @title Perform DEG Analysis in Gene Expression Datasets Derived from MicroArrays
#'
#' @description One function to perform DEG analysis.
#'
#' @param mergeset A matrix of merged expression sets (rows = genes, columns = samples).
#' @param mergedset A large Expression Set. The output of the function `merge_()`.
#' @param studyMetadata A data frame. The data frame holding the study metadata.
#' @param sampleMetadata A data frame. The data frame holding the sample metadata.
#' @param targetcol A character string. Columname of `sampleMetadata` holding the targets. Default: "target". Caution: this should not be changed.
#' @param controlname A character string. Name of the the controls, specified in the 'target' column of `sampleMetadata`.
#' @param design A object. The output of the function `createDiagnosisDesignBatch_()`.
#' @param idtype A character string. The type of ID used to name the genes. One of 'entrez' or 'affy' intended to use either entrez IDs or
#'   affy IDs. Caution: when using entrez IDs, missing and duplicated IDs are being removed!
#' @param FDR A positive numeric value between (max. 0.05) indicating the desired q-Value during DEG analysis (Default: 0.01).
#' @param csvdir A character string. Path to the folder to store output tables. Default: "./tables/".
#' @param plotdir A character string. Path to the folder to store resulting plots. Default: "./plots/".
#'
#' @import data.table
#' @importFrom magrittr "%>%"
#'
#'
#' @export

sigidentDEG <- function(mergeset,
                        mergedset,
                        studyMetadata,
                        sampleMetadata,
                        targetcol,
                        controlname,
                        design,
                        idtype,
                        FDR,
                        plotdir = "./plots/",
                        csvdir = "./tables/"){
  stopifnot(
    class(mergeset) == "matrix",
    is.data.frame(sampleMetadata),
    is.data.frame(studyMetadata),
    is.character(targetcol),
    is.character(controlname),
    is.character(idtype),
    idtype %in% c("entrez", "affy"),
    is.numeric(FDR),
    FDR > 0 | FDR <= 0.05,
    is.character(plotdir),
    is.character(csvdir)
  )

  # create internal list for storage
  rv <- list()

  # store names
  rv$controlname <- controlname
  rv$targetname <- targetname
  rv$targetcol <- targetcol

  # store other variables
  rv$deg_q <- FDR
  rv$design <- design
  rv$idtype <- idtype

  # store dirs
  rv$plotdir <- cleanPathName_(plotdir)
  rv$csvdir <- cleanPathName_(csvdir)

  # add mergeset to list
  rv$mergeset <- mergeset
  rv$mergedset <- mergedset

  ### DEG Analysis ###
  rv$genes <- identifyDEGs_(mergeset = rv$mergeset,
                            design = rv$design,
                            qValue = rv$deg_q)

  # heatmap creation
  # create colors for map
  ht_colors <- colorHeatmap_(sampleMetadata = sampleMetadata,
                             studyMetadata = studyMetadata,
                             targetcol = rv$targetcol,
                             controlname = rv$controlname) # cancer = red
  createDEGheatmap_(mergeset = rv$mergeset,
                    genes = rv$genes,
                    patientcolors = ht_colors,
                    filename = paste0(rv$plotdir, "DEG_heatmap.png"))

  # Export a table with DEGs and annotations
  rv$deg_info <- sigident::exportDEGannotations_(mergedset = rv$mergedset,
                                              genes = rv$genes,
                                              idtype = rv$idtype)
  data.table::fwrite(rv$deg_info, paste0(rv$csvdir, "DEG_info.csv"))

  # export table with differential expression parameters and annotations
  rv$deg_results <- limmaTopTable_(mergeset = rv$mergeset,
                                design = rv$design,
                                qValue = rv$deg_q)
  data.table::fwrite(rv$deg_results, paste0(rv$csvdir, "DEG_results.csv"))

  # return genes
  return(genes)
}
