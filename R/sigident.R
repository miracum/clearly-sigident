#' @title Perform Signature Analyses in Gene Expression Datasets
#'
#' @description
#'
#' @import data.table
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export

sigident <- function(mergedset, plotdir, csvdir, controlname, deg.q.selection = NULL){
  #TODO only for debugging
  plotdir <- "./tests/testthat/plots"


  stopifnot(
    is.data.frame(mergedset) | data.table::is.data.table(mergedset),
    is.character(plotdir),
    is.numeric(deg.q.selection)
  )

  # create internal list for storage
  rv <- list()

  # store dirs
  rv$plotdir <- cleanPathName_(plotdir)
  rv$csvdir <- cleanPathName_(plcsvdirotdir)

  # create output directories
  dir.create(rv$plotdir)
  dir.create(rv$csvdir)

  #TODO only for debugging
  load("./tests/testthat/testdata/esets.RData")
  esets <- c(eset1b, eset2b, eset3b)
  mergedset <- mergeEsets_(esets)

  ### Fileimport ###
  # add mergedset to list
  rv$mergeset <- mergedset

  # visualize log2 transformed expression values of the merged data set
  createImportHistogram_(mergeset = rv$mergeset, filename = paste0(rv$plotdir, "import_histogram.png"))


  ### Batchcorrection ###
  # get diagnosis and design
  dd <- createDiagnosisDesign_(mergeset = rv$mergeset, controlname = "Control", targetname = "Lung Cancer")
  rv$diagnosis <- dd[["diagnosis"]]
  rv$design <- dd[["design"]]

  rv$DF <- createDataMatrix_(rv$mergeset)
  rv$batch <- createBatch_(eset1b, eset2b, eset3b)
  rv$combat <- createCombat_(rv$DF, rv$batch, rv$design)

  rv$gPCA_before <- batchCorrection_(rv$DF, rv$batch)
  filename <- paste0(rv$plotdir, "PCplot_before.png")
  createBatchPlot_(rv$gPCA_before, filename, "before")

  rv$gPCA_after <- batchCorrection_(rv$combat, rv$batch)
  filename <- paste0(rv$plotdir, "PCplot_after.png")
  createBatchPlot_(rv$gPCA_after, filename, "after")


  ### DEG Analysis ###
  rv$deg_q <- qSelection_(rv$mergeset, deg.q.selection)

  rv$genes <- identify.DEGs(BatchRemovedExprs = rv$combat, design = rv$design, qValue = rv$deg_q)

  # heatmap creation
  filename <- paste0(rv$plotdir, "DEG_heatmap.png")
  # create colors for map
  ht_colors <- unlist(lapply(rv$mergeset$Tumor, colorHeatmap_, controlname)) # cancer = red
  createDEGheatmap_(rv$combat, rv$genes, ht_colors, filename)


  rv$genelist <- geneMapping_(rv$mergeset, rv$genes)
  rv$deg_info <- exportDEGannotations_(rv$mergeset, rv$genes)
  data.table::fwrite(rv$deg_info, paste0(csvdir, "DEG_info.csv"))

  rv$deg_results <- limmaTopTable_(rv$combat, rv$design, rv$deg_p)
  data.table::fwrite(rv$deg_results, paste0(csvdir, "DEG_results.csv"))
}
