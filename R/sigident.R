#' @title Perform Signature Analyses in Gene Expression Datasets Derived from MicroArrays
#'
#' @description
#'
#' @param mergedset An expression set merged by the method of Hughey and Butte.
#'
#' @import data.table
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export

sigidentMicroarray <- function(mergedset, plotdir, csvdir, targetcol, controlname, targetname, species, deg.q.selection = NULL, seed = 111, traintest.split = 0.8){
  #TODO only for debugging
  plotdir <- "./tests/testthat/plots"
  csvdir <- "./tests/testthat/csvs"
  deg.q.selection <- NULL
  controlname <- "Control"
  targetname <- "Lung Cancer"
  targetcol <- "Tumor"
  species <- "Hs"
  OrgDb <- "org.Hs.eg.db"
  organism <- "hsa"
  pathwayid <- "hsa04110"
  seed <- 111
  traintest.split <- 0.8

  #TODO only for debugging
  load("./tests/testthat/testdata/esets.RData")
  esets <- c(eset1b, eset2b, eset3b)
  mergedset <- mergeEsets_(esets)


  stopifnot(
    class(mergedset) == "ExpressionSet",
    is.character(plotdir),
    is.character(csvdir),
    is.character(controlname),
    is.character(targetname),
    is.character(targetcol),
    is.character(species),
    is.numeric(deg.q.selection) | is.null(deg.q.selection),
    is.numeric(seed),
    is.numeric(traintest.split),
    traintest.split < 1 & traintest.split > 0
  )

  # create internal list for storage
  rv <- list()

  # store names
  rv$controlname <- controlname
  rv$targetname <- targetname
  rv$targetcol <- targetcol

  # store species, orgdb and orgamism
  rv$species <- species
  rv$orgdb <- OrgDb
  rv$organism <- organism
  rv$pathwayid <- pathwayid

  # store dirs
  rv$plotdir <- cleanPathName_(plotdir)
  rv$csvdir <- cleanPathName_(csvdir)

  # create output directories
  dir.create(rv$plotdir)
  dir.create(rv$csvdir)

  # store seed, traintest.split
  rv$seed <- seed
  rv$traintest.split <- traintest.split

  # add mergedset to list
  rv$mergeset <- mergedset


  ### Fileimport ###
  # visualize log2 transformed expression values of the merged data set
  createImportHistogram_(mergeset = rv$mergeset, filename = paste0(rv$plotdir, "import_histogram.png"))


  ### Batchcorrection ###
  # get diagnosis and design
  dd <- createDiagnosisDesign_(mergeset = rv$mergeset, controlname = rv$controlname, targetname = rv$targetname, targetcol = rv$targetcol)
  rv$diagnosis <- dd$diagnosis
  rv$design <- dd$design

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

  rv$genes <- identify.DEGs_(BatchRemovedExprs = rv$combat, design = rv$design, qValue = rv$deg_q)

  # heatmap creation
  filename <- paste0(rv$plotdir, "DEG_heatmap.png")
  # create colors for map
  ht_colors <- colorHeatmap_(mergeset = rv$mergeset, targetcol = rv$targetcol, controlname = rv$controlname) # cancer = red
  createDEGheatmap_(rv$combat, rv$genes, ht_colors, filename)


  # create genelist
  rv$genelist <- geneMapping_(rv$mergeset, rv$genes)
  rv$deg_info <- exportDEGannotations_(rv$mergeset, rv$genes)
  data.table::fwrite(rv$deg_info, paste0(rv$csvdir, "DEG_info.csv"))

  rv$deg_results <- limmaTopTable_(rv$combat, rv$design, rv$deg_q)
  data.table::fwrite(rv$deg_results, paste0(rv$csvdir, "DEG_results.csv"))

  # gene enrichment
  rv$deg_entrez <- enrichmentDataSelection_(mergeset = rv$mergeset, genes = rv$genes)
  # test for over-representation of gene ontology terms
  rv$enr_topgo <- extractGOterms_(entrez = rv$deg_entrez, species = rv$species)
  data.table::fwrite(rv$enr_topgo, paste0(rv$csvdir, "Top_GO.csv"))

  # test for over-representation of KEGG pathways
  rv$enr_topkegg <- extractKEGGterms_(entrez = rv$deg_entrez, species = rv$species)
  data.table::fwrite(rv$enr_topkegg, paste0(rv$csvdir, "Top_KEGG.csv"))

  # take differential regulation between two groups (design) into account
  rv$enr_fitlm <- goDiffReg_(mergeset = rv$mergeset, BatchRemovedExprs = rv$combat, design = rv$design)
  # test for over-representation of gene ontology terms
  rv$enr_fitlm_topgo <- extractGOterms_(entrez = rv$enr_fitlm, species = rv$species, FDR = 0.01)
  data.table::fwrite(rv$enr_fitlm_topgo, paste0(rv$csvdir, "Top_GO_fitlm.csv"))

  # test for over-representation of KEGG pathways
  rv$enr_fitlm_topkegg <- extractKEGGterms_(entrez = rv$enr_fitlm, species = rv$species)
  data.table::fwrite(rv$enr_fitlm_topkegg, paste0(rv$csvdir, "Top_KEGG_fitlm.csv"))

  # perform enrichment analysis
  rv$enr_analysis <- goEnrichmentAnalysis_(entrez = rv$deg_entrez,
                                           OrgDB = rv$orgdb,
                                           organism = rv$organism,
                                           fitlm = rv$enr_fitlm,
                                           pathwayid = rv$pathwayid,
                                           species = rv$organism,
                                           plotdir = rv$plotdir)

  # plotting enrichmentanalysis
  createEnrichtedBarplot_(enrichmentobj = rv$enr_analysis$go, type = "GO", filename = paste0(rv$plotdir, "Enriched_GO.png"))
  createEnrichtedBarplot_(enrichmentobj = rv$enr_analysis$kegg, type = "KEGG", filename = paste0(rv$plotdir, "Enriched_KEGG.png"))


  # identification of Diagnostic Signature
  # first, create training_list
  rv$training_list <- createTrainingTest_(rv$diagnosis, rv$combat, split = rv$traintest.split, seed = rv$seed)


  # Lasso regression
  rv$diagnostic_lasso <- glmnetSignature_(traininglist = rv$training_list, type = "lasso", nfolds = 10, seed = rv$seed)
  createCVPlot_(cv_obj = rv$diagnostic_lasso$fitCV, filename = paste0(rv$plotdir, "CV_lasso.png"))
  createROCplot_(roc = rv$diagnostic_lasso$roc.min, filename = paste0(rv$plotdir, "ROC_Lasso.min.png"))
  createROCplot_(roc = rv$diagnostic_lasso$roc.1se, filename = paste0(rv$plotdir, "ROC_Lasso.1se.png"))


  # Elastic net regression
  rv$diagnostic_elasticnet <- glmnetSignature_(traininglist = rv$training_list, type = "elastic", alpha = 0.9, nfolds = 10, seed = rv$seed)
  createCVPlot_(cv_obj = rv$diagnostic_elasticnet$fitCV, filename = paste0(rv$plotdir, "CV_elasticNet.png"))
  createROCplot_(roc = rv$diagnostic_elasticnet$roc.min, filename = paste0(rv$plotdir, "ROC_elasticNet.min.png"))
  createROCplot_(roc = rv$diagnostic_elasticnet$roc.1se, filename = paste0(rv$plotdir, "ROC_elasticNet.1se.png"))

  # with both calculated hyperparameters alpha and lambda applying grid search
  rv$diagnostic_glmGrid <- glmnetSignature_(traininglist = rv$training_list, type = "grid", nfolds = 10, seed = rv$seed)
  # plot model of gridsearch
  createGridModelPlot_(model = rv$diagnostic_glmGrid$caret.train, filename = paste0(rv$plotdir, "Gridsearch_model.png"))
  # plot variable importance of gridsearch
  createGridVarImpPlot_(model = rv$diagnostic_glmGrid$caret.train, filename = paste0(rv$plotdir, "Gridsearch_variable_importance.png"))
  # create roc plot
  createROCplot_(roc = rv$diagnostic_glmGrid$roc.elasticNet, filename = paste0(rv$plotdir, "ROC_elasticNet.grid.png"))


  # compare aucs
  model.list <- list(
    lasso.min = list(model = rv$diagnostic_lasso$lambda.min,
                     confmat = rv$diagnostic_lasso$confmat.min,
                     auc = as.numeric(rv$diagnostic_lasso$roc.min$auc)),
    lasso.1se = list(model = rv$diagnostic_lasso$lambda.1se,
                     confmat = rv$diagnostic_lasso$confmat.1se,
                     auc = as.numeric(rv$diagnostic_lasso$roc.1se$auc)),
    elastic.min = list(model = rv$diagnostic_elasticnet$lambda.min,
                       confmat = rv$diagnostic_elasticnet$confmat.min,
                       auc = as.numeric(rv$diagnostic_elasticnet$roc.min$auc)),
    elastic.1se = list(model = rv$diagnostic_elasticnet$lambda.1se,
                       confmat = rv$diagnostic_elasticnet$confmat.1se,
                       auc = as.numeric(rv$diagnostic_elasticnet$roc.1se$auc)),
    elastic.grid = list(model = rv$diagnostic_glmGrid$elasticNet.auto,
                        confmat = rv$diagnostic_glmGrid$confmat.elasticNet,
                        auc = as.numeric(rv$diagnostic_glmGrid$roc.elasticNet$auc))
  )
  return(model.list)
}
