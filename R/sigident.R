#' @title Perform Signature Analyses in Gene Expression Datasets Derived from MicroArrays
#'
#' @description
#'
#' @param mergeset An ExpressionSet resulting from the \code{mergeGEO::mergeGEO()} function.
#' @param controlname A character string. Name of the the controls, specified in the 'target' column of `sampleMetadata`.
#' @param targetname A character string. Name of the the targets, specified in the 'target' column of `sampleMetadata`.
#' @param studyMetadata A data frame. The data frame holding the study metadata (output of the function \code{mergeGEO::readStudyMetadata_()}).
#' @param readSampleMetadata_ A data frame. The data frame holding the sample metadata (output of the function \code{mergeGEO::readSampleMetadata_()}).
#' @param species A character string indicating the sample's species. Currently supported: "Hs".
#' @param OrgDB A character string indicating the OrgDb. Currently supported: "org.Hs.eg.db".
#' @param organism A character string indicating the organism. Currently supported: "hsa".
#' @param pathwayid A character string indicating the pathway to show in the enrichment analysis. Currently supported: "hsa04110".
#' @param deg.q.selection A numeric value between 0 and 1 indicating the desired q-Value during DEG analysis. Default: NULL, which means, that q is
#'   calculated the following: \emph{1/length(mergeset@featureData@data$ID)}.
#' @param seed A integer value. Seed to make machine learning algorithms reproducible. Default: 111.
#' @param nfolds A integer. The number of folds used for cross validation. Default: 10.
#' @param traintest.split A numeric value between 0 and 1. The proportion of the data to be integrated into the training set for machine learning. Default: 0.8.
#' @param csvdir A character string. Path to the folder to store output tables. Default: "./tables/".
#' @param plotdir A character string. Path to the folder to store resulting plots. Default: "./plots/".
#' @param targetcol A character string. Columname of `sampleMetadata` holding the targets. Default: "target". Caution: this should not be changed.
#'
#' @import data.table
#' @importFrom magrittr "%>%"
#'
#' @seealso mergeGEO
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export

sigidentMicroarray <- function(mergeset,
                               controlname,
                               targetname,
                               studyMetadata,
                               sampleMetadata,
                               species,
                               OrgDB,
                               organism,
                               pathwayid,
                               deg.q.selection = NULL,
                               seed = 111,
                               nfolds = 10,
                               traintest.split = 0.8,
                               plotdir = "./plots/",
                               csvdir = "./tables/",
                               targetcol = "target"){
  #TODO only for debugging
  studymetadata = "lungcancer_study_metadata.csv"
  samplemetadata = "lungcancer_sample_metadata.csv"
  studyname = "lungcancer"
  metadatadir <- "./metadata/"
  dir.create(metadatadir)
  file.copy(from=system.file("./demofiles/lungcancer_study_metadata.csv", package = "mergeGEO"), to="./metadata/lungcancer_study_metadata.csv")
  file.copy(from=system.file("./demofiles/lungcancer_sample_metadata.csv", package = "mergeGEO"), to="./metadata/lungcancer_sample_metadata.csv")
  studyMetadata <- mergeGEO::readStudyMetadata_(studymetadataFilename = paste0(metadatadir, studymetadata))
  sampleMetadata <- mergeGEO::readSampleMetadata_(samplemetadataFilename = paste0(metadatadir, samplemetadata),
                                                  studyMetadata = studyMetadata)

  plotdir <- "./vignettes/plots"
  csvdir <- "./vignettes/csv"
  deg.q.selection <- NULL
  denovo = T
  controlname <- "Control"
  targetname <- "Lung Cancer"
  targetcol <- "target"
  species <- "Hs"
  OrgDb <- "org.Hs.eg.db"
  organism <- "hsa"
  pathwayid <- "hsa04110"
  seed <- 111
  traintest.split <- 0.8

  mergeset <- mergeGEO::mergeGEO(studymetadata = studymetadata,
                                  samplemetadata = samplemetadata,
                                  studyname = studyname,
                                  denovo = denovo,
                                  metadatadir = metadatadir)

  stopifnot(
    class(mergeset) == "matrix",
    is.character(plotdir),
    is.character(csvdir),
    is.character(controlname),
    is.character(targetname),
    is.character(targetcol),
    is.character(species),
    is.numeric(deg.q.selection) | is.null(deg.q.selection),
    is.numeric(seed),
    is.numeric(traintest.split),
    is.numeric(nfolds),
    traintest.split < 1 & traintest.split > 0
  )

  if (!is.null(deg.q.selection)){
    stopifnot(
      deg.q.selection > 0 | deg.q.selection < 1
    )
  }

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
  rv$nfolds <- nfolds
  rv$traintest.split <- traintest.split

  # add mergedset to list
  rv$mergeset <- mergeset


  ### Fileimport ###
  # visualize log2 transformed expression values of the merged data set
  createImportBoxplot_(mergeset = rv$mergeset, filename = paste0(rv$plotdir, "import_boxplot.png"))


  ### Batchcorrection ###
  # get diagnosis and design
  dd <- createDiagnosisDesign_(sampleMetadata = sampleMetadata,
                               studyMetadata = studyMetadata,
                               controlname = rv$controlname,
                               targetname = rv$targetname,
                               targetcol = rv$targetcol)
  rv$diagnosis <- dd$diagnosis
  rv$design <- dd$design

  rv$batch <- createBatch_(studyMetadata = studyMetadata,
                           sampleMetadata = sampleMetadata)

  rv$gPCA_after <- batchCorrection_(mergeset = rv$mergeset, batch = rv$batch)
  filename <- paste0(rv$plotdir, "PCplot_after.png")
  createBatchPlot_(correction_obj = rv$gPCA_after, filename = filename, time = "after")


  ### DEG Analysis ###
  rv$deg_q <- qSelection_(sampleMetadata = sampleMetadata,
                          deg.q.selection = deg.q.selection)

  rv$genes <- identify.DEGs_(mergeset = rv$mergeset, design = rv$design, qValue = rv$deg_q)

  # heatmap creation
  filename <- paste0(rv$plotdir, "DEG_heatmap.png")
  # create colors for map
  ht_colors <- colorHeatmap_(sampleMetadata = sampleMetadata, targetcol = rv$targetcol, controlname = rv$controlname) # cancer = red
  createDEGheatmap_(mergeset = rv$mergeset, genes = rv$genes, patientcolors = ht_colors, filename = filename)


  # TODO start here 13.9. Lorenz
  # rv$deg_info <- exportDEGannotations_(mergeset = rv$mergeset,
  #                                      genes = rv$genes)
  # data.table::fwrite(rv$deg_info, paste0(rv$csvdir, "DEG_info.csv"))
  #
  # rv$deg_results <- limmaTopTable_(mergeset = rv$combat,
  #                                  design = rv$design,
  #                                  qValue = rv$deg_q)
  # data.table::fwrite(rv$deg_results, paste0(rv$csvdir, "DEG_results.csv"))

  # gene enrichment
  rv$deg_entrez <- unique(rv$genes)
  # test for over-representation of gene ontology terms
  rv$enr_topgo <- extractGOterms_(entrez = rv$deg_entrez,
                                  species = rv$species)
  data.table::fwrite(rv$enr_topgo, paste0(rv$csvdir, "Top_GO.csv"))

  # test for over-representation of KEGG pathways
  rv$enr_topkegg <- extractKEGGterms_(entrez = rv$deg_entrez,
                                      species = rv$species)
  data.table::fwrite(rv$enr_topkegg, paste0(rv$csvdir, "Top_KEGG.csv"))

  # take differential regulation between two groups (design) into account
  rv$enr_fitlm <- goDiffReg_(mergeset = rv$mergeset,
                             design = rv$design)
  # test for over-representation of gene ontology terms
  rv$enr_fitlm_topgo <- extractGOterms_(entrez = rv$enr_fitlm,
                                        species = rv$species,
                                        FDR = 0.01)
  data.table::fwrite(rv$enr_fitlm_topgo, paste0(rv$csvdir, "Top_GO_fitlm.csv"))

  # test for over-representation of KEGG pathways
  rv$enr_fitlm_topkegg <- extractKEGGterms_(entrez = rv$enr_fitlm,
                                            species = rv$species)
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
  createEnrichtedBarplot_(enrichmentobj = rv$enr_analysis$go,
                          type = "GO",
                          filename = paste0(rv$plotdir, "Enriched_GO.png"))
  createEnrichtedBarplot_(enrichmentobj = rv$enr_analysis$kegg,
                          type = "KEGG",
                          filename = paste0(rv$plotdir, "Enriched_KEGG.png"))


  # identification of Diagnostic Signature
  # first, create training_list
  rv$training_list <- createTrainingTest_(diagnosis = rv$diagnosis,
                                          mergeset = rv$mergeset,
                                          split = rv$traintest.split,
                                          seed = rv$seed)


  # Lasso regression
  rv$diagnostic_lasso <- glmnetSignature_(traininglist = rv$training_list,
                                          type = "lasso",
                                          nfolds = rv$nfolds,
                                          seed = rv$seed)
  createCVPlot_(cv_obj = rv$diagnostic_lasso$fitCV,
                filename = paste0(rv$plotdir, "CV_lasso.png"))
  createROCplot_(roc = rv$diagnostic_lasso$roc.min,
                 filename = paste0(rv$plotdir, "ROC_Lasso.min.png"))
  createROCplot_(roc = rv$diagnostic_lasso$roc.1se,
                 filename = paste0(rv$plotdir, "ROC_Lasso.1se.png"))


  # Elastic net regression
  rv$diagnostic_elasticnet <- glmnetSignature_(traininglist = rv$training_list,
                                               type = "elastic",
                                               alpha = 0.9,
                                               nfolds = rv$nfolds,
                                               seed = rv$seed)
  createCVPlot_(cv_obj = rv$diagnostic_elasticnet$fitCV,
                filename = paste0(rv$plotdir, "CV_elasticNet.png"))
  createROCplot_(roc = rv$diagnostic_elasticnet$roc.min,
                 filename = paste0(rv$plotdir, "ROC_elasticNet.min.png"))
  createROCplot_(roc = rv$diagnostic_elasticnet$roc.1se,
                 filename = paste0(rv$plotdir, "ROC_elasticNet.1se.png"))

  # with both calculated hyperparameters alpha and lambda applying grid search
  rv$diagnostic_glmGrid <- glmnetSignature_(traininglist = rv$training_list,
                                            type = "grid",
                                            nfolds = rv$nfolds,
                                            seed = rv$seed)
  # plot model of gridsearch
  createGridModelPlot_(model = rv$diagnostic_glmGrid$caret.train,
                       filename = paste0(rv$plotdir, "Gridsearch_model.png"))
  # plot variable importance of gridsearch
  createGridVarImpPlot_(model = rv$diagnostic_glmGrid$caret.train,
                        filename = paste0(rv$plotdir, "Gridsearch_variable_importance.png"))
  # create roc plot
  createROCplot_(roc = rv$diagnostic_glmGrid$roc.elasticNet,
                 filename = paste0(rv$plotdir, "ROC_elasticNet.grid.png"))


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
