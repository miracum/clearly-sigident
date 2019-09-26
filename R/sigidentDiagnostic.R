#' @title Perform Diagnostic Signature Analyses in Gene Expression Datasets Derived from MicroArrays
#'
#' @description One function to perform diagnostic signature analysis.
#'
#' @param diagnosis An object. The output of the function `createDiagnosisDesignBatch_()`.
#' @param seed A integer value. Seed to make machine learning algorithms reproducible. Default: 111.
#' @param nfolds A integer. The number of folds used for cross validation. Default: 10.
#' @param split A numeric value between 0 and 1. The proportion of the data to be integrated into the training set for machine learning. Default: 0.8.
#'
#' @inheritParams sigidentDEG
#'
#'
#' @export
sigidentDiagnostic <- function(mergeset,
                               diagnosis,
                               seed = 111,
                               nfolds = 10,
                               split = 0.8,
                               plotdir = "./plots/"){

  stopifnot(
    class(mergeset) == "matrix",
    is.numeric(diagnosis),
    is.numeric(seed),
    is.numeric(split),
    is.numeric(nfolds),
    split < 1 & split > 0
  )

  # create internal list for storage
  rv <- list()

  # store dirs
  rv$plotdir <- cleanPathName_(plotdir)

  # create output directories
  dir.create(rv$plotdir)

  # set diagnosis
  rv$diagnosis <- diagnosis

  # store seed, traintest.split
  rv$seed <- seed
  rv$nfolds <- nfolds
  rv$traintest.split <- split

  # add mergeset to list
  rv$mergeset <- mergeset

  # identification of Diagnostic Signature
  # first, create training_list
  rv$training_list <- createTrainingTest_(diagnosis = rv$diagnosis,
                                          mergeset = rv$mergeset,
                                          split = rv$traintest.split,
                                          seed = rv$seed)


  # Lasso regression
  rv$diagnostic_lasso <- signature_(traininglist = rv$training_list,
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
  rv$diagnostic_elasticnet <- signature_(traininglist = rv$training_list,
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
  rv$diagnostic_glmGrid <- signature_(traininglist = rv$training_list,
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
  diagnosticModels <- list(
    "lasso" = list("CV" = rv$diagnostic_lasso$fitCV,
                   "min" = list("model" = rv$diagnostic_lasso$lambda.min,
                                "confmat" = rv$diagnostic_lasso$confmat.min,
                                "prediction" = rv$diagnostic_lasso$predicted.min,
                                "auc" = as.numeric(rv$diagnostic_lasso$roc.min$auc)),
                   "1se" = list("model" = rv$diagnostic_lasso$lambda.1se,
                                "confmat" = rv$diagnostic_lasso$confmat.1se,
                                "prediction" = rv$diagnostic_lasso$predicted.1se,
                                "auc" = as.numeric(rv$diagnostic_lasso$roc.1se$auc))
    ),
    "elasticnet" = list("CV" = rv$diagnostic_elasticnet$fitCV,
                        "min" = list("model" = rv$diagnostic_elasticnet$lambda.min,
                                     "confmat" = rv$diagnostic_elasticnet$confmat.min,
                                     "prediction" = rv$diagnostic_elasticnet$predicted.min,
                                     "auc" = as.numeric(rv$diagnostic_elasticnet$roc.min$auc)),
                        "1se" = list("model" = rv$diagnostic_elasticnet$lambda.1se,
                                     "confmat" = rv$diagnostic_elasticnet$confmat.1se,
                                     "prediction" = rv$diagnostic_elasticnet$predicted.1se,
                                     "auc" = as.numeric(rv$diagnostic_elasticnet$roc.1se$auc))
    ),
    "grid" = list("CV" = rv$diagnostic_glmGrid$caret.train,
                  "model" = rv$diagnostic_glmGrid$elasticNet.auto,
                  "confmat" = rv$diagnostic_glmGrid$confmat.elasticNet,
                  "prediction" = rv$diagnostic_elasticnet$predicted.elasticNet,
                  "auc" = as.numeric(rv$diagnostic_glmGrid$roc.elasticNet$auc))
  )
  return(diagnosticModels = diagnosticModels)
}
