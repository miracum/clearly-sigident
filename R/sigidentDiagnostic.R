#' @title Perform Diagnostic Signature Analyses in Gene Expression Datasets
#'   Derived from MicroArrays
#'
#' @description One function to perform diagnostic signature analysis.
#'
#' @param diagnosis An object. The output of the function
#'   `create_diagnosisdesignbatch()`.
#' @param seed A integer value. Seed to make machine learning algorithms
#'   reproducible. Default: 111.
#' @param nfolds A integer. The number of folds used for cross validation
#'   (default: 10).
#' @param split A numeric value between 0 and 1. The proportion of the data
#'   to be integrated into the training set for machine learning
#'   (default: 0.8).
#' @param repeats An integer. The number of repeated cross validations
#'   (default = 5).
#' @param tunelength An integer. The caret tuning length (default = 10).
#' @param ncores An integers. The number of CPU cores used for running
#'   the algorithms (default = 4).
#'
#' @inheritParams sigidentPrognostic
#' @inheritParams validate_diagnostic_signatures
#'
#'
#' @export
sigidentDiagnostic <- function(mergeset, # nolint
                               diagnosis,
                               seed = 111,
                               nfolds = 10,
                               repeats = 5,
                               tunelength = 10,
                               split = 0.8,
                               plotdir = paste0(tempdir(), "/plots/"),
                               colindices = NULL,
                               ncores = 4) {

  stopifnot(
    class(mergeset) == c("matrix", "array"),
    is.numeric(diagnosis),
    is.numeric(seed),
    seed > 0,
    is.numeric(nfolds),
    nfolds > 0,
    is.numeric(split),
    split < 1 & split > 0,
    is.numeric(repeats),
    repeats > 0,
    is.numeric(tunelength),
    tunelength > 0,
    is.numeric(ncores),
    ncores > 0
  )

  # create internal list for storage
  rv <- list()

  # store dirs
  rv$plotdir <- sigident.preproc::clean_path_name(plotdir)

  # create output directories
  dir.create(rv$plotdir)

  # set diagnosis
  rv$diagnosis <- diagnosis

  # store seed, traintest_split
  rv$seed <- seed
  rv$nfolds <- nfolds
  rv$repeats <- repeats
  rv$tunelength <- tunelength
  rv$traintest_split <- split

  rv$ncores <- ncores

  # add mergeset to list
  rv$mergeset <- mergeset

  # identification of Diagnostic Signature
  # first, create training_list
  rv$training_list <-
    create_training_test_split(
      diagnosis = rv$diagnosis,
      mergeset = rv$mergeset,
      split = rv$traintest_split,
      seed = rv$seed
    )

  if (!is.null(colindices)) {
    rv$training_list$train$x <- rv$training_list$train$x[, colindices]
    rv$training_list$test$x <- rv$training_list$test$x[, colindices]
  }

  # Lasso regression
  rv$diagnostic_lasso <- sigident_signature(
    traininglist = rv$training_list,
    type = "lasso",
    nfolds = rv$nfolds,
    seed = rv$seed,
    ncores = rv$ncores
  )
  plot_cvplot(
    cv_obj = rv$diagnostic_lasso$fit_cv,
    filename = paste0(rv$plotdir, "CV_lasso.png")
  )
  plot_rocplot(
    roc = rv$diagnostic_lasso$roc_min,
    filename = paste0(rv$plotdir, "ROC_Lasso.min.png")
  )
  plot_rocplot(
    roc = rv$diagnostic_lasso$roc_1se,
    filename = paste0(rv$plotdir, "ROC_Lasso.1se.png")
  )


  # Elastic net regression
  rv$diagnostic_elasticnet <-
    sigident_signature(
      traininglist = rv$training_list,
      type = "elastic",
      a = 0.9,
      nfolds = rv$nfolds,
      seed = rv$seed,
      ncores = rv$ncores
    )
  plot_cvplot(
    cv_obj = rv$diagnostic_elasticnet$fit_cv,
    filename = paste0(rv$plotdir, "CV_elasticNet.png")
  )
  plot_rocplot(
    roc = rv$diagnostic_elasticnet$roc_min,
    filename = paste0(rv$plotdir, "ROC_elasticNet.min.png")
  )
  plot_rocplot(
    roc = rv$diagnostic_elasticnet$roc_1se,
    filename = paste0(rv$plotdir, "ROC_elasticNet.1se.png")
  )

  # with both calculated hyperparameters alpha and
  # lambda applying grid search
  rv$diagnostic_glmgrid <-
    sigident_signature(
      traininglist = rv$training_list,
      type = "glmnet",
      nfolds = rv$nfolds,
      repeats = rv$repeats,
      tunelength = rv$tunelength,
      seed = rv$seed,
      ncores = rv$ncores
    )
  # plot model of gridsearch
  plot_grid_model_plot(
    model = rv$diagnostic_glmgrid$model,
    filename = paste0(rv$plotdir, "glmnet_gridsearch_model.png")
  )
  # plot variable importance of gridsearch
  plot_grid_varimp_plot(
    model = rv$diagnostic_glmgrid$model,
    filename = paste0(rv$plotdir, "glmnet_gridsearch_variable_importance.png")
  )
  # create roc plot
  plot_rocplot(
    roc = rv$diagnostic_glmgrid$roc,
    filename = paste0(rv$plotdir, "ROC_elasticNet.grid.png")
  )


  # SVM
  rv$diagnostic_svm <-
    sigident_signature(
      traininglist = rv$training_list,
      type = "svm",
      nfolds = rv$nfolds,
      repeats = rv$repeats,
      tunelength = rv$tunelength,
      seed = rv$seed,
      ncores = rv$ncores
    )
  # plot model of gridsearch
  plot_grid_model_plot(
    model = rv$diagnostic_svm$model,
    filename = paste0(rv$plotdir, "SVM_model.png")
  )
  # plot variable importance of gridsearch
  plot_grid_varimp_plot(
    model = rv$diagnostic_svm$model,
    filename = paste0(rv$plotdir, "SVM_variable_importance.png")
  )
  # create roc plot
  plot_rocplot(
    roc = rv$diagnostic_svm$roc,
    filename = paste0(rv$plotdir, "ROC_SVM.png")
  )

  # RF
  rv$diagnostic_rf <-
    sigident_signature(
      traininglist = rv$training_list,
      type = "rf",
      nfolds = rv$nfolds,
      repeats = rv$repeats,
      tunelength = rv$tunelength,
      seed = rv$seed,
      ncores = rv$ncores
    )
  # plot model of gridsearch
  plot_grid_model_plot(
    model = rv$diagnostic_rf$model,
    filename = paste0(rv$plotdir, "RF_model.png")
  )
  # plot variable importance of gridsearch
  plot_grid_varimp_plot(
    model = rv$diagnostic_rf$model,
    filename = paste0(rv$plotdir, "RF_variable_importance.png")
  )
  # create roc plot
  plot_rocplot(
    roc = rv$diagnostic_rf$roc,
    filename = paste0(rv$plotdir, "ROC_RF.png")
  )

  # KNN
  rv$diagnostic_knn <-
    sigident_signature(
      traininglist = rv$training_list,
      type = "knn",
      nfolds = rv$nfolds,
      repeats = rv$repeats,
      tunelength = rv$tunelength,
      seed = rv$seed,
      ncores = rv$ncores
    )
  # plot model of gridsearch
  plot_grid_model_plot(
    model = rv$diagnostic_knn$model,
    filename = paste0(rv$plotdir, "KNN_model.png")
  )
  # plot variable importance of gridsearch
  plot_grid_varimp_plot(
    model = rv$diagnostic_knn$model,
    filename = paste0(rv$plotdir, "KNN_variable_importance.png")
  )
  # create roc plot
  plot_rocplot(
    roc = rv$diagnostic_knn$roc,
    filename = paste0(rv$plotdir, "ROC_KNN.png")
  )

  # GBM
  rv$diagnostic_gbm <-
    sigident_signature(
      traininglist = rv$training_list,
      type = "gbm",
      nfolds = rv$nfolds,
      repeats = rv$repeats,
      tunelength = rv$tunelength,
      seed = rv$seed,
      ncores = rv$ncores
    )
  # plot model of gridsearch
  plot_grid_model_plot(
    model = rv$diagnostic_gbm$model,
    filename = paste0(rv$plotdir, "GBM_model.png")
  )
  # plot variable importance of gridsearch
  plot_grid_varimp_plot(
    model = rv$diagnostic_gbm$model,
    filename = paste0(rv$plotdir, "GBM_variable_importance.png")
  )
  # create roc plot
  plot_rocplot(
    roc = rv$diagnostic_gbm$roc,
    filename = paste0(rv$plotdir, "ROC_GBM.png")
  )

  # compare aucs
  diagnostic_models <- list(

    "lasso" = list(
      "CV" = rv$diagnostic_lasso$fit_cv,
      "min" = list(
        "model" = rv$diagnostic_lasso$lambda_min,
        "confmat" = rv$diagnostic_lasso$confmat_min,
        "prediction" = rv$diagnostic_lasso$prediction_min,
        "auc" = as.numeric(rv$diagnostic_lasso$roc_min$auc)
      ),
      "1se" = list(
        "model" = rv$diagnostic_lasso$lambda_1se,
        "confmat" = rv$diagnostic_lasso$confmat_1se,
        "prediction" = rv$diagnostic_lasso$prediction_1se,
        "auc" = as.numeric(rv$diagnostic_lasso$roc_1se$auc)
      )
    ),

    "elasticnet" = list(
      "CV" = rv$diagnostic_elasticnet$fit_cv,
      "min" = list(
        "model" = rv$diagnostic_elasticnet$lambda_min,
        "confmat" = rv$diagnostic_elasticnet$confmat_min,
        "prediction" = rv$diagnostic_elasticnet$prediction_min,
        "auc" = as.numeric(rv$diagnostic_elasticnet$roc_min$auc)
      ),
      "1se" = list(
        "model" = rv$diagnostic_elasticnet$lambda_1se,
        "confmat" = rv$diagnostic_elasticnet$confmat_1se,
        "prediction" = rv$diagnostic_elasticnet$prediction_1se,
        "auc" = as.numeric(rv$diagnostic_elasticnet$roc_1se$auc)
      )
    ),

    "glmnet" = list(
      "model" = rv$diagnostic_glmgrid$model,
      "confmat" = rv$diagnostic_glmgrid$confmat,
      "prediction" = rv$diagnostic_glmgrid$prediction,
      "auc" = as.numeric(rv$diagnostic_glmgrid$roc$auc)
    ),

    "svm" = list(
      "model" = rv$diagnostic_svm$model,
      "confmat" = rv$diagnostic_svm$confmat,
      "prediction" = rv$diagnostic_svm$prediction,
      "auc" = as.numeric(rv$diagnostic_svm$roc$auc)
    ),

    "rf" = list(
      "model" = rv$diagnostic_rf$model,
      "confmat" = rv$diagnostic_rf$confmat,
      "prediction" = rv$diagnostic_rf$prediction,
      "auc" = as.numeric(rv$diagnostic_rf$roc$auc)
    ),

    "knn" = list(
      "model" = rv$diagnostic_knn$model,
      "confmat" = rv$diagnostic_knn$confmat,
      "prediction" = rv$diagnostic_knn$prediction,
      "auc" = as.numeric(rv$diagnostic_knn$roc$auc)
    ),

    "gbm" = list(
      "model" = rv$diagnostic_gbm$model,
      "confmat" = rv$diagnostic_gbm$confmat,
      "prediction" = rv$diagnostic_gbm$prediction,
      "auc" = as.numeric(rv$diagnostic_gbm$roc$auc)
    )
  )
  return(diagnostic_models)
}
