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
#'
#' @inheritParams sigidentPrognostic
#'
#'
#' @export
sigidentDiagnostic <- function(mergeset, # nolint
                               diagnosis,
                               seed = 111,
                               nfolds = 10,
                               split = 0.8,
                               plotdir = "./plots/") {

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
  rv$plotdir <- sigident.preproc::clean_path_name(plotdir)

  # create output directories
  dir.create(rv$plotdir)

  # set diagnosis
  rv$diagnosis <- diagnosis

  # store seed, traintest_split
  rv$seed <- seed
  rv$nfolds <- nfolds
  rv$traintest_split <- split

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


  # Lasso regression
  rv$diagnostic_lasso <- signature(
    traininglist = rv$training_list,
    type = "lasso",
    nfolds = rv$nfolds,
    seed = rv$seed
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
    signature(
      traininglist = rv$training_list,
      type = "elastic",
      alpha = 0.9,
      nfolds = rv$nfolds,
      seed = rv$seed
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


  # SVM
  rv$diagnostic_svm <- signature(
    traininglist = rv$training_list,
    type = "svm",
    nfolds = rv$nfolds,
    seed = rv$seed
  )
  plot_cvplot(
    cv_obj = rv$diagnostic_svm$model,
    filename = paste0(rv$plotdir, "Model_SVM.png")
  )
  plot_rocplot(
    roc = rv$diagnostic_svm$roc,
    filename = paste0(rv$plotdir, "ROC_SVM.png")
  )

  # KNN
  rv$diagnostic_knn <- signature(
    traininglist = rv$training_list,
    type = "knn",
    nfolds = rv$nfolds,
    seed = rv$seed
  )
  plot_cvplot(
    cv_obj = rv$diagnostic_knn$model,
    filename = paste0(rv$plotdir, "Model_KNN.png")
  )
  plot_rocplot(
    roc = rv$diagnostic_knn$roc,
    filename = paste0(rv$plotdir, "ROC_KNN.png")
  )

  # RF
  rv$diagnostic_rf <- signature(
    traininglist = rv$training_list,
    type = "rf",
    nfolds = rv$nfolds,
    seed = rv$seed
  )
  plot_cvplot(
    cv_obj = rv$diagnostic_rf$model,
    filename = paste0(rv$plotdir, "Model_RF.png")
  )
  plot_rocplot(
    roc = rv$diagnostic_rf$roc,
    filename = paste0(rv$plotdir, "ROC_RF.png")
  )

  # with both calculated hyperparameters alpha and
  # lambda applying grid search
  rv$diagnostic_glmgrid <-
    signature(
      traininglist = rv$training_list,
      type = "elasticnet_grid",
      nfolds = rv$nfolds,
      seed = rv$seed
    )
  # plot model of gridsearch
  plot_grid_model_plot(
    model = rv$diagnostic_glmgrid$caret_train,
    filename = paste0(rv$plotdir, "Gridsearch_model.png")
  )
  # plot variable importance of gridsearch
  plot_grid_varimp_plot(
    model = rv$diagnostic_glmgrid$caret_train,
    filename = paste0(rv$plotdir, "Gridsearch_variable_importance.png")
  )
  # create roc plot
  plot_rocplot(
    roc = rv$diagnostic_glmgrid$roc_elasticnet,
    filename = paste0(rv$plotdir, "ROC_elasticNet.grid.png")
  )


  # compare aucs
  diagnostic_models <- list(
    "lasso" = list(
      "CV" = rv$diagnostic_lasso$fit_cv,
      "min" = list(
        "model" = rv$diagnostic_lasso$lambda_min,
        "confmat" = rv$diagnostic_lasso$confmat_min,
        "prediction" = rv$diagnostic_lasso$predicted_min,
        "auc" = as.numeric(rv$diagnostic_lasso$roc_min$auc)
      ),
      "1se" = list(
        "model" = rv$diagnostic_lasso$lambda_1se,
        "confmat" = rv$diagnostic_lasso$confmat_1se,
        "prediction" = rv$diagnostic_lasso$predicted_1se,
        "auc" = as.numeric(rv$diagnostic_lasso$roc_1se$auc)
      )
    ),
    "elasticnet" = list(
      "CV" = rv$diagnostic_elasticnet$fit_cv,
      "min" = list(
        "model" = rv$diagnostic_elasticnet$lambda_min,
        "confmat" = rv$diagnostic_elasticnet$confmat_min,
        "prediction" = rv$diagnostic_elasticnet$predicted_min,
        "auc" = as.numeric(rv$diagnostic_elasticnet$roc_min$auc)
      ),
      "1se" = list(
        "model" = rv$diagnostic_elasticnet$lambda_1se,
        "confmat" = rv$diagnostic_elasticnet$confmat_1se,
        "prediction" = rv$diagnostic_elasticnet$predicted_1se,
        "auc" = as.numeric(rv$diagnostic_elasticnet$roc_1se$auc)
      )
    ),

    "svm" = list(
        "model" = rv$diagnostic_svm$model,
        "confmat" = rv$diagnostic_svm$confusion_matrix,
        "prediction" = rv$diagnostic_svm$prediction,
        "auc" = as.numeric(rv$diagnostic_svm$roc$auc)
    ),

    "knn" = list(
      "model" = rv$diagnostic_knn$model,
      "confmat" = rv$diagnostic_knn$confusion_matrix,
      "prediction" = rv$diagnostic_knn$prediction,
      "auc" = as.numeric(rv$diagnostic_knn$roc$auc)
    ),

    "rf" = list(
      "model" = rv$diagnostic_rf$model,
      "confmat" = rv$diagnostic_rf$confusion_matrix,
      "prediction" = rv$diagnostic_rf$prediction,
      "auc" = as.numeric(rv$diagnostic_rf$roc$auc)
    ),

    "elasticnet_grid" = list(
      "CV" = rv$diagnostic_glmgrid$caret_train,
      "model" = rv$diagnostic_glmgrid$elasticnet_auto,
      "confmat" = rv$diagnostic_glmgrid$confmat_elasticnet,
      "prediction" = rv$diagnostic_glmgrid$predicted_elasticnet,
      "auc" = as.numeric(rv$diagnostic_glmgrid$roc_elasticnet$auc)
    )
  )
  return(diagnostic_models = diagnostic_models)
}
