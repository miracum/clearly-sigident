#' @title sigident_signature
#'
#' @description Helper function to conduct lasso or elastic net regualization
#'   for fitting a GLM in order to identify a multi-gene classifier
#'
#' @param traininglist A list object containing the training data. The output
#'   of the function `create_training_test_split()`.
#' @param type A character string. The algorihm used to perform calculations.
#'   Currently implemented are \emph{"elasticnet_grid", "lasso", "elastic",
#'   "svm", "knn"}.
#'
#' @param alpha A numeric between 0 and 1. The elastic net mixing parameter
#'   passed to `glmnet::glmnet()`.
#'
#' @inheritParams sigidentDiagnostic
#'
#' @export
sigident_signature <- function(traininglist,
                               type,
                               alpha = NULL,
                               nfolds = 10,
                               repeats = 3,
                               seed) {

  stopifnot(
    type %in% c("elasticnet_grid",
                "lasso",
                "elastic",
                "svm",
                "knn",
                "rf"),
    is.numeric(nfolds),
    is.numeric(seed),
    is.list(traininglist)
  )

  if (type == "elasticnet_grid") {
    outlist <- glmnet_gridsearch(traininglist, seed, nfolds, repeats)
  } else if (type == "svm") {
    outlist <- svm_classifier(traininglist, seed, nfolds, repeats)
  } else if (type == "knn") {
    outlist <- knn_classifier(traininglist, seed, nfolds, repeats)
  } else if (type == "rf") {
    outlist <- rf_classifier(traininglist, seed, nfolds, repeats)
  } else if (type %in% c("lasso", "elastic")) {
    outlist <- glmnet_classifier(traininglist, type, seed, nfolds)
  }

  return(outlist)
}
