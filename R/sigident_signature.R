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
#' @param a A numeric between 0 and 1. The elastic net mixing parameter 'alpha'
#'   passed to `glmnet::glmnet()`.
#'
#' @inheritParams sigidentDiagnostic
#'
#' @export
sigident_signature <- function(traininglist,
                               type,
                               a = NULL,
                               nfolds = 10,
                               repeats = 5,
                               tunelength = 10,
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
    outlist <- glmnet_gridsearch(
      traininglist = traininglist,
      seed = seed,
      nfolds = nfolds,
      repeats = repeats,
      tunelength = tunelength
    )
  } else if (type == "svm") {
    outlist <- svm_classifier(
      traininglist = traininglist,
      seed = seed,
      nfolds = nfolds,
      repeats = repeats,
      tunelength = tunelength
    )
  } else if (type == "knn") {
    outlist <- knn_classifier(
      traininglist = traininglist,
      seed = seed,
      nfolds = nfolds,
      repeats = repeats,
      tunelength = tunelength
    )
  } else if (type == "rf") {
    outlist <- rf_classifier(
      traininglist = traininglist,
      seed = seed,
      nfolds = nfolds,
      repeats = repeats,
      tunelength = tunelength
    )
  } else if (type %in% c("lasso", "elastic")) {
    outlist <- glmnet_classifier(traininglist, type, seed, nfolds, a)
  }

  return(outlist)
}
