#' @title sigident_signature
#'
#' @description Helper function to identify a multi-gene classifier
#'
#' @param traininglist A list object containing the training data. The output
#'   of the function `create_training_test_split()`.
#' @param type A character string. The algorihm used to perform calculations.
#'   Currently implemented are \emph{"elasticnet_grid", "lasso", "elastic",
#'   "svm" (support vector machine), "rf" (random forest) and "knn" (k-nearest
#'   neighbors)}.
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
                               seed,
                               ncores = 4) {

  stopifnot(
    type %in% c("elasticnet_grid",
                "lasso",
                "elastic",
                "svm",
                "knn",
                "rf"),
    is.list(traininglist),
    is.numeric(seed),
    seed > 0,
    is.numeric(nfolds),
    nfolds > 0,
    is.numeric(repeats),
    repeats > 0,
    is.numeric(tunelength),
    tunelength > 0,
    is.numeric(ncores),
    ncores > 0
  )

  message(paste0("signature identification using type = ", type))

  # go parallel
  available_cores <- parallel::detectCores()

  # reset cores, if user-provided ncores > available_cores
  if (ncores > available_cores) {
    cores <- available_cores
  } else {
    cores <- ncores
  }

  message(paste0("### parallel computation using ", cores, " cores ###"))

  cl <- parallel::makeCluster(cores)
  set.seed(seed)
  doParallel::registerDoParallel(cl)

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

  # stop parallel computation
  parallel::stopCluster(cl)
  gc()

  return(outlist)
}
