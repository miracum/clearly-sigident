#' @title svm_classifier
#'
#' @description Helper function to create a svm classifier model
#'
#' @param traininglist A list object containing the training data. The output
#'   of the function `create_training_test_split()`.
#' @param seed Intilization state of random number generator
#'
#' @inheritParams sigidentDiagnostic
#'
svm_classifier <- function(traininglist, seed, nfolds) {
  # initialize outlist
  outlist <- list()

  fit_cv <- caret::trainControl(
    method = "repeatedcv",
    number = nfolds,
    repeats = 5,
    savePredictions = TRUE,
    search = "random"
  )

  outlist$svm_model <- build_predictive_svm(
    train_x = traininglist$train$x,
    train_y = traininglist$train$y,
    fit_cv = fit_cv
  )

  outlist$predicted_svm <- predict_svm(
    model = outlist$svm_model,
    test_x = traininglist$test$x
  )
  outlist$confmat_svm <- caret::confusionMatrix(
    data = outlist$predicted_svm,
    reference = traininglist$test$y,
    positive = "1"
  )
  outlist$roc_svm <- calc_roc(
    test_y = traininglist$test$y,
    prediction = outlist$predicted_svm
  )

  return(outlist)
}

#' @title build_predictive_svm
#'
#' @description Function builds a svm classifier model based on the given data.
#'
#' @param train_x The learning data values.
#' @param train_y The learning data classes.
#' @param fit_cv Options for the cross validation
#'
#' @export
build_predictive_svm <- function(train_x,
                                 train_y,
                                 fit_cv) {

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 4, 4, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  model <- caret::train(
    x = train_x,
    y = as.factor(as.factor(train_y)),
    method = "svmLinear",
    trControl = fit_cv,
    preProc = c("center", "scale"),
    tuneLength = 5,
    allowParallel = T)

  # stop parallel computation
  parallel::stopCluster(cl)
  gc()


  return(model)
}

#' @title predict_svm
#'
#' @description Function to classify given data based on the created
#'   model.
#'
#' @param model A list object containing the training data. The output
#'   of the function `build_predictive_svm()`.
#'
#' @param test_x The data to be classified.
#'
#' @export
predict_svm <- function(model,
                        test_x) {

  outdat <- caret::predict.train(
    model,
    test_x,
    type = "raw"
  )

  return(outdat)
}
