#' @title svm_classifier
#'
#' @description Helper function to create a svm classifier model
#'
#' @param traininglist A list object containing the training data. The output
#'   of the function `create_training_test_split()`.
#'
#' @param seed Intilization state of random number generator
#'
#' @export
svm_classifier <- function(traininglist, seed, nfolds) {
  # initialize outlist
  outlist <- list()

  outlist$fit_cv <- caret::trainControl(
    method = "repeatedcv", number = nfolds,
    repeats = 5, savePredictions = TRUE, search = "random"
  )

  outlist$model = build_predictive_svm(traininglist$train$x,
                               traininglist$train$y,
                               outlist$fit_cv)


  outlist$importance <- caret::varImp(outlist$model, scale = FALSE)
  outlist$prediction <- predict(outlist$model, traininglist$test$x)
  outlist$confusion_matrix <- caret::confusionMatrix(
   outlist$prediction, as.factor(traininglist$test$y)
  )
  outlist$roc <- calc_roc(outlist$prediction, as.factor(traininglist$test$y))

 return(outlist)
}

#' @title build_predictive_svm
#'
#' @description Function builds a svm classifier model based on the given data.
#'
#' @param train_x The learning data values.
#'
#' @param train_y The learning data classes.
#'
#' @param fit_cv Options for the cross validation
#'
#' @export
build_predictive_svm <- function(train_x,
                                 train_y,
                                 fit_cv) {

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 2, 2, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  model <- caret::train(
    x = train_x,
    y = as.factor(as.factor(train_y)),
    method = "svmLinear",
    trControl = fit_cv,
    preProc = c("center", "scale"),
    tuneLength = 10,
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

  outdat <- as.factor(predict(model, test_x))

  return(outdat)
}
