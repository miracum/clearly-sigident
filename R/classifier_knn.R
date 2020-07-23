#' @title knn_classifier
#'
#' @description Helper function to create a knn classifier model
#'
#' @param traininglist A list object containing the training data. The output
#'   of the function `create_training_test_split()`.
#' @param seed Intilization state of random number generator
#'
#' @inheritParams sigidentDiagnostic
knn_classifier <- function(traininglist, seed, nfolds) {
  # initialize outlist
  outlist <- list()

  fit_cv <- caret::trainControl(
    method = "repeatedcv",
    repeats = 5
  )

  outlist$knn_model <- build_predictive_knn(
    train_x = traininglist$train$x,
    train_y = traininglist$train$y,
    fit_cv = fit_cv
  )

  outlist$predicted_knn <- predict_knn(
    model = outlist$knn_model,
    test_x = traininglist$test$x
  )
  outlist$confmat_knn <- caret::confusionMatrix(
    data = outlist$predicted_knn,
    reference = traininglist$test$y,
    positive = "1"
  )
  outlist$roc_knn <- calc_roc(
    test_y = traininglist$test$y,
    prediction = outlist$predicted_knn
  )

  return(outlist)
}

#' @title build_predictive_knn
#'
#' @description Function builds a knn classifier model based on the given data.
#'
#' @param train_x The learning data values.
#' @param train_y The learning data classes.
#' @param fit_cv Options for the cross validation
build_predictive_knn <- function(train_x,
                                 train_y,
                                 fit_cv
) {

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 4, 4, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  model <- caret::train(
    x = train_x,
    y = as.factor(as.factor(train_y)),
    method = "knn",
    trControl = fit_cv,
    preProc = c("center", "scale"),
    tuneLength = 5
  )


  # stop parallel computation
  parallel::stopCluster(cl)
  gc()


  return(model)
}

#' @title predict_knn
#'
#' @description Function to classify given data based on the created
#'   model.
#'
#' @param model A list object containing the training data. The output
#'   of the function `build_predictive_knn()`.
#' @param test_x The data to be classified.
#'
predict_knn <- function(model,
                        test_x) {

  outdat <- caret::predict.train(
    model,
    test_x,
    type = "raw"
  )

  return(outdat)
}
