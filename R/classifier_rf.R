#' @title rf_classifier
#'
#' @description Helper function to create a rf classifier model
#'
#' @param traininglist A list object containing the training data. The output
#'   of the function `create_training_test_split()`.
#' @param seed Intilization state of random number generator
#'
#' @inheritParams sigidentDiagnostic
#'
rf_classifier <- function(traininglist, seed, nfolds) {
  # initialize outlist
  outlist <- list()

  fit_cv <- caret::trainControl(
    method = "repeatedcv",
    number = nfolds,
    repeats = 5
  )
  metric <- "Accuracy"

  outlist$rf_model <- build_predictive_rf(
    train_x = traininglist$train$x,
    train_y = traininglist$train$y,
    fit_cv = fit_cv,
    metric = metric
  )

  outlist$predicted_rf <- predict_rf(
    model = outlist$rf_model,
    test_x = traininglist$test$x
  )
  outlist$confmat_rf <- caret::confusionMatrix(
    data = outlist$predicted_rf,
    reference = traininglist$test$y,
    positive = "1"
  )
  outlist$roc_rf <- calc_roc(
    test_y = traininglist$test$y,
    prediction = outlist$predicted_rf
  )

  return(outlist)
}


#' @title build_predictive_rf
#'
#' @description Function builds a random forest classifier model based on
#' the given data.
#'
#' @param train_x The learning data values.
#' @param train_y The learning data classes.
#' @param fit_cv Options for the cross validation
#' @param metric A character. The metric to be used during caret::train.
#'
build_predictive_rf <- function(train_x,
                                train_y,
                                fit_cv,
                                metric) {

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 4, 4, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  model <- caret::train(
    x = train_x,
    y = train_y,
    method = "rf",
    trControl = fit_cv,
    metric = metric,
    preProc = c("center", "scale"),
    tuneLength = 5,
    allowParallel = T
  )
  # stop parallel computation
  parallel::stopCluster(cl)
  gc()


  return(model)
}

#' @title predict_rf
#'
#' @description Function to classify given data based on the created
#'   model.
#'
#' @param model A list object containing the training data. The output
#'   of the function `build_predictive_rf()`.
#'
#' @param test_x The data to be classified.
#'
predict_rf <- function(model,
                       test_x) {

  outdat <- caret::predict.train(
    model,
    test_x,
    type = "raw"
  )

  return(outdat)
}
