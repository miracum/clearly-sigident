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
rf_classifier <- function(traininglist, seed, nfolds, repeats) {

  stopifnot(
    unique(traininglist$train$y) %in% c(0, 1)
  )

  # initialize outlist
  outlist <- list()

  trn_ctrl <- caret::trainControl(
    method = "repeatedcv",
    number = nfolds,
    repeats = repeats,
    classProbs = TRUE
  )

  outlist$model <- build_predictive_rf(
    train_x = traininglist$train$x,
    train_y = paste0("X", traininglist$train$y),
    trn_ctrl = trn_ctrl
  )

  outlist$prediction <- predict_caret(
    model = outlist$model,
    test_x = traininglist$test$x
  )
  outlist$confmat <- caret::confusionMatrix(
    data = factor(ifelse(as.numeric(
      as.character(outlist$prediction)
    ) < 0.5, 0, 1)),
    reference = traininglist$test$y,
    positive = "1"
  )
  outlist$roc <- calc_roc(
    test_y = traininglist$test$y,
    prediction = outlist$prediction
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
#' @param trn_ctrl Options for the cross validation.
#'
build_predictive_rf <- function(train_x,
                                train_y,
                                trn_ctrl) {

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 4, 4, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  model <- caret::train(
    x = train_x,
    y = as.factor(train_y),
    method = "rf",
    trControl = trn_ctrl,
    preProc = c("center", "scale"),
    tuneLength = 5,
    allowParallel = T
  )
  # stop parallel computation
  parallel::stopCluster(cl)
  gc()


  return(model)
}

#' @title predict_caret
#'
#' @description Function to classify given data based on the created
#'   model.
#'
#' @param model The caret model.
#' @param test_x The data to be classified.
#'
predict_caret <- function(model,
                       test_x) {

  outdat <- caret::predict.train(
    model,
    test_x,
    type = "prob"
  )[, "X1"]

  return(outdat)
}
