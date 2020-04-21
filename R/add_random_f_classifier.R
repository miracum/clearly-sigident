random_forest <- function(traininglist, seed, nfolds) {
  # initialize outlist
  outlist <- list()

  train_rf <- caret::trainControl(
    method = "repeatedcv", number = nfolds,
    repeats = 3
  )
  metric <- "Accuracy"

  outlist$model <- build_predictive_rf(traininglist$train$x,
                                       traininglist$train$y,
                                       train_rf,
                                       metric
                                       )

  outlist$importance <- caret::varImp(outlist$model, scale = FALSE)
  outlist$prediction <- predict(outlist$model, traininglist$test$x)
  outlist$confusion_matrix <- caret::confusionMatrix(
   outlist$prediction, as.factor(traininglist$test$y)
  )
  outlist$roc <- calc_roc(outlist$prediction, as.factor(traininglist$test$y))

  return(outlist)
 }


#' @title build_predictive_rf
#'
#' @description Function builds a random forest classifier model based on
#' the given data.
#'
#' @param train_x The learning data values.
#'
#' @param train_y The learning data classes.
#'
#' @param fit_cv Options for the cross validation
#'
#' @export
build_predictive_rf <- function(train_x,
                                train_y,
                                train_rf,
                                metric) {

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 2, 2, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  model <- caret::train(
    x = train_x,
    y = train_y,
    method = "rf",
    trControl = train_rf,
    metric = metric,
    preProc = c("center", "scale"),
    tuneLength = 10,
    allowParallel = T)
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
#' @export
predict_rf <- function(model,
                        test_x) {

  outdat <- as.factor(predict(model, test_x, type = "class"))

  return(outdat)
}

