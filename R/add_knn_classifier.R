#' @title knn_classifier
#'
#' @description Helper function to create a knn classifier model
#'
#' @param traininglist A list object containing the training data. The output
#'   of the function `create_training_test_split()`.
#'
#' @param seed Intilization state of random number generator
#'
#' @export
knn_classifier <- function(traininglist, seed) {
  # initialize outlist
  outlist <- list()

  train_knn <- caret::trainControl(
    method = "repeatedcv", number = 10
  )


  outlist$model = build_predictive_knn(traininglist$train$x,
                                       traininglist$train$y,
                                       train_knn,
                                       4)

  #outlist$importance <- caret::varImp(outlist$model, scale = FALSE)
  outlist$prediction <- predict_knn(outlist$model, traininglist$test$x)
  outlist$confusion_matrix <- caret::confusionMatrix(
    outlist$prediction, as.factor(traininglist$test$y)
  )
  outlist$roc <- calc_roc(outlist$prediction, as.factor(traininglist$test$y))

 return(outlist)
}


#' @title build_predictive_knn
#'
#' @description Function builds a knn classifier model based on the given data.
#'
#' @param train_x The learning data values.
#'
#' @param train_y The learning data classes.
#'
#' @param fit_cv Options for the cross validation
#'
#' @param neighbors Compared neighbors for knn
#'
#' @export
build_predictive_knn <- function(train_x,
                                 train_y,
                                 fit_cv,
                                 neighbors) {

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 2, 2, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  # https://daviddalpiaz.github.io/stat432sp18/supp/knn_class_r.html
  model <- caret::knn3(
    x = train_x,
    y = as.factor(train_y),
    k = neighbors
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
#'
#' @param test_x The data to be classified.
#'
#' @export
predict_knn <- function(model,
                        test_x) {

  outdat <- as.factor(predict(model, test_x, type = "class"))

  return(outdat)
}
