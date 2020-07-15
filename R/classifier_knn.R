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

  outlist$fit_cv <- caret::trainControl(
    method = "repeatedcv",
    repeats = 3
  )

  outlist$knn_model <- build_predictive_knn(
    traininglist$train$x,
    traininglist$train$y,
    outlist$fit_cv
  )

  outlist$predicted_knn <- predict_knn(
    model = outlist$knn_model,
    test_x = traininglist$test$x
  )
  outlist$confmat_knn <- caret::confusionMatrix(
    outlist$predicted_knn,
    as.factor(traininglist$test$y)
  )
  outlist$roc_knn <- calc_roc(
    outlist$predicted_knn,
    as.factor(traininglist$test$y)
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
  cores <- ifelse(ncores > 2, ncores - 1, ncores)
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

  outdat <- as.factor(
    stats::predict(
      model,
      test_x,
      type = "class"
    )
  )

  return(outdat)
}
