kknn_classifier <- function(traininglist, seed) {
  # initialize outlist
  outlist <- list()

  train_kknn <- caret::trainControl(
    method = "repeatedcv", number = 10
  )
  #   repeats = 5, savePredictions = TRUE, search = "random"
  # )

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 2, 2, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  # https://daviddalpiaz.github.io/stat432sp18/supp/knn_class_r.html
  # outlist$model <- caret::knn3(
  #   x = traininglist$train$x,
  #   y = as.factor(traininglist$train$y),
  #   k = 25
  # )

  outlist$model <- caret::train(
    x = traininglist$train$x,
    y = as.factor(traininglist$train$y),
    method = "kknn",
    #family = "binomial",
    #tuneGrid = gr_init$srchGrd,
    trControl = train_kknn,
    preProcess = c("center", "scale"),
    tuneLength = 10
  )

  # stop parallel computation
  parallel::stopCluster(cl)
  gc()

  outlist$importance <- caret::varImp(outlist$model, scale = FALSE)

  outlist$pred <- predict(outlist$model, traininglist$test$x)


# prediction for knn3
 # outlist$confusion_kknn <- caret::confusionMatrix(
 # as.factor(ifelse(outlist$pred[,1] > 0.5, 1, 0)), as.factor(traininglist$test$y)
 # )
  outlist$confusion_kknn <- caret::confusionMatrix(
    outlist$pred, as.factor(traininglist$test$y)
  )

 return(outlist)
}
