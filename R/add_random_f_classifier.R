random_forest <- function(traininglist, seed, nfolds) {
  # initialize outlist
  outlist <- list()

  train_rf <- caret::trainControl(
    method = "repeatedcv", number = nfolds,
    repeats = 3
  )

  #TODO
  library(randomForest)

  metric <- "Accuracy"

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 2, 2, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  outlist$model <- caret::train(
               x = traininglist$train$x,
               y = traininglist$train$y,
               method = "rf",
               trControl = train_rf,
               metric = metric,
               preProc = c("center", "scale"),
               tuneLength = 10,
               allowParallel = T)
  # stop parallel computation
  parallel::stopCluster(cl)
  gc()

  outlist$importance <- caret::varImp(outlist$model, scale = FALSE)
  outlist$predictons <- predict(outlist$model, traininglist$test$x)
  outlist$confusion_matrix <- caret::confusionMatrix(
   outlist$predictons, as.factor(traininglist$test$y)
  )
  #roc


  return(outlist)
}
