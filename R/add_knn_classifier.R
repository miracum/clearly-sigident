kknn_classifier <- function(traininglist, seed) {
  # initialize outlist
  outlist <- list()
  train_kknn <- caret::trainControl(
    method = "repeatedcv", number = 10,
    repeats = 5, savePredictions = TRUE, search = "random"
  )

  print("step1")

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 2, 2, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  print("step2")

  outlist$model <- caret::train(
               x = traininglist$train$x,
               y = as.factor(traininglist$train$y),
               # wts = 0,
               method = "kknn",
               #family = "binomial",
               #tuneGrid = gr_init$srchGrd,
               trctrl = train_kknn, 
               preProc = c("center", "scale"),
               tuneLength = 10,
               allowParallel = T)
  # stop parallel computation
  parallel::stopCluster(cl)
  gc()

  print(outlist$kknn_model)

  outlist$importance_kknn <- caret::varImp(outlist$kknn_model, scale = FALSE)

  outlist$kknn_pred <- predict(outlist$kknn_model, traininglist$test$x)


 outlist$confusion_kknn <- caret::confusionMatrix(
 outlist$kknn_pred, as.factor(traininglist$test$y)
 )

 return(outlist)
}
