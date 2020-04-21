svm_classifier <- function(traininglist, seed, nfolds) {
  # initialize outlist
  outlist <- list()

  train_svm <- caret::trainControl(
    method = "repeatedcv", number = nfolds,
    repeats = 5, savePredictions = TRUE, search = "random"
  )

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 2, 2, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  outlist$svm_model <- caret::train(
               x = traininglist$train$x,
               y = as.factor(as.factor(traininglist$train$y)),
               method = "svmLinear",
               #family = "binomial",
               #tuneGrid = gr_init$srchGrd,
               trControl = train_svm, # or train.svm
               preProc = c("center", "scale"),
               tuneLength = 10,
               allowParallel = T)
  # stop parallel computation
  parallel::stopCluster(cl)
  gc()

  outlist$importance_svm <- caret::varImp(outlist$svm_model, scale = FALSE)
  outlist$svm_pred <- predict(outlist$svm_model, traininglist$test$x)
  outlist$confusion_svm <- caret::confusionMatrix(
   outlist$svm_pred, as.factor(traininglist$test$y)
  )

 return(outlist)
}
