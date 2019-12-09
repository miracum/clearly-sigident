svm_classifier <- function(traininglist, seed) {
  # initialize outlist
  outlist <- list()

  # # initialize gridserach parameters
  # gr_init <- init_grid_search()

  # perform cv forecasting
  #set.seed(seed)
  #outlist$caret_train <- caret::train(
  #  x = traininglist$train$x,
  #  y = as.factor(traininglist$train$y),
  #  method = "glmnet",
  #  family = "binomial",
  #  tuneGrid = gr_init$srchGrd,
  #  trControl = gr_init$trnCtrl,
  #  standardize = TRUE,
  #  maxit = 1e7
  #)

  train_svm <- caret::trainControl(
    method = "repeatedcv", number = 10,
    repeats = 5, savePredictions = TRUE, search = "random"
  )

  print("step1")

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 4, 4, ncores)
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

  print(outlist$svm_model)

  outlist$importance_svm <- caret::varImp(outlist$svm_model, scale = FALSE)

  outlist$svm_pred <- predict(outlist$svm_model, traininglist$test$x)


 outlist$confusion_svm <- caret::confusionMatrix(
   outlist$svm_pred, as.factor(traininglist$test$y)
   )

 return(outlist)
}
