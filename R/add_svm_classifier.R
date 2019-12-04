svm_classifier <- function(traininglist, seed) {
  # initialize outlist
  outlist <- list()

  # # initialize gridserach parameters
  # gr_init <- init_grid_search()

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 4, 4, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

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

  train.svm <- caret::trainControl(method = "repeatedcv", number = 10,
                                    repeats = 5, savePredictions = TRUE, search = "random")

  print("step1")

  outlist$svm_train <- caret::train(
               x = traininglist$train$x,
               y = as.factor(traininglist$train$y),
               method = "svmLinear",
               #family = "binomial",
               #tuneGrid = gr_init$srchGrd,
               trControl = train.svm, # or train.svm
               preProc = c("center", "scale"),
               tuneLength = 10,
               allowParallel = T)


  # importance.svm <- varImp(svm, scale = FALSE)
  # print(importance.svm)
  # plot(importance.svm, top = 20)

  # svm.pred <- predict(svm, test)

  # confusion.svm <- confusionMatrix(svm.pred, test$Therapy)
  # confusion.svm

  # svm.outcome <- test$Therapy


  # stop parallel computation
  parallel::stopCluster(cl)
  gc()

  ########################################
}
