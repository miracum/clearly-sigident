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

  #train.svm <- caret::trainControl(method = "repeatedcv", number = 10,
  #                                  repeats = 5, savePredictions = TRUE, search = "random")

  outlist$svm_train <- caret::train(
               x = traininglist$train$x,
               y = as.factor(traininglist$train$y),
               method = "svmLinear",
               #family = "binomial",
               #tuneGrid = gr_init$srchGrd,
               trControl = gr_init$trnCtrl, # or train.svm
               preProc = c("center", "scale"),
               tuneLength = 10,
               allowParallel = T)


  importance.svm <- varImp(svm, scale = FALSE)
  print(importance.svm)
  plot(importance.svm, top = 20)

  svm.pred <- predict(svm, test)

  confusion.svm <- confusionMatrix(svm.pred, test$Therapy)
  confusion.svm

  svm.outcome <- test$Therapy


  # stop parallel computation
  parallel::stopCluster(cl)
  gc()

  ########################################

  outlist$svm_auto <-
    perform_glmnet(
      train_x = traininglist$train$x,
      train_y = traininglist$train$y,
      alpha = outlist$svm_train$bestTune$alpha,
      lambda = outlist$svm_train$bestTune$lambda
    )

  # prediction
  pred_svm <- glm_prediction(
    model = outlist$svm_auto,
    test_x = traininglist$test$x,
    test_y = traininglist$test$y,
    s = NULL
  )

  outlist$predicted_svm <- pred_svm$prediction
  outlist$confmat_svm <- pred_svm$confmat
  outlist$roc_svm <- pred_svm$roc

  ########################################

  outlist$elasticnet_auto <-
    perform_glmnet(
      train_x = traininglist$train$x,
      train_y = traininglist$train$y,
      alpha = outlist$caret_train$bestTune$alpha,
      lambda = outlist$caret_train$bestTune$lambda
    )

  # prediction
  pred_elasticnet <- glm_prediction(
    model = outlist$elasticnet_auto,
    test_x = traininglist$test$x,
    test_y = traininglist$test$y,
    s = NULL
  )

  outlist$predicted_elasticnet <- pred_elasticnet$prediction
  outlist$confmat_elasticnet <- pred_elasticnet$confmat
  outlist$roc_elasticnet <- pred_elasticnet$roc

  return(outlist)
}
