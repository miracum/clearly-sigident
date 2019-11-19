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
  set.seed(seed)
  outlist$caret_train <- caret::train(
    x = traininglist$train$x,
    y = as.factor(traininglist$train$y),
    method = "glmnet",
    family = "binomial",
    tuneGrid = gr_init$srchGrd,
    trControl = gr_init$trnCtrl,
    standardize = TRUE,
    maxit = 1e7
  )
  # stop parallel computation
  parallel::stopCluster(cl)
  gc()

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
