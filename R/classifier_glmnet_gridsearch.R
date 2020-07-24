glmnet_gridsearch <- function(traininglist, seed, nfolds, repeats) {
  # initialize outlist
  outlist <- list()


  # set up search grid for alpha and lambda parameters
  # initialize gridserach parameters
  # set up alpha and lambda grid to search for pair that minimizes CV erros
  lambda_grid <- 10 ^ seq(0, -4, length = 100)
  alpha_grid <- seq(0.1, 1, length = 10)

  # set up cross validation method for train function
  trn_ctrl <- caret::trainControl(
    method = "repeatedcv",
    number = nfolds,
    repeats = repeats,
    classProbs = TRUE
  )

  srch_grd <- expand.grid(.alpha = alpha_grid,
                          .lambda = lambda_grid)

  gr_init <- list(srchGrd = srch_grd, trnCtrl = trn_ctrl)

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

  outlist$model <-
    perform_glmnet(
      train_x = traininglist$train$x,
      train_y = traininglist$train$y,
      alpha = outlist$caret_train$bestTune$alpha,
      lambda = outlist$caret_train$bestTune$lambda
    )

  # prediction
  pred_elasticnet <- glm_prediction(
    model = outlist$model,
    test_x = traininglist$test$x,
    test_y = traininglist$test$y,
    s = NULL
  )

  outlist$prediction <- pred_elasticnet$predicted
  outlist$confmat <- pred_elasticnet$confmat
  outlist$roc <- pred_elasticnet$roc

  return(outlist)
}
