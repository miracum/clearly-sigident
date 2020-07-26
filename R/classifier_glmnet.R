perform_glmnet <- function(train_x,
                           train_y,
                           alpha,
                           lambda,
                           seed) {

  set.seed(seed)
  outdat <- glmnet::glmnet(
    x = train_x,
    y = train_y,
    family = "binomial",
    alpha = alpha,
    lambda = lambda
  )
  gc()
  return(outdat)
}


build_predictive_glm <- function(train_x,
                                 train_y,
                                 alpha,
                                 fit_cv) {

  outlist <- list()
  for (i in c("lambda.min", "lambda.1se")) {
    outlist[[i]] <- perform_glmnet(
      train_x = train_x,
      train_y = train_y,
      alpha = alpha,
      lambda = fit_cv[[i]],
      seed = seed
    )
  }
  gc()
  return(outlist)
}


predict_glm <- function(model,
                        test_x,
                        s = NULL,
                        type) {

  outdat <- as.factor(
    glmnet::predict.glmnet(
      object = model,
      newx = test_x,
      s = s,
      type = type
    )
  )
  return(outdat)
}

glm_prediction <- function(model,
                           test_x,
                           test_y,
                           s = NULL) {

  # initialize outlist
  outlist <- list()

  # Calculate prediction
  outlist$predicted <- predict_glm(
    model = model,
    test_x = test_x,
    s = s,
    type = "class"
  )

  # Generate Confusion Matrix
  outlist$confmat <-
    caret::confusionMatrix(
      data = factor(ifelse(as.numeric(
        as.character(outlist$predicted)
      ) < 0.5, 0, 1)),
      reference = test_y,
      positive = "1"
    ) # determine the true case with the 'positive' argument
  # Calculate Roc
  outlist$roc <- calc_roc(test_y = test_y,
                          prediction = outlist$predicted)

  return(outlist)
}

glmnet_classifier <- function(traininglist, type, seed, nfolds, a) {
  # use provided alpha only in elastic
  if (type == "lasso") {
    a <- 1
  } else if (type == "elastic") {
    stopifnot(
      !is.null(a),
      a >= 0 || a <= 1
    )
  }

  # initialize outlist
  outlist <- list()

  set.seed(seed)
  outlist$fit_cv <- glmnet::cv.glmnet(
    x = traininglist$train$x,
    y = traininglist$train$y,
    family = "binomial",
    type.measure = "mse",
    nfolds = nfolds,
    alpha = a,
    parallel = TRUE
  )

  # build the predictive models utilizing calculated lambda values
  glmpred <- build_predictive_glm(
    train_x = traininglist$train$x,
    train_y = traininglist$train$y,
    alpha = a,
    fit_cv = outlist$fit_cv
  )
  outlist$lambda_min <- glmpred$lambda.min
  outlist$lambda_1se <- glmpred$lambda.1se


  # predict the response variable (diagnosis) of the test data
  # min
  pred_min <- glm_prediction(
    model = outlist$lambda_min,
    test_x = traininglist$test$x,
    test_y = traininglist$test$y,
    s = outlist$fit_cv$lambda.min
  )

  outlist$prediction_min <- pred_min$predicted
  outlist$confmat_min <- pred_min$confmat
  outlist$roc_min <- pred_min$roc


  # 1se
  pred_1se <- glm_prediction(
    model = outlist$lambda_1se,
    test_x = traininglist$test$x,
    test_y = traininglist$test$y,
    s = outlist$fit_cv$lambda.1se
  )

  outlist$prediction_1se <- pred_1se$predicted
  outlist$confmat_1se <- pred_1se$confmat
  outlist$roc_1se <- pred_1se$roc

  return(outlist)
}
