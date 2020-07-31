glmnet_gridsearch <- function(
  traininglist,
  seed,
  nfolds,
  repeats,
  tunelength
) {

  stopifnot(
    unique(traininglist$train$y) %in% c(0, 1)
  )

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

  set.seed(seed)
  outlist$model <- caret::train(
    x = traininglist$train$x,
    y = paste0("X", traininglist$train$y),
    method = "glmnet",
    family = "binomial",
    tuneGrid = gr_init$srchGrd,
    trControl = gr_init$trnCtrl,
    standardize = TRUE,
    tuneLength = tunelength,
    maxit = 1e7
  )

  outlist$prediction <- predict_caret(
    model = outlist$model,
    test_x = traininglist$test$x
  )
  outlist$confmat <- caret::confusionMatrix(
    data = factor(ifelse(as.numeric(
      as.character(outlist$prediction)
    ) < 0.5, 0, 1)),
    reference = traininglist$test$y,
    positive = "1"
  )
  outlist$roc <- calc_roc(
    test_y = traininglist$test$y,
    prediction = outlist$prediction
  )

  return(outlist)
}
