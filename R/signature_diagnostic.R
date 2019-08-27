createTrainingTest_ <- function(diagnosis, combat, split = 0.8, seed = 111){

  data.for.calculation <- data.table::data.table(cbind(diagnosis, t(combat)))
  colnames(data.for.calculation)[1] <- "diagnosis"

  # randomly split the data into training set (80% for building a predictive model) and test set (20% for evaluating the model)
  set.seed(seed)
  training.samples <- caret::createDataPartition(y = data.for.calculation[,get("diagnosis")], p = split, list = FALSE)

  train.x  <- as.matrix(data.for.calculation[training.samples, -1])
  test.x <- as.matrix(data.for.calculation[-training.samples, -1])

  y <- data.for.calculation[,factor(get("diagnosis"))]
  train.y <- y[training.samples]
  test.y <- y[-training.samples]

  return(list(train = list(x = train.x, y = train.y),
              test = list(x = test.x, y = test.y)))
}

performGLMnet_ <- function(train.x, train.y, alpha, lambda){
  outdat <- glmnet::glmnet(train.x,
                           train.y,
                           family = "binomial",
                           alpha = alpha,
                           lambda = lambda)
  return(outdat)
}


buildPredictiveGLM_ <- function(train.x, train.y, alpha, fitCV){
  outlist <- list()
  for (i in c("lambda.min", "lambda.1se")){
    outlist[[i]] <- performGLMnet_(train.x, train.y, alpha, fitCV[[i]])
  }
  return(outlist)
}


predictGLM_ <- function(model, test.x, s, type){
  outdat <- as.factor(glmnet::predict.glmnet(
    model,
    newx=test.x,
    s=s,
    type=type))
  return(outdat)
}


calcROC_ <- function(test.y, prediction){
  return(pROC::roc(response = as.numeric(as.character(test.y)), predictor = as.numeric(as.character(prediction))))
}



glmnetSignature_ <- function(traininglist, alpha = 0.9, nfolds = 10, seed){
  # initialize outlist
  outlist <- list()

  set.seed(seed)
  outlist$fitCV <- glmnet::cv.glmnet(traininglist$train$x,
                                     traininglist$train$y,
                                     family = "binomial",
                                     type.measure = "mse",
                                     nfolds = nfolds,
                                     alpha = alpha)

  # build the predictive models utilizing calculated lambda values
  glmpred <- buildPredictiveGLM_(traininglist$train$x,
                                 traininglist$train$y,
                                 alpha = alpha,
                                 fitCV = outlist$fitCV)
  outlist$lambda.min <- glmpred$lambda.min
  outlist$lambda.1se <- glmpred$lambda.1se


  # predict the response variable (diagnosis) of the test data
  # min
  outlist$predicted.min <- predictGLM_(model = outlist$lambda.min,
                                       test.x=traininglist$test$x,
                                       s=outlist$fitCV$lambda.min,
                                       type = "class")
  outlist$confmat.min <- caret::confusionMatrix(data = factor(ifelse(as.numeric(as.character(outlist$predicted.min)) < 0.5, 0, 1)),
                                                reference = traininglist$test$y,
                                                positive = "1") # determine the true case with the 'positive' argument
  # Calculate Roc (min)
  outlist$roc.min <- calcROC_(test.y = traininglist$test$y,
                              prediction = outlist$predicted.min)


  # 1se
  outlist$predicted.1se <- predictGLM_(model = outlist$lambda.1se,
                                       test.x=traininglist$test$x,
                                       s=outlist$fitCV$lambda.1se,
                                       type = "class")
  outlist$confmat.1se <- caret::confusionMatrix(data = factor(ifelse(as.numeric(as.character(outlist$predicted.1se)) < 0.5, 0, 1)),
                                                reference = traininglist$test$y,
                                                positive = "1") # determine the true case with the 'positive' argument
  # Calculate Roc (1se)
  outlist$roc.1se <- calcROC_(test.y = traininglist$test$y,
                              prediction = outlist$predicted.1se)

  return(outlist)
}
