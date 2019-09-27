#' @title createTrainingTest_
#'
#' @description Helper function to split data into training and test set
#'
#' @inheritParams sigidentDEG
#' @inheritParams sigidentDiagnostic
#'
#' @export
createTrainingTest_ <- function(diagnosis, mergeset, split = 0.8, seed = 111){

  data.for.calculation <- as.data.frame(cbind(diagnosis, t(mergeset)))

  # randomly split the data into training set (80% for building a predictive model) and test set (20% for evaluating the model)
  set.seed(seed)
  training.samples <- caret::createDataPartition(y = data.for.calculation$diagnosis, p = split, list = FALSE)

  train.x  <- as.matrix(data.for.calculation[training.samples, -1])
  test.x <- as.matrix(data.for.calculation[-training.samples, -1])

  y <- factor(data.for.calculation$diagnosis)
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
  gc()
  return(outdat)
}


buildPredictiveGLM_ <- function(train.x, train.y, alpha, fitCV){
  outlist <- list()
  for (i in c("lambda.min", "lambda.1se")){
    outlist[[i]] <- performGLMnet_(train.x, train.y, alpha, fitCV[[i]])
  }
  gc()
  return(outlist)
}


predictGLM_ <- function(model, test.x, s = NULL, type){
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

glmPrediction_ <- function(model, test.x, test.y, s = NULL){

  # initialize outlist
  outlist <- list()

  # Calculate prediction
  outlist$predicted <- predictGLM_(model = model,
                                   test.x = test.x,
                                   s = s,
                                   type = "class")

  # Generate Confusion Matrix
  outlist$confmat <- caret::confusionMatrix(data = factor(ifelse(as.numeric(as.character(outlist$predicted)) < 0.5, 0, 1)),
                                            reference = test.y,
                                            positive = "1") # determine the true case with the 'positive' argument
  # Calculate Roc
  outlist$roc <- calcROC_(test.y = test.y,
                          prediction = outlist$predicted)

  return(outlist)
}


#' @title signature_
#'
#' @description Helper function to conduct lasso or elastic net regualization for fitting a GLM in order to identify a multi-gene classifier
#'
#' @param traininglist A list object containing the training data. The output of the function `createTrainingTest_()`.
#' @param type A character string. The algorihm used to perform calculations. Currently implemented are \emph{"grid", "lasso", "elastic"}.
#'
#' @param alpha A numeric between 0 and 1. The elastic net mixing parameter passed to `glmnet::glmnet()`.
#'
#' @inheritParams sigidentDiagnostic
#'
#' @export
signature_ <- function(traininglist, type, alpha = NULL, nfolds = 10, seed){

  stopifnot(
    type %in% c("grid", "lasso", "elastic"),
    is.numeric(nfolds),
    is.numeric(seed),
    is.list(traininglist)
  )

  if (type == "grid"){
    outlist <- glmnetGridSearch_(traininglist, seed)
  } else {

    # use provided alpha only in elastic
    if (type == "lasso"){
      alpha <- 1
    } else if (type == "elastic"){
      stopifnot(
        !is.null(alpha),
        alpha >= 0 | alpha <= 1
      )
    }

    # initialize outlist
    outlist <- list()


    # go parallel
    ncores <- parallel::detectCores()
    cores <- ifelse(ncores > 4, 4, ncores)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)

    set.seed(seed)
    outlist$fitCV <- glmnet::cv.glmnet(traininglist$train$x,
                                       traininglist$train$y,
                                       family = "binomial",
                                       type.measure = "mse",
                                       nfolds = nfolds,
                                       alpha = alpha,
                                       parallel = TRUE)
    # stop parallel computation
    parallel::stopCluster(cl)
    gc()

    # build the predictive models utilizing calculated lambda values
    glmpred <- buildPredictiveGLM_(traininglist$train$x,
                                   traininglist$train$y,
                                   alpha = alpha,
                                   fitCV = outlist$fitCV)
    outlist$lambda.min <- glmpred$lambda.min
    outlist$lambda.1se <- glmpred$lambda.1se


    # predict the response variable (diagnosis) of the test data
    # min
    pred.min <- glmPrediction_(model = outlist$lambda.min,
                               test.x = traininglist$test$x,
                               test.y = traininglist$test$y,
                               s = outlist$fitCV$lambda.min)

    outlist$predicted.min <- pred.min$predicted
    outlist$confmat.min <- pred.min$confmat
    outlist$roc.min <- pred.min$roc


    # 1se
    pred.1se <- glmPrediction_(model = outlist$lambda.1se,
                               test.x = traininglist$test$x,
                               test.y = traininglist$test$y,
                               s = outlist$fitCV$lambda.1se)

    outlist$predicted.1se <- pred.1se$predicted
    outlist$confmat.1se <- pred.1se$confmat
    outlist$roc.1se <- pred.1se$roc
  }

  return(outlist)
}


initGridSearch_ <- function(){
  # set up alpha and lambda grid to search for pair that minimizes CV erros
  lambda.grid <- 10^seq(0, -4, length = 100)
  alpha.grid <- seq(0.1, 1, length=10)

  # set up cross validation method for train function
  trnCtrl=caret::trainControl(
    method = "repeatedcv",
    number=10
  )

  srchGrd=expand.grid(
    .alpha = alpha.grid,
    .lambda = lambda.grid
  )

  # set up search grid for alpha and lambda parameters
  return(list(srchGrd=srchGrd, trnCtrl=trnCtrl))
}


glmnetGridSearch_ <- function(traininglist, seed){

  # initialize outlist
  outlist <- list()

  # initialize gridserach parameters
  gr.init <- initGridSearch_()

  # go parallel
  ncores <- parallel::detectCores()
  cores <- ifelse(ncores > 4, 4, ncores)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  # perform cv forecasting
  set.seed(seed)
  outlist$caret.train <- caret::train(x = traininglist$train$x,
                                      y = as.factor(traininglist$train$y),
                                      method = "glmnet",
                                      family = "binomial",
                                      tuneGrid = gr.init$srchGrd,
                                      trControl = gr.init$trnCtrl,
                                      standardize = TRUE,
                                      maxit = 1e7)
  # stop parallel computation
  parallel::stopCluster(cl)
  gc()

  outlist$elasticNet.auto <- performGLMnet_(train.x = traininglist$train$x,
                                            train.y = traininglist$train$y,
                                            alpha = outlist$caret.train$bestTune$alpha,
                                            lambda = outlist$caret.train$bestTune$lambda)

  # prediction
  pred.elasticNet <- glmPrediction_(model = outlist$elasticNet.auto,
                                    test.x = traininglist$test$x,
                                    test.y = traininglist$test$y,
                                    s = NULL)

  outlist$predicted.elasticNet <- pred.elasticNet$prediction
  outlist$confmat.elasticNet <- pred.elasticNet$confmat
  outlist$roc.elasticNet <- pred.elasticNet$roc

  return(outlist)
}

#' @title geneMapSig_
#'
#' @description Helper function to map relevant input variables of a diagnostic model to corresponding IDs.
#'
#' @inheritParams createGridModelPlot_
#' @inheritParams sigidentDEG
#'
#' @export
geneMapSig_ <- function(mergeset, model){
  id <- rownames(mergeset)
  # TODO warum i+1?
  index <- model[["beta"]]@i+1
  # TODO map entrez_id on gene symbol here and include as second columen to ouput
  return(as.data.frame(x = cbind("ID" = id[index])))
}

#' @title validateDiagnosticSignature_
#'
#' @description Helper function to validate diagnostic signatures
#'
#' @param validationstudylist A list, containing metainformation on the study used for validation of the diagnostic signature.
#' @param models A prediction model to be used for validation.
#' @param datadir A character string. Path to the data-folder inside the metadata folder.
#'
#' @inheritParams sigidentDEG
#' @inheritParams createDEGheatmap_
#' @inheritParams batchCorrection_
#' @inheritParams createDiagnosisDesignBatch_
#'
#' @export
validateDiagnosticSignature_ <- function(validationstudylist,
                                         models,
                                         genes,
                                         idtype,
                                         targetname,
                                         controlname,
                                         targetcol,
                                         datadir){
  stopifnot(
    is.character(targetcol),
    is.character(targetname),
    is.character(controlname),
    is.list(validationstudylist),
    is.character(validationstudylist$studyname),
    is.character(validationstudylist$targetcolname),
    is.character(validationstudylist$targetlevelname),
    is.character(validationstudylist$controllevelname),
    is.numeric(validationstudylist$setid)
  )

  outlist <- list()

  eset <- loadEset_(name = validationstudylist$studyname,
                    datadir = datadir,
                    targetcolname = validationstudylist$targetcolname,
                    targetlevelname = validationstudylist$targetlevelname,
                    controllevelname = validationstudylist$controllevelname,
                    targetcol = targetcol,
                    targetname = targetname,
                    controlname = controlname,
                    setid = validationstudylist$setid)



  diagnosis <- diagnosis_(vector = eset[[targetcol]],
                          targetname = targetname,
                          controlname = controlname)


  expr <- createExpressionSet_(eset = eset,
                               idtype = idtype)

  # TODO why is this code in the original script?
  # # creating data frame, selecting only DEGs
  # v.data.DEG <- base::subset(t(expr), select = genes)

  # creating data frame including all genes
  v.data.all <- t(expr)

  for (i in names(models)){

    if (i %in% c("lasso", "elasticnet")){

      for (j in c("min", "1se")){

        predicted <- predict(model[[i]][[j]]$model, v.data.all, type = "response")

        confmat <- caret::confusionMatrix(data = factor(ifelse(as.numeric(as.character(predicted)) < 0.5, 0, 1)),
                                          reference = factor(diagnosis),
                                          positive = "1") # determine the true case with the 'positive' argument

        # Calculate Roc
        roc <- calcROC_(test.y = diagnosis,
                        prediction = predicted)

        outlist[[i]][[j]] <- list(predicted = predicted, confmat = confmat, roc = roc)
      }
    }
  }
  return(outlist)
}
