#' @title create_training_test_split
#'
#' @description Helper function to split data into training and test set
#'
#' @inheritParams sigidentDEG
#' @inheritParams sigidentDiagnostic
#'
#' @export
create_training_test_split <- function(diagnosis,
                                       mergeset,
                                       split = 0.8,
                                       seed = 111) {

  data_for_calculation <- as.data.frame(cbind(diagnosis, t(mergeset)))

  # randomly split the data into training set (80% for building a predictive
  # model) and test set (20% for evaluating the model)
  set.seed(seed)
  training_samples <- caret::createDataPartition(
    y = data_for_calculation$diagnosis,
    p = split,
    list = FALSE
  )

  train_x  <- as.matrix(data_for_calculation[training_samples, -1])
  test_x <- as.matrix(data_for_calculation[-training_samples, -1])

  y <- factor(data_for_calculation$diagnosis)
  train_y <- y[training_samples]
  test_y <- y[-training_samples]

  return(list(
    train = list(x = train_x,
                 y = train_y),
    test = list(x = test_x,
                y = test_y)
  ))
}

perform_glmnet <- function(train_x,
                           train_y,
                           alpha,
                           lambda) {

  outdat <- glmnet::glmnet(
    train_x,
    train_y,
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
    outlist[[i]] <- perform_glmnet(train_x,
                                   train_y,
                                   alpha,
                                   fit_cv[[i]])
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
      model,
      newx = test_x,
      s = s,
      type = type
    )
  )
  return(outdat)
}


calc_roc <- function(test_y,
                     prediction) {

  return(pROC::roc(
    response = as.numeric(as.character(test_y)),
    predictor = as.numeric(as.character(prediction))
  ))
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
    caret::confusionMatrix(data = factor(ifelse(as.numeric(
      as.character(outlist$predicted)
    ) < 0.5, 0, 1)),
    reference = test_y,
    positive = "1") # determine the true case with the 'positive' argument
  # Calculate Roc
  outlist$roc <- calc_roc(test_y = test_y,
                          prediction = outlist$predicted)

  return(outlist)
}


#' @title signature
#'
#' @description Helper function to conduct lasso or elastic net regualization
#'   for fitting a GLM in order to identify a multi-gene classifier
#'
#' @param traininglist A list object containing the training data. The output
#'   of the function `create_training_test_split()`.
#' @param type A character string. The algorihm used to perform calculations.
#'   Currently implemented are \emph{"grid", "lasso", "elastic"}.
#'
#' @param alpha A numeric between 0 and 1. The elastic net mixing parameter
#'   passed to `glmnet::glmnet()`.
#'
#' @inheritParams sigidentDiagnostic
#'
#' @export
signature <- function(traininglist,
                      type,
                      alpha = NULL,
                      nfolds = 10,
                      seed) {

  stopifnot(
    type %in% c("grid", "lasso", "elastic"),
    is.numeric(nfolds),
    is.numeric(seed),
    is.list(traininglist)
  )

  if (type == "grid") {
    outlist <- glmnet_gridsearch(traininglist, seed)
  } else {
    # use provided alpha only in elastic
    if (type == "lasso") {
      alpha <- 1
    } else if (type == "elastic") {
      stopifnot(!is.null(alpha),
                alpha >= 0 | alpha <= 1)
    }

    # initialize outlist
    outlist <- list()


    # go parallel
    ncores <- parallel::detectCores()
    cores <- ifelse(ncores > 4, 4, ncores)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)

    set.seed(seed)
    outlist$fit_cv <- glmnet::cv.glmnet(
      traininglist$train$x,
      traininglist$train$y,
      family = "binomial",
      type.measure = "mse",
      nfolds = nfolds,
      alpha = alpha,
      parallel = TRUE
    )
    # stop parallel computation
    parallel::stopCluster(cl)
    gc()

    # build the predictive models utilizing calculated lambda values
    glmpred <- build_predictive_glm(
      traininglist$train$x,
      traininglist$train$y,
      alpha = alpha,
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

    outlist$predicted_min <- pred_min$predicted
    outlist$confmat_min <- pred_min$confmat
    outlist$roc_min <- pred_min$roc


    # 1se
    pred_1se <- glm_prediction(
      model = outlist$lambda_1se,
      test_x = traininglist$test$x,
      test_y = traininglist$test$y,
      s = outlist$fit_cv$lambda.1se
    )

    outlist$predicted_1se <- pred_1se$predicted
    outlist$confmat_1se <- pred_1se$confmat
    outlist$roc_1se <- pred_1se$roc
  }

  return(outlist)
}


init_grid_search <- function() {
  # set up alpha and lambda grid to search for pair that minimizes CV erros
  lambda_grid <- 10 ^ seq(0, -4, length = 100)
  alpha_grid <- seq(0.1, 1, length = 10)

  # set up cross validation method for train function
  trn_ctrl <- caret::trainControl(method = "repeatedcv",
                                  number = 10)

  srch_grd <- expand.grid(.alpha = alpha_grid,
                          .lambda = lambda_grid)

  # set up search grid for alpha and lambda parameters
  return(list(srchGrd = srch_grd, trnCtrl = trn_ctrl))
}


glmnet_gridsearch <- function(traininglist, seed) {
  # initialize outlist
  outlist <- list()

  # initialize gridserach parameters
  gr_init <- init_grid_search()

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

#' @title gene_map_sig
#'
#' @description Helper function to map relevant input variables of a diagnostic
#'   model to corresponding IDs.
#'
#' @inheritParams plot_grid_model_plot
#' @inheritParams sigidentDEG
#'
#' @export
gene_map_sig <- function(mergeset, model) {
  id <- rownames(mergeset)
  # TODO warum i+1?
  index <- model[["beta"]]@i + 1
  # TODO map entrez_id on gene symbol here and include as second
  # column to ouput
  return(as.data.frame(x = cbind("ID" = id[index])))
}

#' @title validate_diagnostic_signature
#'
#' @description Helper function to validate diagnostic signatures
#'
#' @param validationstudylist A list, containing metainformation on the
#'   study used for validation of the diagnostic signature.
#' @param models A list of prediction models. Usually the output of the
#'   function `sigidentDiagnostic`.
#' @param datadir A character string. Path to the data-folder inside the
#'   metadata folder.
#'
#' @inheritParams sigidentDEG
#' @inheritParams plot_deg_heatmap
#' @inheritParams batch_correction
#' @inheritParams create_diagnosisdesignbatch
#'
#' @export
validate_diagnostic_signature <- function(validationstudylist,
                                          models,
                                          genes,
                                          idtype,
                                          targetname,
                                          controlname,
                                          targetcol,
                                          datadir) {
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

  eset <- load_eset(
    name = validationstudylist$studyname,
    datadir = datadir,
    targetcolname = validationstudylist$targetcolname,
    targetlevelname = validationstudylist$targetlevelname,
    controllevelname = validationstudylist$controllevelname,
    targetcol = targetcol,
    targetname = targetname,
    controlname = controlname,
    setid = validationstudylist$setid
  )



  diagnosis <- create_diagnosis(vector = eset[[targetcol]],
                                targetname = targetname,
                                controlname = controlname)


  expr <- create_expressionset(eset = eset,
                               idtype = idtype)

  # TODO why is this code in the original script?
  # # creating data frame, selecting only DEGs
  #% v.data.DEG <- base::subset(t(expr), select = genes)

  # creating data frame including all genes
  v_data_all <- t(expr)

  for (i in names(models)) {
    if (i %in% c("lasso", "elasticnet")) {
      for (j in c("min", "1se")) {
        predicted <- predict_glm(model = models[[i]][[j]]$model,
                                 test_x = v_data_all,
                                 type = "response")

        confmat <-
          caret::confusionMatrix(
            data = factor(ifelse(
              as.numeric(as.character(predicted)) < 0.5, 0, 1
            )),
            reference = factor(diagnosis),
            positive = "1"
          ) # determine the true case with the 'positive' argument

        # Calculate Roc
        roc <- calc_roc(test_y = diagnosis,
                        prediction = predicted)

        outlist[[i]][[j]] <-
          list(predicted = predicted,
               confmat = confmat,
               roc = roc)
      }
    }
  }
  return(outlist)
}
