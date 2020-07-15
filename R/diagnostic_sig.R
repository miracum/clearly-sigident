#' @title signature
#'
#' @description Helper function to conduct lasso or elastic net regualization
#'   for fitting a GLM in order to identify a multi-gene classifier
#'
#' @param traininglist A list object containing the training data. The output
#'   of the function `create_training_test_split()`.
#' @param type A character string. The algorihm used to perform calculations.
#'   Currently implemented are \emph{"elasticnet_grid", "lasso", "elastic",
#'   "svm", "knn"}.
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
    type %in% c("elasticnet_grid",
                "lasso",
                "elastic",
                "svm",
                "knn",
                "rf"),
    is.numeric(nfolds),
    is.numeric(seed),
    is.list(traininglist)
  )

  if (type == "elasticnet_grid") {
    outlist <- glmnet_gridsearch(traininglist, seed, nfolds)
  } else if (type == "svm") {
    outlist <- svm_classifier(traininglist, seed, nfolds)
  } else if (type == "knn") {
    outlist <- knn_classifier(traininglist, seed, nfolds)
  } else if (type == "rf") {
    outlist <- rf_classifier(traininglist, seed, nfolds)
  } else {
    outlist <- glmnet_classifier(traininglist, type, seed, nfolds)
  }

  return(outlist)
}


#' @title gene_map_sig
#'
#' @description Helper function to map relevant input variables of a diagnostic
#'   model to corresponding IDs.
#'
#' @inheritParams plot_grid_model_plot
#' @inheritParams sigidentPrognostic
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

#' @title validate_diagnostic_signatures
#'
#' @description Helper function to validate diagnostic signatures
#'
#' @param validationstudylist A list, containing metainformation on the
#'   study used for validation of the diagnostic signature.
#' @param models A list of prediction models. Usually the output of the
#'   function `sigidentDiagnostic`.
#'
#' @inheritParams sigidentPrognostic
#'
#' @export
validate_diagnostic_signatures <- function(validationstudylist,
                                           models,
                                           genes,
                                           idtype,
                                           datadir) {
  stopifnot(
    is.list(validationstudylist)
  )

  targetcol <- "target"
  controlname <- "Control"
  targetname <- "Target"

  outlist <- list()

  for (st in names(validationstudylist)) {


    stopifnot(
      is.character(validationstudylist[[st]]$targetcolname),
      is.character(validationstudylist[[st]]$targetlevelname),
      is.character(validationstudylist[[st]]$controllevelname),
      is.numeric(validationstudylist[[st]]$setid)
    )

    # setd use_raw, if not provided with function arguments
    use_raw <- ifelse(
      is.null(validationstudylist[[st]]$use_rawdata),
      FALSE,
      validationstudylist[[st]]$use_rawdata
    )

    eset <- tryCatch(
      expr = {
        eset <- eval(parse(text = st),
                     envir = 1L)
        cat(paste0("\nLoaded ",
                   st,
                   " from .Globalenv...\n"))
        eset
      }, error = function(e) {
        eset <- sigident.preproc::geo_load_eset(
          name = st,
          datadir = datadir,
          targetcolname = validationstudylist[[st]]$targetcolname,
          targetcol = targetcol,
          targetname = targetname,
          controlname = controlname,
          targetlevelname = validationstudylist[[st]]$targetlevelname,
          controllevelname = validationstudylist[[st]]$controllevelname,
          use_rawdata = use_raw,
          setid = validationstudylist[[st]]$setid
        )
        cat(paste0("\nLoaded ",
                   st,
                   " from URL\n"))
        eset
      }, finally = function(f) {
        return(eset)
      }
    )

    diagnosis <- sigident.preproc::geo_create_diagnosis(
      vector = eset[[targetcol]],
      targetname = targetname,
      controlname = controlname
    )

    expr <- sigident.preproc::geo_create_expressionset(
      eset = eset,
      idtype = idtype
    )

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

          outlist[[st]][[i]][[j]] <- list(
            predicted = predicted,
            confmat = confmat,
            roc = roc
          )
        }
      } else if (i %in% "elasticnet_grid") {
        # TODO add other caret-baased models here
        predicted <- predict_glm(model = models[[i]]$model,
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

        outlist[[st]][[i]] <- list(
          predicted = predicted,
          confmat = confmat,
          roc = roc
        )
      } else if (i %in% c("svm", "knn", "rf")) {

        predicted <- stats::predict(
          models[[i]]$model,
          v_data_all
        )

        confmat <-
          caret::confusionMatrix(
            predicted,
            as.factor(diagnosis),
            positive = "1"
          )
        # determine the true case with the 'positive' argument

        # Calculate Roc
        roc <- calc_roc(test_y = as.factor(diagnosis),
                        prediction = predicted)

        outlist[[st]][[i]] <- list(
          predicted = predicted,
          confmat = confmat,
          roc = roc
        )
      }
    }
  }
  return(outlist)
}
