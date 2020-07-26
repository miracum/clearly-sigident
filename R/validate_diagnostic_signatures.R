#' @title validate_diagnostic_signatures
#'
#' @description Helper function to validate diagnostic signatures
#'
#' @param validationstudylist A list, containing metainformation on the
#'   study used for validation of the diagnostic signature.
#' @param models A list of prediction models. Usually the output of the
#'   function `sigidentDiagnostic`.
#' @param colindices A vector of integers. Indices of columns to run the
#'   validation on (default = NULL).
#'
#' @inheritParams sigidentPrognostic
#'
#' @export
validate_diagnostic_signatures <- function(validationstudylist,
                                           models,
                                           genes,
                                           idtype,
                                           datadir,
                                           colindices = NULL) {
  stopifnot(
    is.list(validationstudylist),
    is.null(colindices) || is.numeric(colindices)
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

    # creating data frame including all genes
    v_data_all <- t(expr)

    if (!is.null(colindices)) {
      v_data_all <- v_data_all[, colindices]
    }

    for (i in names(models)) {
      if (i %in% c("lasso", "elasticnet")) {

        message(paste0("Validating signature using model: ", i))

        for (j in c("min", "1se")) {
          prediction <- predict_glm(model = models[[i]][[j]]$model,
                                    test_x = v_data_all,
                                    type = "response")

          confmat <-
            caret::confusionMatrix(
              data = factor(ifelse(
                as.numeric(as.character(prediction)) < 0.5, 0, 1
              )),
              reference = factor(diagnosis),
              positive = "1"
            ) # determine the true case with the 'positive' argument

          # Calculate Roc
          roc <- calc_roc(test_y = diagnosis,
                          prediction = prediction)

          outlist[[st]][[i]][[j]] <- list(
            predicted = prediction,
            confmat = confmat,
            roc = roc
          )
        }
      } else {
        if (i %in% "elasticnet_grid") {

          message(paste0("Validating signature using model: ", i))

          prediction <- predict_glm(model = models[[i]]$model,
                                    test_x = v_data_all,
                                    type = "response")

        } else if (i %in% c("svm", "knn", "rf")) {

          message(paste0("Validating signature using model: ", i))

          prediction <- predict_caret(
            model = models[[i]]$model,
            test_x = v_data_all
          )
        }

        # create confmat
        confmat <-
          caret::confusionMatrix(
            data = factor(ifelse(
              as.numeric(as.character(prediction)) < 0.5, 0, 1
            )),
            reference = factor(diagnosis),
            positive = "1"
          ) # determine the true case with the 'positive' argument

        # Calculate Roc
        roc <- calc_roc(test_y = diagnosis,
                        prediction = prediction)

        outlist[[st]][[i]] <- list(
          prediction = prediction,
          confmat = confmat,
          roc = roc
        )
      }
    }
  }
  return(outlist)
}
