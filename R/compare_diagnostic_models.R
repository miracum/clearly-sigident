#' @title compare_diagnostic_models
#'
#' @description Print an overview of all diagnostic models included
#'   in `modellist`.
#'
#' @param modellist A list object containing a list of diagnostic models.
#'   The output of the function `sigidentDiagnostic()`.
#'
#' @export
compare_diagnostic_models <- function(modellist) {
  outdat <- data.table::data.table("Model" = character(0),
                                   "AUC" = numeric(0))

  for (n in names(modellist)) {
    if (n %in% c("lasso", "elasticnet")) {
      for (m in names(modellist[[n]])) {
        if (m %in% c("min", "1se")) {
          outdat <-
            rbind(outdat,
                  data.table::data.table(
                    "Model" = paste(n, m, sep = "."),
                    "AUC" = round(modellist[[n]][[m]]$auc, 3)
                  ))
        }
      }
    } else if (n %in% c("glmnet", "svm", "knn", "rf", "gbm")) {
      outdat <-
        rbind(outdat,
              data.table::data.table(
                "Model" = n,
                "AUC" = round(modellist[[n]]$auc, 3)
              ))
    }
  }
  return(outdat)
}


#' @title get_diagnostic_lambda_values
#'
#' @description Print an overview of all lambda values of the diagnostic models
#'   in `modellist`.
#'
#' @inheritParams compare_diagnostic_models
#'
#' @export
get_diagnostic_lambda_values <- function(modellist) {
  outdat <-
    data.table::data.table(
      "Model" = character(0),
      "Lambda min" = numeric(0),
      "Lambda 1se" = numeric(0)
    )

  for (n in names(modellist)) {
    if (n %in% c("lasso", "elasticnet")) {
      outdat <- rbind(
        outdat,
        data.table::data.table(
          "Model" = n,
          "Lambda min" = round(modellist[[n]]$CV$lambda.min, 6),
          "Lambda 1se" = round(modellist[[n]]$CV$lambda.1se, 6)
        )
      )
    } else if (n %in% c("glmnet")) {
      outdat <-
        rbind(
          outdat,
          data.table::data.table(
            "Model" = paste(n, "(best tune: alpha, lambda)"),
            "Lambda min" = round(modellist[[n]]$model$bestTune$alpha, 6),
            "Lambda 1se" = round(modellist[[n]]$model$bestTune$lambda, 6)
          )
        )
    }
  }
  return(outdat)
}
