
#' @title predict_caret
#'
#' @description Function to classify given data based on the created
#'   model.
#'
#' @param model The caret model.
#' @param test_x The data to be classified.
#'
predict_caret <- function(model,
                          test_x) {

  outdat <- caret::predict.train(
    model,
    test_x,
    type = "prob"
  )[, "X1"]

  return(outdat)
}
