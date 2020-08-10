calc_roc <- function(test_y,
                     prediction) {

  return(pROC::roc(
    response = as.numeric(as.character(test_y)),
    predictor = as.numeric(as.character(prediction))
  ))
}
