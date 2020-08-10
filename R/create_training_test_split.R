#' @title create_training_test_split
#'
#' @description Helper function to split data into training and test set
#'
#' @inheritParams sigidentPrognostic
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
