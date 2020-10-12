#' @title plot_rocplot
#'
#' @description Helper function to create roc plots
#'
#' @param filename A character string. The filename.
#' @param roc An object containing the roc information.
#'
#' @export
plot_rocplot <- function(roc,
                         filename) {

  grDevices::png(
    filename = filename,
    res = 150,
    height = 1000,
    width = 1500
  )
  print({
    graphics::plot(roc)
    graphics::text(0.4,
                   0,
                   paste0("AUC: ", round(roc$auc, 4)))
  })
  grDevices::dev.off()
}


#' @title plot_cvplot
#'
#' @description Helper function to create cross-validation plots
#'
#' @param cv_obj A object. The result of the function `sigident_signature()`.
#'
#' @inheritParams plot_rocplot
#'
#' @export
plot_cvplot <- function(cv_obj,
                        filename) {

  grDevices::png(
    filename = filename,
    res = 150,
    height = 1000,
    width = 1500
  )
  print({
    graphics::plot(cv_obj)
  })
  grDevices::dev.off()
}


#' @title plot_grid_model_plot
#'
#' @description Helper function to create the grid-model plot
#'
#' @param filename A character string. The filename.
#' @param model An object containing the model information.
#'
#' @export
plot_grid_model_plot <- function(model,
                                 filename) {

  grDevices::png(
    filename = filename,
    res = 150,
    height = 1000,
    width = 1500
  )
  print({
    graphics::plot(model)
  })
  grDevices::dev.off()
}


#' @title plot_grid_varimp_plot
#'
#' @description Helper function to create the variable importance plot
#'
#' @inheritParams plot_grid_model_plot
#'
#' @importFrom gbm relative.influence
#'
#' @export
plot_grid_varimp_plot <- function(model,
                                  filename) {

  assign(
    x = "relative.influence",
    value = gbm::relative.influence,
    pos = -1
  )

  var_imp <- caret::varImp(model)

  grDevices::png(
    filename = filename,
    res = 150,
    height = 1000,
    width = 1500
  )
  print({
    graphics::plot(var_imp, top = 20)
  })
  grDevices::dev.off()
}


#' @title plot_survplot
#'
#' @description Helper function to create survival plot
#'
#' @param fit A cox proportional hazards model. The output of the
#'   function `fit_kaplan_estimator()` or `prognostic_classifier()`.
#' @param risktable A data.frame. The output of the function
#'   `prognostic_classifier()`.
#'
#' @inheritParams plot_grid_model_plot
#'
#' @export
plot_survplot <- function(fit,
                          risktable,
                          filename) {

  grDevices::png(
    filename = filename,
    res = 150,
    height = 1000,
    width = 1500
  )
  print({
    survminer::ggsurvplot(
      fit,
      risktable,
      conf.int = TRUE,
      legend.labs = c("Low-Risk", "High-Risk")
    )
  })
  grDevices::dev.off()
}
