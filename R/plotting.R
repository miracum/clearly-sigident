#' @title plot_deg_heatmap
#'
#' @description Helper function to create DEG heatmap
#'
#' @param filename A character string indicating the filename. If default
#'   (\code{NULL}) a plot named `DEG_heatmap.png` will be created
#'   inside the directory "./plots".
#' @param genes A object. The output of the function `identify_degs()`.
#' @param patientcolors A object. The ouput of the function `color_heatmap()`.
#'
#' @inheritParams sigidentDEG
#'
#' @export
plot_deg_heatmap <- function(mergeset,
                             genes,
                             patientcolors,
                             filename = NULL) {
  if (is.null(filename)) {
    filename <- "./plots/DEG_heatmap.png"
    if (!dir.exists("./plots/")) {
      dir.create("./plots/")
    }
  }

  grDevices::png(
    filename = filename,
    res = 150,
    height = 2000,
    width = 3000
  )
  print({
    gplots::heatmap.2(
      mergeset[genes, ],
      ColSideColors = patientcolors,
      key = TRUE,
      symkey = FALSE,
      density.info = "none",
      scale = "none",
      trace = "none",
      col = grDevices::topo.colors(100),
      cexRow = 0.4,
      cexCol = 0.4
    )
  })
  grDevices::dev.off()
}


#' @title plot_enriched_barplot
#'
#' @description Helper function to create enrichted barplots
#'
#' @param filename A character string indicating the filename. If default
#'   (\code{NULL}) a plot named `Enriched_{type}.png` will be created
#'   inside the directory "./plots".
#' @param enrichmentobj An object resulting from the function
#'   `go_enrichment_analysis()`.
#' @param type A character string. One of eiter "GO" or "KEGG".
#' @param show_category An integer. Indicating the number of maximum categories
#'   to show in barplot.
#'
#' @export
plot_enriched_barplot <- function(enrichmentobj,
                                  type,
                                  filename = NULL,
                                  show_category = 20) {

  stopifnot(is.character(type),
            type %in% c("GO", "KEGG"),
            is.numeric(show_category))

  if (is.null(filename)) {
    filename <- paste0("./plots/Enriched_", type, ".png")
    if (!dir.exists("./plots/")) {
      dir.create("./plots/")
    }
  }

  grDevices::png(
    filename = filename,
    res = 150,
    height = 1000,
    width = 2000
  )
  print({
    graphics::barplot(
      enrichmentobj,
      showCategory = show_category) +
      ggplot2::ggtitle(paste0("Enriched ", type, " terms")) +
      ggplot2::ylab("Gene count")
  })
  grDevices::dev.off()
}


#' @title color_heatmap
#'
#' @description Helper function to color the heatmap
#'
#' @inheritParams sigidentDEG
#'
#' @export
color_heatmap <-
  function(sample_metadata) {

    targetcol <- "target"
    controlname <- "Control"

    discoverydata <- sample_metadata[[targetcol]]

    return(
      unlist(
        lapply(
          discoverydata, function(tumor) {
            #% healthy=blue
            return(
              ifelse(tumor == controlname,
                     "#0000FF",
                     "#FF0000")
            )
          }
        )
      )
    )
  }


#' @title plot_rocplot
#'
#' @description Helper function to create roc plots
#'
#' @param filename A character string. The filename.
#' @param roc An object containing the roc information.
#'   The ouput of the function `glmnetSignature_()`
#'   with \code{type = "elastic"} or \code{type = "lasso"}.
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
#' @param cv_obj A object. The result of the function `signature()`.
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
#' @param model An object containing the model information. The
#'   ouput of the function `glmnetSignature_()` with \code{type = "grid"}.
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
#' @export
plot_grid_varimp_plot <- function(model,
                                  filename) {


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
#' @param risk_table A data.frame. The output of the function
#'   `prognostic_classifier()`.
#'
#' @inheritParams plot_grid_model_plot
#'
#' @export
plot_survplot <- function(fit,
                          risk_table,
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
      risk_table,
      conf.int = TRUE,
      legend.labs = c("Low-Risk", "High-Risk")
    )
  })
  grDevices::dev.off()
}
