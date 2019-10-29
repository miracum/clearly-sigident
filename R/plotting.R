#' @title plot_import_boxplot
#'
#' @description Helper function to create boxplot of mergeset
#'
#' @param filename A character string indicating the filename. If default
#'   (\code{NULL}) a plot named `import_boxplot.png` will be created
#'   inside the directory "./plots".
#'
#' @inheritParams sigidentDEG
#'
#' @export
plot_import_boxplot <- function(mergeset,
                                 filename = NULL) {
  if (is.null(filename)) {
    filename <- "./plots/import_boxplot.png"
    if (!dir.exists("./plots/")) {
      dir.create("./plots/")
    }
  }

  outplot <- graphics::boxplot(
    mergeset,
    main = "Merged microarray data",
    xlab = "Samples",
    ylab = "Expression"
  )

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = 10,
    width = 20
  )
}


#' @title plot_batchplot
#'
#' @description Helper function to create batchplot
#'
#' @param filename A character string indicating the filename. If default
#'   (\code{NULL}) a plot named `PCplot{time}.png` will be created
#'   inside the directory "./plots".
#' @param time A character string indicating if the plot is "before" or
#'   "after" batch correction. This information is integrated
#'   into the filename.
#' @param correction_obj An object. The output of the function
#'   `batch_detection()`.
#'
#' @export
plot_batchplot <- function(correction_obj,
                             filename = NULL,
                             time) {

  if (is.null(filename)) {
    filename <- paste0("./plots/PCplot_", time, ".png")
    if (!dir.exists("./plots/")) {
      dir.create("./plots/")
    }
  }

  # time == "before" or "after"
  outplot <- gPCA::PCplot(
    correction_obj,
    ug = "guided",
    type = "1v2",
    main = paste("gPCA",
                 time,
                 "batch correction")
  )

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = 15,
    width = 15
  )
}


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
plot_deg_heatmap <-
  function(mergeset, genes, patientcolors, filename = NULL) {
    if (is.null(filename)) {
      filename <- "./plots/DEG_heatmap.png"
      if (!dir.exists("./plots/")) {
        dir.create("./plots/")
      }
    }

    outplot <- gplots::heatmap.2(
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

    ggplot2::ggsave(
      filename = filename,
      plot = outplot,
      device = "png",
      height = 20,
      width = 30
    )
  }


#' @title plot_enrichted_barplot
#'
#' @description Helper function to create enrichted barplots
#'
#' @param filename A character string indicating the filename. If default
#'   (\code{NULL}) a plot named `Enriched_{type}.png` will be created
#'   inside the directory "./plots".
#' @param enrichmentobj An object resulting from the function
#'   `go_enrichment_analysis()`.
#' @param type A character string. One of eiter "GO" or "KEGG".
#' @param showCategory An integer. Indicating the number of maximum categories
#'   to show in barplot.
#'
#' @export
plot_enrichted_barplot <- function(enrichmentobj,
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

  outplot <- graphics::barplot(
    enrichmentobj,
    showCategory = show_category) +
    ggplot2::ggtitle(paste0("Enriched ", type, " terms")) +
    ggplot2::ylab("Gene count")

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = 10,
    width = 20
  )
}


#' @title color_heatmap
#'
#' @description Helper function to color the heatmap
#'
#' @inheritParams sigidentDEG
#'
#' @export
color_heatmap <-
  function(sample_metadata,
           study_metadata,
           targetcol,
           controlname) {
    discovery <- discovery_func(
      sample_metadata = sample_metadata,
      study_metadata = study_metadata
    )

    discoverydata <- sample_metadata[which(
      sample_metadata$study %in% discovery
    ), ][[targetcol]]

    return(
      unlist(
        lapply(
          discoverydata, function(tumor) {
            #% healthy=blue
            return(
              ifelse(tumor == controlname, "#0000FF", "#FF0000")
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

  outplot <- graphics::plot(roc) +
    graphics::text(0.4,
                   0,
                   paste0("AUC: ", round(roc$auc, 4)))

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = 10,
    width = 15
  )
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

  outplot <- graphics::plot(cv_obj)

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = 10,
    width = 15
  )
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

  outplot <- graphics::plot(model)

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = 10,
    width = 15
  )
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
  outplot <- graphics::plot(var_imp, top = 20)

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = 10,
    width = 15
  )
}


#' @title plot_survplot
#'
#' @description Helper function to create survival plot
#'
#' @param fit A cox proportional hazards model. The output of the
#'   function `fit_kaplan_estimator()` or `prognostic_classifier()`.
#' @param RiskTable A data.frame. The output of the function
#'   `prognostic_classifier()`.
#'
#' @inheritParams plot_grid_model_plot
#'
#' @export
plot_survplot <- function(fit,
                          risk_table,
                          filename) {

  outplot <- survminer::ggsurvplot(
    fit,
    risk_table,
    conf.int = TRUE,
    legend.labs = c("Low-Risk", "High-Risk")
  )

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = 10,
    width = 15
  )
}
