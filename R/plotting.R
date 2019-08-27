createImportHistogram_ <- function(mergeset, filename = NULL){
  if (is.null(filename)){
    filename <- "./plots/import_histogram.png"
    if (!dir.exists("./plots/")){
      dir.create("./plots/")
    }
  }
  shiny::plotPNG({
    # add plot and plot statistics here, "j" is necessary to get values for curve in equations
    return(print(graphics::boxplot(mergeset@assayData$exprs, main = "Merged data before batch correction",
                                   xlab = "Samples", ylab ="Expression value")
    ))
  },
  filename = filename,
  height = 400,
  width = 450)
}


createBatchPlot_ <- function(correction_obj, filename = NULL, time){
  if (is.null(filename)){
    filename <- paste0("./plots/PCplot_", time, ".png")
    if (!dir.exists("./plots/")){
      dir.create("./plots/")
    }
  }
  # time == "before" or "after"
  shiny::plotPNG({
    # add plot and plot statistics here, "j" is necessary to get values for curve in equations
    return(print(
      gPCA::PCplot(
        correction_obj, ug="guided", type="1v2", main = paste("gPCA", time, "batch correction")
      )
    ))
  },
  filename = filename,
  height = 400,
  width = 450)
}


createDEGheatmap_ <- function(combat, genes, patientcolors, filename = NULL){
  if (is.null(filename)){
    filename <- "./plots/DEG_heatmap.png"
    if (!dir.exists("./plots/")){
      dir.create("./plots/")
    }
  }
  # time == "before" or "after"
  shiny::plotPNG({
    # add plot and plot statistics here, "j" is necessary to get values for curve in equations
    return(print(
      gplots::heatmap.2(combat[genes,],
                        ColSideColors= patientcolors,
                        key= TRUE,
                        symkey= FALSE,
                        density.info= "none",
                        scale = "none",
                        trace= "none",
                        col = topo.colors(100),
                        cexRow = 0.4,
                        cexCol = 0.4)
    ))
  },
  filename = filename,
  height = 400,
  width = 450)
}

colorHeatmap_ <- function(tumor, controlname){
  # healthy lung=blue
  if(Tumor==controlname) "#0000FF" else "#FF0000"
}
