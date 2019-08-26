createImportHistogram_ <- function(mergeset, filename){
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


createBatchPlot_ <- function(correction_obj, filename, time){
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


createDEGheatmap_ <- function(combat, genes, patientcolors, filename){
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
