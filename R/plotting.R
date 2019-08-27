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
  height = 1000,
  width = 1500)
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
  height = 1000,
  width = 1500)
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
  height = 1000,
  width = 1500)
}


createEnrichtedBarplot_ <- function(enrichmentobj, type, filename = NULL, showCategory = 20){
  stopifnot(
    is.character(type),
    type %in% c("GO", "KEGG")
  )
  if (is.null(filename)){
    filename <- paste0("./plots/Enriched_", type, ".png")
    if (!dir.exists("./plots/")){
      dir.create("./plots/")
    }
  }
  shiny::plotPNG({
    # add plot and plot statistics here, "j" is necessary to get values for curve in equations
    return(print(graphics::barplot(enrichmentobj, showCategory = showCategory) +
                   ggplot2::ggtitle(paste0("Enriched ", type, " terms")) +
                   ggplot2::ylab("Gene count")
                 ))
  },
  filename = filename,
  height = 1000,
  width = 1500)
}

colorHeatmap_ <- function(mergeset, targetcol, controlname){
  return(
    unlist(lapply(mergeset[[targetcol]], function(tumor){
      # healthy=blue
      return(ifelse(tumor==controlname, "#0000FF", "#FF0000"))
    }))
  )
}

createROCplot_ <- function(roc, filename){
  shiny::plotPNG({
    # add plot and plot statistics here, "j" is necessary to get values for curve in equations
    return(print({graphics::plot(roc)
                 graphics::text(0.4, 0, paste0("AUC: ", round(roc$auc, 4)))
                 }
    ))
  },
  filename = filename,
  height = 1000,
  width = 1500)
}


createCVPlot_ <- function(cv_obj, filename){
  # time == "before" or "after"
  shiny::plotPNG({
    # add plot and plot statistics here, "j" is necessary to get values for curve in equations
    return(print(
      graphics::plot(cv_obj)
    ))
  },
  filename = filename,
  height = 400,
  width = 450)
}


