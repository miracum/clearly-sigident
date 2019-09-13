#' @title createImportBoxplot_
#'
#' @description Helper function to create boxplot of mergeset
#'
#' @param filename A character string indicating the filename. If default (\code{NULL}) a plot named `import_boxplot.png` will be created
#'   inside the directory "./plots".
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
createImportBoxplot_ <- function(mergeset, filename = NULL){
  if (is.null(filename)){
    filename <- "./plots/import_boxplot.png"
    if (!dir.exists("./plots/")){
      dir.create("./plots/")
    }
  }
  shiny::plotPNG({
    return(print(
      graphics::boxplot(mergeset, main = "Merged microarray data recieved by mergeGEO()",
                                   xlab = "Samples", ylab ="Expression")
    ))
  },
  filename = filename,
  height = 1000,
  width = 1500)
}


#' @title createBatchPlot_
#'
#' @description Helper function to create batchplot
#'
#' @param filename A character string indicating the filename. If default (\code{NULL}) a plot named `PCplot{time}.png` will be created
#'   inside the directory "./plots".
#' @param time A character string indicating if the plot is "before" or "after" batch correction. This information is integrated
#'   into the filename.
#' @param correction_obj An object. The output of the function `batchCorrection_()`.
#'
#' @export
createBatchPlot_ <- function(correction_obj, filename = NULL, time){
  if (is.null(filename)){
    filename <- paste0("./plots/PCplot_", time, ".png")
    if (!dir.exists("./plots/")){
      dir.create("./plots/")
    }
  }
  # time == "before" or "after"
  shiny::plotPNG({
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


#' @title createDEGheatmap_
#'
#' @description Helper function to create DEG heatmap
#'
#' @param filename A character string indicating the filename. If default (\code{NULL}) a plot named `DEG_heatmap.png` will be created
#'   inside the directory "./plots".
#' @param genes A object. The output of the function `identifyDEGs_()`.
#' @param patientcolors A object. The ouput of the function `colorHeatmap_()`.
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
createDEGheatmap_ <- function(mergeset, genes, patientcolors, filename = NULL){
  if (is.null(filename)){
    filename <- "./plots/DEG_heatmap.png"
    if (!dir.exists("./plots/")){
      dir.create("./plots/")
    }
  }
  shiny::plotPNG({
    return(print(
      gplots::heatmap.2(mergeset[genes,],
                        ColSideColors= patientcolors,
                        key= TRUE,
                        symkey= FALSE,
                        density.info= "none",
                        scale = "none",
                        trace= "none",
                        col = grDevices::topo.colors(100),
                        cexRow = 0.4,
                        cexCol = 0.4)
    ))
  },
  filename = filename,
  height = 1000,
  width = 1500)
}


#' @title createEnrichtedBarplot_
#'
#' @description Helper function to create enrichted barplots
#'
#' @param filename A character string indicating the filename. If default (\code{NULL}) a plot named `Enriched_{type}.png` will be created
#'   inside the directory "./plots".
#' @param enrichmentobj An object resulting from the function `goEnrichmentAnalysis_()`.
#' @param type A character string. One of eiter "GO" or "KEGG".
#' @param showCategory An integer. Indicating the number of maximum categories to show in barplot.
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
createEnrichtedBarplot_ <- function(enrichmentobj, type, filename = NULL, showCategory = 20){
  stopifnot(
    is.character(type),
    type %in% c("GO", "KEGG"),
    is.numeric(showCategory)
  )
  if (is.null(filename)){
    filename <- paste0("./plots/Enriched_", type, ".png")
    if (!dir.exists("./plots/")){
      dir.create("./plots/")
    }
  }
  shiny::plotPNG({
    return(print(graphics::barplot(enrichmentobj, showCategory = showCategory) +
                   ggplot2::ggtitle(paste0("Enriched ", type, " terms")) +
                   ggplot2::ylab("Gene count")
                 ))
  },
  filename = filename,
  height = 1000,
  width = 1500)
}


#' @title colorHeatmap_
#'
#' @description Helper function to color the heatmap
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
colorHeatmap_ <- function(sampleMetadata, studyMetadata, targetcol, controlname){

  discovery <- studyMetadata[which(studyMetadata$discovery), "study"]
  discoverydata <- sampleMetadata[which(sampleMetadata$study %in% discovery),][[targetcol]]

  return(
    unlist(lapply(discoverydata, function(tumor){
      # healthy=blue
      return(ifelse(tumor==controlname, "#0000FF", "#FF0000"))
    }))
  )
}


#' @title createROCplot_
#'
#' @description Helper function to create roc plots
#'
#' @param filename A character string. The filename.
#' @param roc An object containing the roc information. The ouput of the function `glmnetSignature_()`
#'   with \code{type = "elastic"} or \code{type = "lasso"}.
#'
#' @export
createROCplot_ <- function(roc, filename){
  shiny::plotPNG({
    return(print({graphics::plot(roc)
                 graphics::text(0.4, 0, paste0("AUC: ", round(roc$auc, 4)))
                 }
    ))
  },
  filename = filename,
  height = 1000,
  width = 1500)
}


#' @title createCVPlot_
#'
#' @description Helper function to create cross-validation plots
#'
#' @param cv_obj A object. The result of the function `signature_()`.
#'
#' @inheritParams createROCplot_
#'
#' @export
createCVPlot_ <- function(cv_obj, filename){
  shiny::plotPNG({
    return(print(
      graphics::plot(cv_obj)
    ))
  },
  filename = filename,
  height = 400,
  width = 450)
}


#' @title createGridModelPlot_
#'
#' @description Helper function to create the grid-model plot
#'
#' @param filename A character string. The filename.
#' @param model An object containing the model information. The ouput of the function `glmnetSignature_()` with \code{type = "grid"}.
#'
#' @export
createGridModelPlot_ <- function(model, filename){
  shiny::plotPNG({
    return(print({graphics::plot(model)}))
  },
  filename = filename,
  height = 1000,
  width = 1500)
}


#' @title createGridVarImpPlot_
#'
#' @description Helper function to create the variable importance plot
#'
#' @inheritParams createGridModelPlot_
#'
#' @export
createGridVarImpPlot_ <- function(model, filename){
  shiny::plotPNG({
    varImp <- caret::varImp(model)
    return(print({graphics::plot(varImp, top = 20)}))
  },
  filename = filename,
  height = 1000,
  width = 1500)
}
