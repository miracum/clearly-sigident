#' @title Perform Prognostic Signature Analyses in Gene Expression Datasets Derived from MicroArrays
#'
#' @description One function to perform prognostic signature analysis.
#'
#' @inheritParams getSurvivalTime_
#' @inheritParams sigidentDiagnostic
#' @inheritParams generateExpressionPattern_
#' @inheritParams prognosticClassifier_
#'
#' @import data.table
#' @importFrom magrittr "%>%"
#'
#'
#' @export
sigidentPrognostic <- function(mergeset,
                               studyMetadata,
                               sampleMetadata,
                               idtype,
                               discoverystudies.w.timedata,
                               classifier_studies,
                               validationstudiesinfo,
                               targetname,
                               controlname,
                               datadir,
                               plotdir = "./plots/",
                               csvdir = "./tables/",
                               targetcol = "target"){



  stopifnot(
    class(mergeset) == "matrix",
    is.character(plotdir),
    is.character(csvdir),
    is.character(controlname),
    is.character(targetname),
    is.character(targetcol),
    is.list(discoverystudies.w.timedata),
    is.list(validationstudiesinfo),
    is.character(classifier_studies),
    is.character(datadir),
    dir.exists(datadir),
    idtype %in% c("entrez", "affy")
  )


  # create internal list for storage
  rv <- list()

  # store names
  rv$controlname <- controlname
  rv$targetname <- targetname
  rv$targetcol <- targetcol

  # store dirs
  rv$plotdir <- cleanPathName_(plotdir)
  rv$csvdir <- cleanPathName_(csvdir)
  rv$datadir <- cleanPathName_(datadir)

  # create output directories
  dir.create(rv$plotdir)
  dir.create(rv$csvdir)

  # add mergedset to list
  rv$mergeset <- mergeset

  # initialize storage lists
  rv$survTable <- list()
  rv$entrezIDs <- list()
  rv$surv_correlated <- list()
  rv$exprPattern <- list()
  rv$pC <- list()
  rv$results <- list()

  # get survivalData
  rv$survivalData <- getSurvivalTime_(studyMetadata = studyMetadata,
                                      sampleMetadata = sampleMetadata,
                                      discoverystudies.w.timedata = discoverystudies.w.timedata,
                                      targetname = rv$targetname,
                                      controlname = rv$controlname,
                                      targetcol = rv$targetcol,
                                      datadir = rv$datadir)

  # outlist
  outlist <- list()
  for (i in names(discoverystudies.w.timedata)){
    # extract data
    rv$survTable[[i]] <- rv$survivalData[[i]]$survtable
    rv$ids[[i]] <- rv$survivalData[[i]]$ids

    # compute univariat cox regression
    rv$surv_correlated[[i]] <- univCox_(survtable = rv$survTable[[i]],
                                        ids = rv$ids[[i]])
    # export table with survival correlated genes
    data.table::fwrite(rv$surv_correlated[[i]], paste0(rv$csvdir, i, "_survival_correlated_genes.csv"))

    # evaluate expression pattern
    rv$exprPattern[[i]] <- generateExpressionPattern_(classifier_studies = classifier_studies,
                                                      sigCov = rv$surv_correlated[[i]],
                                                      mergeset = rv$mergeset,
                                                      studyMetadata = studyMetadata,
                                                      sampleMetadata = sampleMetadata,
                                                      controlname = rv$controlname,
                                                      targetname = rv$targetname,
                                                      targetcol = rv$targetcol)

    # apply prognostic classifier
    rv$pC[[i]] <- prognosticClassifier_(PatternCom = rv$exprPattern[[i]],
                                        validationstudiesinfo = validationstudiesinfo,
                                        datadir = rv$datadir,
                                        controlname = rv$controlname,
                                        targetname = rv$targetname,
                                        targetcol = rv$targetcol)

    for (n in names(validationstudiesinfo)){
      fit <- rv$pC[[i]][[n]]$kaplan.estimator$fit
      RiskTable <- rv$pC[[i]][[n]]$risktable

      rv$results[[i]][[n]] <- list(fit = fit, risktable = RiskTable)

      filename <- paste0(rv$plotdir, n, "_Prognostic_Kaplan-Meier_Plot.png")
      createSurvPlot_(fit = fit,
                      RiskTable = RiskTable,
                      filename = filename)
    }
  }
  return(rv$results)
}

