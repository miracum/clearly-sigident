#' @title Perform Prognostic Signature Analyses in Gene Expression Datasets Derived from MicroArrays
#'
#' @description One function to perform prognostic signature analysis.
#'
#' @inheritParams get_survival_time
#' @inheritParams sigidentDEG
#' @inheritParams generate_expression_pattern
#' @inheritParams prognostic_classifier
#'
#' @import data.table
#' @importFrom magrittr "%>%"
#'
#'
#' @export
sigidentPrognostic <- function(mergeset,
                               study_metadata,
                               sample_metadata,
                               idtype,
                               genes,
                               discoverystudies_w_timedata,
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
    is.list(discoverystudies_w_timedata),
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

  # store other variables
  rv$idtype <- idtype
  rv$genes <- genes

  # store dirs
  rv$plotdir <- clean_path_name(plotdir)
  rv$csvdir <- clean_path_name(csvdir)
  rv$datadir <- clean_path_name(datadir)

  # create output directories
  dir.create(rv$plotdir)
  dir.create(rv$csvdir)

  # add mergedset to list
  rv$mergeset <- mergeset

  # initialize storage lists
  rv$survival_table <- list()
  rv$entrezIDs <- list()
  rv$surv_correlated <- list()
  rv$exprPattern <- list()
  rv$pC <- list()
  rv$results <- list()

  # get survivalData
  rv$survivalData <- get_survival_time(study_metadata = study_metadata,
                                      sample_metadata = sample_metadata,
                                      discoverystudies_w_timedata = discoverystudies_w_timedata,
                                      idtype = rv$idtype,
                                      genes = rv$genes,
                                      targetname = rv$targetname,
                                      controlname = rv$controlname,
                                      targetcol = rv$targetcol,
                                      datadir = rv$datadir)

  # outlist
  outlist <- list()
  for (i in names(discoverystudies_w_timedata)){
    # extract data
    rv$survival_table[[i]] <- rv$survivalData[[i]]

    # compute univariat cox regression
    rv$surv_correlated[[i]] <- univ_cox(survtable = rv$survival_table[[i]],
                                        genes = rv$genes)
    # export table with survival correlated genes
    data.table::fwrite(rv$surv_correlated[[i]], paste0(rv$csvdir, i, "_survival_correlated_genes.csv"))

    # evaluate expression pattern
    rv$exprPattern[[i]] <- generate_expression_pattern(classifier_studies = classifier_studies,
                                                      sig_cov = rv$surv_correlated[[i]],
                                                      mergeset = rv$mergeset,
                                                      study_metadata = study_metadata,
                                                      sample_metadata = sample_metadata,
                                                      controlname = rv$controlname,
                                                      targetname = rv$targetname,
                                                      targetcol = rv$targetcol)

    # apply prognostic classifier
    rv$pC[[i]] <- prognostic_classifier(pattern_com = rv$exprPattern[[i]],
                                        validationstudiesinfo = validationstudiesinfo,
                                        idtype = rv$idtype,
                                        datadir = rv$datadir,
                                        controlname = rv$controlname,
                                        targetname = rv$targetname,
                                        targetcol = rv$targetcol)

    for (n in names(validationstudiesinfo)){
      fit <- rv$pC[[i]][[n]]$kaplan_estimator$fit
      RiskTable <- rv$pC[[i]][[n]]$risktable

      rv$results[[i]][[n]] <- list(fit = fit, risktable = RiskTable)

      filename <- paste0(rv$plotdir, n, "_Prognostic_Kaplan-Meier_Plot.png")
      plot_survplot(fit = fit,
                      risk_table =RiskTable,
                      filename = filename)
    }
  }
  return(rv$results)
}

