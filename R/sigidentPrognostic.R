#' @title Perform Prognostic Signature Analyses in Gene Expression
#'   Datasets Derived from MicroArrays
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
                               sample_metadata,
                               idtype,
                               genes,
                               discoverystudies_w_timedata,
                               classifier_studies,
                               validationstudiesinfo,
                               datadir,
                               plotdir = "./plots/",
                               csvdir = "./tables/") {
  stopifnot(
    class(mergeset) == "matrix",
    is.character(plotdir),
    is.character(csvdir),
    is.list(discoverystudies_w_timedata),
    is.list(validationstudiesinfo),
    is.character(classifier_studies),
    is.character(datadir),
    dir.exists(datadir),
    idtype %in% c("entrez", "affy")
  )

  targetcol <- "target"
  controlname <- "Control"
  targetname <- "Target"

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
  rv$entrez_ids <- list()
  rv$surv_correlated <- list()
  rv$expr_pattern <- list()
  rv$p_c <- list()
  rv$results <- list()

  # get survival_data
  rv$survival_data <-
    get_survival_time(
      sample_metadata = sample_metadata,
      discoverystudies_w_timedata = discoverystudies_w_timedata,
      idtype = rv$idtype,
      genes = rv$genes,
      datadir = rv$datadir
    )

  for (i in names(discoverystudies_w_timedata)) {
    # extract data
    rv$survival_table[[i]] <- rv$survival_data[[i]]

    # compute univariat cox regression
    rv$surv_correlated[[i]] <-
      univ_cox(survtable = rv$survival_table[[i]],
               genes = rv$genes)
    # export table with survival correlated genes
    data.table::fwrite(rv$surv_correlated[[i]],
                       paste0(rv$csvdir, i, "_survival_correlated_genes.csv"))

    # evaluate expression pattern
    rv$expr_pattern[[i]] <-
      generate_expression_pattern(
        classifier_studies = classifier_studies,
        sig_cov = rv$surv_correlated[[i]],
        mergeset = rv$mergeset,
        sample_metadata = sample_metadata
      )

    # apply prognostic classifier
    rv$p_c[[i]] <-
      prognostic_classifier(
        pattern_com = rv$expr_pattern[[i]],
        validationstudiesinfo = validationstudiesinfo,
        idtype = rv$idtype,
        datadir = rv$datadir
      )

    for (n in names(validationstudiesinfo)) {
      fit <- rv$p_c[[i]][[n]]$kaplan_estimator$fit
      risk_table <- rv$p_c[[i]][[n]]$risktable

      rv$results[[i]][[n]] <- list(fit = fit,
                                   risktable = risk_table)

      filename <- paste0(rv$plotdir,
                         n,
                         "_Prognostic_Kaplan-Meier_Plot.png")
      plot_survplot(fit = fit,
                    risk_table = risk_table,
                    filename = filename)
    }
  }
  return(rv$results)
}
