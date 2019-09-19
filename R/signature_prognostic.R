#' @title geneMapSig_
#'
#' @description Helper function to map relevant input variables of a diagnostic model to Entrez IDs.
#'
#' @inheritParams createGridModelPlot_
#' @inheritParams sigidentMicroarray
#'
#' @export
geneMapSig_ <- function(mergeset, model){
  entrez <- rownames(mergeset)
  # TODO warum i+1?
  index <- model[["beta"]]@i+1
  return(as.data.frame(x = cbind("Entrez_ID" = entrez[index])))
}

#' @title getSurvivalTime_
#'
#' @description Helper function to get survival time.
#'
#' @param discoverystudies.w.timedata A list that contains specifications on the study/studies that contain(s) survival time information.
#' @param timecol A character string. The name of the column containing the survival time information.
#' @param statuscol A character string. The name of the column containing the survival status information.
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
getSurvivalTime_ <- function(studyMetadata, sampleMetadata, discoverystudies.w.timedata, targetname, controlname, targetcol, datadir){

  discovery <- discovery_(sampleMetadata = sampleMetadata,
                          studyMetadata = studyMetadata)
  stopifnot(
    is.character(targetcol),
    is.character(targetname),
    is.list(discoverystudies.w.timedata),
    length(setdiff(names(discoverystudies.w.timedata), discovery) == 0)
  )

  for (st in names(discoverystudies.w.timedata)){
    whichsamples <- sampleMetadata[get("study")==st,]
    whichtumorsamples <- whichsamples[eval(parse(text=paste0("whichsamples$",targetcol)))==targetname,]

    # original GEO data
    eset <- GEOquery::getGEO(st, destdir = datadir)[[1]]
    colnames(Biobase::pData(eset))[which(colnames(Biobase::pData(eset)) == discoverystudies.w.timedata[[st]]$targetcolname)] <- targetcol
    levels(eset[[targetcol]]) = c(controlname, targetname)
    # filter only data of tumor samples
    esetTargets <- eset[which(eset$geo_accession %in% whichtumorsamples$sample),]
    time <- eval(parse(text=paste0("esetTargets$", discoverystudies.w.timedata[[st]]$timecol)))
    time <- as.character(time)
    time <- unlist(lapply(time, strsplit, ": ")) # string is split
    time <- time[c(FALSE,TRUE)] # first part of string ("overall survival") is discarded
    time <- as.numeric(time) # results in a vector of survival times (in sample order)

    status <- eval(parse(text=paste0("esetTargets$", discoverystudies.w.timedata[[st]]$statuscol))) # coding if patient survived or not
    levels(status) <- c(1,2,NA) # 1 = no event/alive; 2 = event/death
    status <- as.numeric(status)
    return(data.frame(time,status)) # generating a table: columns: survivaltime, status; rows: individual samples
  }
}

