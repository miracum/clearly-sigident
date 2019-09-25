#' @title geneMapSig_
#'
#' @description Helper function to map relevant input variables of a diagnostic model to corresponding IDs.
#'
#' @inheritParams createGridModelPlot_
#' @inheritParams sigidentDEG
#'
#' @export
geneMapSig_ <- function(mergeset, model){
  id <- rownames(mergeset)
  # TODO warum i+1?
  index <- model[["beta"]]@i+1
  # TODO map entrez_id on gene symbol here and include as second columen to ouput
  return(as.data.frame(x = cbind("ID" = id[index])))
}



#' @title getSurvivalTime_
#'
#' @description Helper function to get survival time.
#'
#' @param discoverystudies.w.timedata A list that contains specifications on the study/studies that contain(s) survival time information.
#' @param datadir A character string. Path to the data-folder inside the metadata folder.
#'
#' @inheritParams sigidentDEG
#' @inheritParams createDEGheatmap_
#' @inheritParams batchCorrection_
#' @inheritParams createDiagnosisDesignBatch_
#'
#' @export
getSurvivalTime_ <- function(studyMetadata,
                             sampleMetadata,
                             genes,
                             idtype,
                             discoverystudies.w.timedata,
                             targetname,
                             controlname,
                             targetcol,
                             datadir){

  discovery <- discovery_(sampleMetadata = sampleMetadata,
                          studyMetadata = studyMetadata)
  stopifnot(
    is.character(targetcol),
    is.character(targetname),
    is.list(discoverystudies.w.timedata),
    length(setdiff(names(discoverystudies.w.timedata), discovery)) == 0
  )

  outlist <- list()

  for (st in names(discoverystudies.w.timedata)){
    samples <- sampleMetadata[sampleMetadata$study==st,]
    tumorsamples <- samples[eval(parse(text=paste0("samples$",targetcol)))==targetname,]


    eset <- loadEset_(name = st,
                      datadir = datadir,
                      targetcolname = discoverystudies.w.timedata[[st]]$targetcolname,
                      targetlevelname = discoverystudies.w.timedata[[st]]$targetlevelname,
                      controllevelname = discoverystudies.w.timedata[[st]]$controllevelname,
                      targetcol = targetcol,
                      targetname = targetname,
                      controlname = controlname)

    # filter only data of tumor samples
    if (!is.null(discoverystudies.w.timedata[[st]]$targetlevelname)){
      esetTargets <- eset[,which(eset$geo_accession %in% tumorsamples$sample)]
    } else {
      # use the whole dataset
      esetTargets <- eset
    }

    timestatus <- extractTimeStatus_(esetTargets = esetTargets,
                                     timecol = discoverystudies.w.timedata[[st]]$timecol,
                                     statuscol = discoverystudies.w.timedata[[st]]$status$statuscol,
                                     statuslevels = discoverystudies.w.timedata[[st]]$status$levels)
    time <- timestatus$time
    status <- timestatus$status

    survTable <- data.frame(time, status) # generating a table: columns: survivaltime, status; rows: individual samples

    expr <- createExpressionSet_(eset = esetTargets,
                                 idtype = idtype)

    ids <- rownames(expr)

    for (DEG in genes){
      exp <- exprsVector_(expr[DEG,]) # loop to bind a new column filled with exprs-values for each DEG
      DF <- data.frame(exp)
      colnames(DF) <- DEG
      survTable <- cbind(survTable, DF) # DF without batch effect removal, alternative would be rmBatch
    }

    rownames(survTable) = colnames(expr)

    outlist[[st]] <- list(survtable = survTable,
                          ids = ids)

  }
  return(outlist)
}


loadEset_ <- function(name, datadir, targetcolname, targetcol, targetname, controlname, targetlevelname, controllevelname){
  # original GEO data
  eset <- GEOquery::getGEO(name, destdir = datadir)[[1]]
  # rename targetcol
  colnames(Biobase::pData(eset))[which(colnames(Biobase::pData(eset)) == targetcolname)] <- targetcol

  # rename levels of targetcol
  if (!is.null(targetlevelname)){
    levelnames <- c(targetname, controlname)
    names(levelnames) <- c(targetlevelname,
                           controllevelname)
    eset[[targetcol]] <- plyr::revalue(eset[[targetcol]], levelnames)
  }
  return(eset)
}


#' @title univCox_
#'
#' @description Helper function to compute univariate cox regression and determine significance of each gene through separate univariate Cox regressions
#'
#' @param survtable A data.frame. Output of the function `getSurvivalTime_()`.
#' @param ids A character string. Output of the function `getSurvivalTime_()`.
#'
#' @export
univCox_ <- function(survtable, ids){

  covariates <- colnames(survtable)[-c(1:2)]
  covariates <- gsub("[[:punct:]]","",covariates)
  covariates <- paste("ID",covariates)
  covariates <- gsub(" ","",covariates)

  colnames(survtable)[-c(1:2)] <- covariates

  univ_formulas <- sapply(covariates, function(x) stats::as.formula(paste('survival::Surv(time, status) ~', x))) # build separate formula for each variable
  univ_models <- lapply(univ_formulas, function(f){survival::coxph(formula = f, data = survtable)}) # build cox model for each variable separatly

  # # parallel:
  # univ_models <- mcmapply(univ_formulas, FUN =  function(f) {
  #   survival::coxph(formula = f, data = survtable)},
  #   mc.cores = parallel::detectCores()-1, mc.preschedule = TRUE)

  # extract results
  univ_results <- lapply(univ_models, #
                         function(x){ #
                           x <- summary(x) #
                           p.value <- signif(x$wald["pvalue"], digits=2) #
                           wald.test <- signif(x$wald["test"], digits=2) #
                           beta <- signif(x$coef[1], digits=2);#coeficient beta # this block applies a function to every DEG,
                           HR <- signif(x$coef[2], digits=2);#exp(beta) # executing a cox regression analysis
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2) # DEG ~ survivaltime
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2) #
                           HR <- paste0(HR, " (", #
                                        HR.confint.lower, "-", HR.confint.upper, ")") #
                           res<-c(beta, HR, wald.test, p.value) #
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", #
                                         "p.value") #
                           return(res) #
                           #return(exp(cbind(coef(x),confint(x)))) #
                         }) #

  res <- t(as.data.frame(univ_results, check.names = FALSE)) #formatting of result
  result <- as.data.frame(res)
  result$p.value <- as.character(result$p.value) # reforamtting of p.value-column
  result$p.value <- as.numeric(result$p.value) #
  indices <- which(result$p.value <= 0.05) # selecting indices with a significant p-value
  IDs <- ids[indices] # IDs now contains Entrez-IDs of all genes with significant p-values

  # TODO symbols?

  sigCov <- result[with(result, which(result$p.value <= 0.05)),]
  rownames(sigCov) <- IDs
  sigCov <- cbind(id = IDs, sigCov)

  return(sigCov)
}


exprsVector_ <- function(DEG){
  vector <- c() # creates empty vector
  for(i in 1:length(DEG)){
    vector <- append(vector, DEG[[i]]) # vector is filled by columns of exprs (in sample order)
  }
  return(vector) # vector now contains one row (of the called gene) of exprs in sample/column order
}


#' @title generateExpressionPattern_
#'
#' @description Helper function to generate expression pattern
#'
#' @param classifier_studies A character vector. Names of the studies used to train the classifier.
#' @param sigCov A data.frame. Output of the function `univCox_()`.
#'
#' @inheritParams sigidentDEG
#' @inheritParams createDiagnosisDesignBatch_
#'
#' @export
generateExpressionPattern_ <- function(classifier_studies,
                                       sigCov,
                                       mergeset,
                                       studyMetadata,
                                       sampleMetadata,
                                       controlname,
                                       targetname,
                                       targetcol){

  stopifnot(
    length(setdiff(classifier_studies, studyMetadata$study)) == 0
  )

  dd <- createDiagnosisDesignBatch_(sampleMetadata = sampleMetadata[sampleMetadata$study %in% classifier_studies,],
                                    studyMetadata = studyMetadata[studyMetadata$study %in% classifier_studies,],
                                    controlname = controlname,
                                    targetname = targetname,
                                    targetcol = targetcol)
  diagnosis <- dd$diagnosis
  design <- dd$design

  # get sample names (= colnames of mergeset)
  classifier_colnames <- sampleMetadata[sampleMetadata$study %in% classifier_studies, get("sample")]

  # reduce mergedata to use only the data specified in classifier_studies
  classifier_data <- mergeset[,classifier_colnames]

  # generate vectors 'control' and 'tumor' for function expressionPattern
  control <- ctrl.vec_(classifier_data, diagnosis)
  tumor <- tum.vec_(classifier_data, diagnosis)

  IDs <- sigCov[,"id"]

  # train the prognostic classifier
  pattern <- c()
  for(i in 1:length(IDs)){
    pattern <- append(pattern, expressionPattern_(mergeset = classifier_data,
                                                  ids = IDs[i],
                                                  tumor = tumor,
                                                  control = control))
  }
  PatternCom <- data.frame(IDs, pattern)

  colnames(PatternCom) = c("Gene","Over/Under")

  return(PatternCom)
}




ctrl.vec_ <- function(mergeset, diagnosis){
  t1 <- colnames(mergeset)
  t1 <- cbind(t1, diagnosis)
  t1 <- as.data.frame(t1)
  ctrlM <- t1[with(t1, which(diagnosis == 0)),]
  control <- as.vector(ctrlM$t)
  return(control)
}

tum.vec_ <- function(mergeset, diagnosis){
  t1 <- colnames(mergeset)
  t1 <- cbind(t1,diagnosis)
  t1 <- as.data.frame(t1)
  tumM <- t1[with(t1, which(diagnosis == 1)),]
  tumor <- as.vector(tumM$t)
  return(tumor)
}



expressionPattern_ <- function(mergeset, ids, tumor, control){
  ctrl <- (mean(mergeset[ids, control])) # selecting control-samples and computing mean
  tum <-(mean(mergeset[ids, tumor])) # selecting tumor-samples and computing mean
  pattern=c() # creating empty vecotr
  if(ctrl < tum){
    pattern = "Over"
  }
  else if (ctrl > tum){ # Over-/Underexpression is checked by comparing group-means
    pattern = "Under" # vector is filled with respective symbols
  }
  else{pattern=NA}
  return(pattern) # pattern vector is returned
}



#' @title prognosticClassifier_
#'
#' @description Helper function to get survival time.
#'
#' @param validationstudiesinfo A list that contains specifications on the study that contains the validation information.
#' @param PatternCom A data.frame. The output of the function `generateExpressionPattern_()`.
#'
#' @inheritParams sigidentDEG
#' @inheritParams getSurvivalTime_
#' @inheritParams batchCorrection_
#' @inheritParams createDiagnosisDesignBatch_
#'
#' @export
prognosticClassifier_ <- function(PatternCom, idtype, validationstudiesinfo, datadir, targetcol, targetname, controlname){

  outlist <- list()

  for (st in names(validationstudiesinfo)){
    eset <- loadEset_(name = st,
                      datadir = datadir,
                      targetcolname = validationstudiesinfo[[st]]$targetcolname,
                      targetlevelname = validationstudiesinfo[[st]]$targetlevelname,
                      controllevelname = validationstudiesinfo[[st]]$controllevelname,
                      targetcol = targetcol,
                      targetname = targetname,
                      controlname = controlname)

    diagnosis <- diagnosis_(vector = eset[[targetcol]],
                            targetname = targetname,
                            controlname = controlname)

    # generate vectors 'control' and 'tumor' for function expressionPattern
    classifier_data <- Biobase::exprs(eset)
    control <- ctrl.vec_(classifier_data, diagnosis)
    tumor <- tum.vec_(classifier_data, diagnosis)

    # filter only data of tumor samples
    if (!is.null(validationstudiesinfo[[st]]$targetlevelname)){
      esetTargets <- eset[,which(eset[[targetcol]] == targetname)]
    } else {
      # use the whole dataset
      esetTargets <- eset
    }

    timestatus <- extractTimeStatus_(esetTargets = esetTargets,
                                     timecol = validationstudiesinfo[[st]]$timecol,
                                     statuscol = validationstudiesinfo[[st]]$status$statuscol,
                                     statuslevels = validationstudiesinfo[[st]]$status$levels)


    time <- timestatus$time
    status <- timestatus$status
    RiskTable <- data.frame(time, status)

    expr <- createExpressionSet_(eset = esetTargets,
                                 idtype = idtype)

    # classification and Kaplan-Meier estimator
    Sig <- sigAnalysis_(expr = expr, PatternCom = PatternCom)
    RiskTable.Sig <- classify_(Sig)

    # prepare classification and required data for Kaplan-Meier estimator
    Groups <- RiskTable.Sig[, "RiskGroup"]
    RiskTable <- cbind(RiskTable, Groups)
    rownames(RiskTable) = colnames(expr)

    kap <- fitKaplanEstimator_(RiskTable = RiskTable)

    # prepare classification and required data for Kaplan-Meier estimator
    outlist[[st]] <- list(kaplan.estimator = kap,
                          risktable = RiskTable)

  }
  return(outlist)
}

fitKaplanEstimator_ <- function(RiskTable){

  # fit proportional hazards regression model
  res.cox <- survival::coxph(survival::Surv(time, status) ~ Groups, data = RiskTable)
  new_df <- with(RiskTable,
                 data.frame(Groups = c(0, 1),
                            survival_time = rep(mean(time, na.rm = TRUE), 2),
                            ph.ecog = c(1, 1)
                 )
  )
  fit <- survival::survfit(res.cox, newdata = new_df)
  return(list(res.cox = res.cox,
              fit = fit))
}


classify_ <- function(Sigtable){
  table <- Sigtable
  Sums <- rowSums(table)
  vector <- c()
  for(i in 1:length(Sums)){
    if(Sums[i] > 0.5*length(colnames(table))){
      a <- 1
    }
    else{
      a <- 0
    }
    vector <- append(vector,a)
  }
  newColumn <- data.frame(vector)
  colnames(newColumn)[1] = "RiskGroup"
  table <- cbind(table,newColumn)
  return(table)
}


extractTimeStatus_ <- function(esetTargets, timecol, statuscol, statuslevels){

  # time <- esetTargets[, get("timecol")] -> does not work with expression sets
  time <- eval(parse(text=paste0("esetTargets$", timecol)))
  time <- as.character(time)
  time[which(time=="")] <- "time: not available"
  time <- unlist(lapply(time, strsplit, ": ")) # string is split
  time <- time[c(FALSE,TRUE)] # first part of string ("overall survival") is discarded
  time <- as.numeric(time) # results in a vector of survival times (in sample order)

  status <- eval(parse(text=paste0("esetTargets$", statuscol))) # coding if patient survived or not
  # 1 = no event/alive; 2 = event/death
  if (is.na(statuslevels$na)){
    levelnames <- c(1, 2)
    names(levelnames) <- c(statuslevels$alive,
                           statuslevels$deceased)
  } else {
    levelnames <- c(1, 2, NA)
    names(levelnames) <- c(statuslevels$alive,
                           statuslevels$deceased,
                           statuslevels$na)
  }
  status <- plyr::revalue(status, levelnames)
  status <- as.numeric(status)
  return(list(time=time, status=status))
}



over <- function(avrg, level){
  if(avrg < level){
    a <- 1
  }
  else{
    a <- 0
  }
  return(a)
}

# defining functions to check for Over-/Underexpression
under <- function(avrg, level){ # if expression pattern is matching the signature, 1 is returned, otherwise 0
  if(avrg>level){
    a <- 1
  }
  else{
    a <- 0
  }
  return(a)
}

sigAnalysis_ <- function(expr, PatternCom){ # Input is an eset, consisting of tumor-samples only and the signature to be tested

  Sigframe <- data.frame(c(1:length(colnames(expr))))

  for(i in 1:length(PatternCom[,"Gene"])){
    id <- as.character(PatternCom[i,"Gene"])
    avrg <- mean(expr[which(rownames(expr)==id),])
    vector <- c()

    for(j in 1:length(colnames(expr))){
      level <- expr[id,j]

      if(PatternCom[i,"Over/Under"] == "Over"){
        a <- over(avrg, level)
      }
      else if(PatternCom[i,"Over/Under"] == "Under"){
        a <- under(avrg, level)
      }
      else{
        print("Wrong formatting of Signature")
      }
      vector <- append(vector,a)
    }

    newColumn <- data.frame(vector)
    colnames(newColumn)[1] = id
    Sigframe <- cbind(Sigframe, newColumn)
    rownames(Sigframe) = colnames(expr)
  }
  Sigframe <- Sigframe[-1]
  return(Sigframe)
}

createExpressionSet_ <- function(eset, idtype){
  expr <- Biobase::exprs(eset)
  expr <- idType_(expr = expr, eset = eset, idtype = idtype)
  return(expr)
}


idType_ <- function(expr, eset, idtype){
  stopifnot(
    idtype %in% c("entrez", "affy")
  )
  if (idtype == "entrez"){
    rownames(expr) <- as.character(eset@featureData@data$ENTREZ_GENE_ID)
    # remove empty characters and replicates in EntrezIDs
    expr <- expr[rownames(expr)!="",]
    expr <- expr[!duplicated(rownames(expr)),]
  } else if (idtype == "affy") {
    rownames(expr) <- as.character(eset@featureData@data$ID)
  }
  return(expr)
}
