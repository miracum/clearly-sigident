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
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
getSurvivalTime_ <- function(studyMetadata,
                             sampleMetadata,
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
    whichsamples <- sampleMetadata[get("study")==st,]
    whichtumorsamples <- whichsamples[eval(parse(text=paste0("whichsamples$",targetcol)))==targetname,]


    eset <- loadEset_(name = st,
                      datadir = datadir,
                      targetcolname = discoverystudies.w.timedata[[st]]$targetcolname,
                      targetlevelname = discoverystudies.w.timedata[[st]]$targetlevelname,
                      controllevelname = discoverystudies.w.timedata[[st]]$controllevelname,
                      targetcol = targetcol,
                      targetname = targetname,
                      controlname = controlname)

    # filter only data of tumor samples
    esetTargets <- eset[which(eset$geo_accession %in% whichtumorsamples$sample),]

    timestatus <- extractTimeStatus_(esetTargets = esetTargets,
                                     timecol = discoverystudies.w.timedata[[st]]$timecol,
                                     statuscol = discoverystudies.w.timedata[[st]]$statuscol)
    time <- timestatus$time
    status <- timestatus$status

    survTable <- data.frame(time,status) # generating a table: columns: survivaltime, status; rows: individual samples

    expr <- Biobase::exprs(esetTargets)
    rownames(expr) <- as.character(esetTargets@featureData@data$ENTREZ_GENE_ID)
    # remove empty characters and replicates in EntrezIDs
    expr <- expr[rownames(expr)!="",]
    expr <- expr[!duplicated(rownames(expr)),]

    entrezIDs <- rownames(expr)

    for (DEG in entrezIDs){
      exp <- exprsVector_(expr[DEG,]) # loop to bind a new column filled with exprs-values for each DEG
      DF <- data.frame(exp)
      colnames(DF) <- DEG
      survTable <- cbind(survTable, DF) # DF without batch effect removal, alternative would be rmBatch
    }

    rownames(survTable) = colnames(Biobase::exprs(esetTargets))

    outlist[[st]] <- list(survTable = survTable,
                          entrezIDs = entrezIDs)

  }
  return(outlist)
}


loadEset_ <- function(name, datadir, targetcolname, targetcol, targetname, controlname, targetlevelname, controllevelname){
  # original GEO data
  eset <- GEOquery::getGEO(name, destdir = datadir)[[1]]
  # rename targetcol
  colnames(Biobase::pData(eset))[which(colnames(Biobase::pData(eset)) == targetcolname)] <- targetcol

  # rename levels of targetcol
  levelnames <- c(targetname, controlname)
  names(levelnames) <- c(targetlevelname,
                         controllevelname)
  eset[[targetcol]] <- plyr::revalue(eset[[targetcol]], levelnames)
  return(eset)
}


#' @title univCox_
#'
#' @description Helper function to compute univariate cox regression and determine significance of each gene through separate univariate Cox regressions
#'
#' @param survTable A data.frame. Output of the function `getSurvivalTime_()`.
#' @param entrezIDs A character string. Output of the function `getSurvivalTime_()`.
#'
#' @export
univCox_ <- function(survTable, entrezIDs){

  covariates <- colnames(survTable)[-c(1:2)]
  covariates <- gsub("_|/","",covariates)
  covariates <- paste("ENTREZID",covariates)
  covariates <- gsub(" ","",covariates)

  colnames(survTable)[-c(1:2)] <- covariates

  univ_formulas <- sapply(covariates, function(x) as.formula(paste('survival::Surv(time, status) ~', x))) # build separate formula for each variable
  univ_models <- lapply(univ_formulas, function(f){survival::coxph(formula = f, data = survTable)}) # build cox model for each variable separatly
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
  IDs <- entrezIDs[indices] # IDs now contains Entrez-IDs of all genes with significant p-values

  sigCov <- result[with(result, which(result$p.value <= 0.05)),]
  rownames(sigCov) <- IDs
  sigCov <- cbind(Entrez_ID = IDs, sigCov)

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
#' @inheritParams sigidentMicroarray
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

  dd <- sigident::createDiagnosisDesign_(sampleMetadata = sampleMetadata[sampleMetadata$study %in% classifier_studies,],
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

  IDs <- sigCov[,"Entrez_ID"]

  # train the prognostic classifier
  pattern <- c()
  for(i in 1:length(IDs)){
    pattern <- append(pattern, expressionPattern_(mergeset = classifier_data,
                                                  entrezID = IDs[i],
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



expressionPattern_ <- function(mergeset, entrezID, tumor, control){
  ctrl <- (mean(mergeset[entrezID, control])) # selecting control-samples and computing mean
  tum <-(mean(mergeset[entrezID, tumor])) # selecting tumor-samples and computing mean
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
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
prognosticClassifier_ <- function(PatternCom, validationstudiesinfo, datatdir, targetcol, targetname, controlname){

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
    esetTargets <- eset[which(eset[[targetcol]] == targetname),]

    timestatus <- extractTimeStatus_(esetTargets = esetTargets,
                                     timecol = validationstudiesinfo[[st]]$timecol,
                                     statuscol = validationstudiesinfo[[st]]$statuscol)


    time <- timestatus$time
    status <- timestatus$status
    RiskTable <- data.frame(time, status)

    # classification and Kaplan-Meier estimator
    Sig <- sigAnalysis_(eset = esetTargets, PatternCom = PatternCom)
  }
}

extractTimeStatus_ <- function(esetTargets, timecol, statuscol){
  time <- eval(parse(text=paste0("esetTargets$", timecol)))
  time <- as.character(time)
  time[which(time=="")] <- "time: not available"
  time <- unlist(lapply(time, strsplit, ": ")) # string is split
  time <- time[c(FALSE,TRUE)] # first part of string ("overall survival") is discarded
  time <- as.numeric(time) # results in a vector of survival times (in sample order)

  status <- eval(parse(text=paste0("esetTargets$", statuscol))) # coding if patient survived or not
  levels(status) <- c(1,2,NA) # 1 = no event/alive; 2 = event/death
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

sigAnalysis_ <- function(eset, PatternCom){ # Input is an eset, consisting of tumor-samples only and the signature to be tested

  expr <- Biobase::exprs(eset)
  rownames(expr) <- as.character(esetTargets@featureData@data$ENTREZ_GENE_ID)
  # remove empty characters and replicates in EntrezIDs
  expr <- expr[rownames(expr)!="",]
  expr <- expr[!duplicated(rownames(expr)),]

  Sigframe <- data.frame(c(1:length(colnames(expr))))

  for(i in 1:length(PatternCom[,"Gene"])){
    entrID <- as.character(PatternCom[i,"Gene"])
    avrg <- mean(Biobase::exprs(eset)[which(rownames(eset)==entrID),])
    vector <- c()

    for(j in 1:length(colnames(expr))){
      level <- expr[entrID,j]

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
    colnames(newColumn)[1] = entrID
    Sigframe <- cbind(Sigframe, newColumn)
    rownames(Sigframe) = colnames(expr)
  }
  Sigframe <- Sigframe[-1]
  return(Sigframe)
}
