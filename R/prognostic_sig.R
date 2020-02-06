#' @title get_survival_time
#'
#' @description Helper function to get survival time.
#'
#' @param discoverystudies_w_timedata A list that contains specifications
#'   on the study/studies that contain(s) survival time information.
#' @param datadir A character string. Path to the data-folder inside the
#'   metadata folder.
#'
#' @inheritParams sigidentPrognostic
#'
#' @export
get_survival_time <- function(sample_metadata,
                              genes,
                              idtype,
                              discoverystudies_w_timedata,
                              datadir) {

  stopifnot(
    is.list(discoverystudies_w_timedata)
  )

  targetcol <- "target"
  controlname <- "Control"
  targetname <- "Target"

  outlist <- list()

  for (st in names(discoverystudies_w_timedata)) {
    stopifnot(
      is.character(discoverystudies_w_timedata[[st]]$timecol),
      is.list(discoverystudies_w_timedata[[st]]$status),
      is.character(discoverystudies_w_timedata[[st]]$status$statuscol),
      is.list(discoverystudies_w_timedata[[st]]$status$levels),
      is.character(discoverystudies_w_timedata[[st]]$status$levels$alive),
      is.character(discoverystudies_w_timedata[[st]]$status$levels$deceased),
      is.character(discoverystudies_w_timedata[[st]]$targetcolname),
      is.logical(discoverystudies_w_timedata[[st]]$use_rawdata),
      is.numeric(discoverystudies_w_timedata[[st]]$setid)
    )

    tumorsamples <- sample_metadata[eval(
      parse(text = paste0("sample_metadata$", targetcol))
    ) == targetname, ]

    # setd use_raw, if not provided with function arguments
    use_raw <- ifelse(
      is.null(discoverystudies_w_timedata[[st]]$use_rawdata),
      FALSE,
      TRUE
    )

    eset <- tryCatch(
      expr = {
        eset <- eval(parse(text = st),
                     envir = 1L)
        cat(paste0("\nLoaded ", st, " from .Globalenv...\n"))
        eset
      }, error = function(e) {
        eset <- sigident.preproc::geo_load_eset(
          name = st,
          datadir = datadir,
          targetcolname = discoverystudies_w_timedata[[st]]$targetcolname,
          targetcol = targetcol,
          targetname = targetname,
          controlname = controlname,
          targetlevelname = discoverystudies_w_timedata[[st]]$targetlevelname,
          controllevelname = discoverystudies_w_timedata[[st]]$controllevelname,
          use_rawdata = use_raw,
          setid = discoverystudies_w_timedata[[st]]$setid
        )
        cat(paste0("\nLoaded ",
                   st,
                   " from URL\n"))
        eset
      }, finally = function(f) {
        return(eset)
      }
    )

    # filter only data of tumor samples
    if (!is.null(discoverystudies_w_timedata[[st]]$targetlevelname)) {
      eset_targets <-
        eset[, which(eset$geo_accession %in% tumorsamples$sample)]
    } else {
      # use the whole dataset
      eset_targets <- eset
    }

    timestatus <- extract_time_status(
      eset_targets = eset_targets,
      timecol = discoverystudies_w_timedata[[st]]$timecol,
      statuscol = discoverystudies_w_timedata[[st]]$status$statuscol,
      statuslevels = discoverystudies_w_timedata[[st]]$status$levels
    )
    time <- timestatus$time
    status <- timestatus$status

    # generating a table: columns: survivaltime, status;
    # rows: individual samples
    survival_table <- data.frame(time, status)

    expr <- sigident.preproc::geo_create_expressionset(
      eset = eset_targets,
      idtype = idtype
    )

    #% ids <- rownames(expr)

    for (deg in genes) {
      # loop to bind a new column filled with exprs-values for each DEG
      exp <- exprs_vector(expr[deg, ])
      df <- data.frame(exp)
      colnames(df) <- deg
      # df without batch effect removal, alternative would be rmBatch
      survival_table <- cbind(survival_table, df)
    }

    rownames(survival_table) <- colnames(expr)

    outlist[[st]] <- survival_table

  }
  return(outlist)
}


#' @title univ_cox
#'
#' @description Helper function to compute univariate cox regression and
#'   determine significance of each gene through separate univariate Cox
#'   regressions
#'
#' @param survtable A data.frame. Output of the function `get_survival_time()`.
#'
#' @inheritParams sigidentPrognostic
#'
#' @export
univ_cox <- function(survtable, genes) {

  covariates <- colnames(survtable)[-c(1:2)]
  covariates <- gsub("[[:punct:]]", "", covariates)
  covariates <- paste("ID", covariates)
  covariates <- gsub(" ", "", covariates)

  colnames(survtable)[-c(1:2)] <- covariates

  univ_formulas <- sapply(covariates, function(x) {
    # build separate formula for each variable
    stats::as.formula(paste("survival::Surv(time, status) ~", x))
  })
  univ_models <- lapply(univ_formulas, function(f) {
    survival::coxph(formula = f, data = survtable)
  }) # build cox model for each variable separatly

  # parallel:
  #% univ_models <- mcmapply(univ_formulas, FUN =  function(f) {
  #%   survival::coxph(formula = f, data = survtable)},
  #%   mc.cores = parallel::detectCores()-1, mc.preschedule = TRUE)

  # extract results
  univ_results <- lapply(
    univ_models,
    function(x) {
      x <- summary(x)
      p_value <- signif(x$wald["pvalue"], digits = 2)
      wald_test <- signif(x$wald["test"], digits = 2)
      beta <- signif(x$coef[1], digits = 2)
      #coeficient beta # this block applies a function to every DEG,
      hr <- signif(x$coef[2], digits = 2)
      #% exp(beta) # executing a cox regression analysis
      # DEG ~ survivaltime
      hr_confint_lower <- signif(x$conf.int[, "lower .95"], 2)
      hr_confint_upper <- signif(x$conf.int[, "upper .95"], 2)

      hr <- paste0(hr,
                   " (",
                   hr_confint_lower,
                   "-",
                   hr_confint_upper,
                   ")")

      res <- c(beta, hr, wald_test, p_value) #
      names(res) <- c("beta",
                      "HR (95% CI for HR)",
                      "wald_test",
                      "p_value")

      return(res) #
      #% return(exp(cbind(coef(x),confint(x)))) #
    }) #

  #formatting of result
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  result <- as.data.frame(res)
  # reforamtting of p.value-column
  result$p_value <- as.character(result$p_value)
  result$p_value <- as.numeric(result$p_value)
  # selecting indices with a significant p-value
  indices <- which(result$p_value <= 0.05)
  # IDs now contains Entrez-IDs of all genes with significant p-values
  ids <- genes[indices]

  # TODO symbols?

  sig_cov <- result[with(result, which(result$p_value <= 0.05)), ]
  rownames(sig_cov) <- ids
  sig_cov <- cbind(id = ids, sig_cov)

  return(sig_cov)
}


exprs_vector <- function(deg) {
  vector <- c() # creates empty vector
  for (i in seq_len(length(deg))) {
    # vector is filled by columns of exprs (in sample order)
    vector <- append(vector, deg[[i]])
  }
  # vector now contains one row (of the called gene) of exprs in
  # sample/column order
  return(vector)
}


#' @title generate_expression_pattern
#'
#' @description Helper function to generate expression pattern
#'
#' @param classifier_studies A character vector. Names of the studies
#'   used to train the classifier.
#' @param sig_cov A data.frame. Output of the function `univ_cox()`.
#'
#' @inheritParams sigidentPrognostic
#'
#' @export
generate_expression_pattern <- function(classifier_studies,
                                        sig_cov,
                                        mergeset,
                                        sample_metadata) {

  dd <- sigident.preproc::geo_create_diagnosisbatch(
    sample_metadata =
      sample_metadata[sample_metadata$study %in% classifier_studies, ]
  )
  diagnosis <- dd$diagnosis


  # get sample names (= colnames of mergeset)
  classifier_colnames <-
    sample_metadata[sample_metadata$study %in%
                      classifier_studies, get("sample")]

  # reduce mergedata to use only the data specified in classifier_studies
  classifier_data <- mergeset[, classifier_colnames]

  # generate vectors 'control' and 'tumor' for function expressionPattern
  control <- ctrl_vec(classifier_data, diagnosis)
  tumor <- tum_vec(classifier_data, diagnosis)

  ids <- sig_cov[, "id"]

  # train the prognostic classifier
  pattern <- c()
  for (id in ids) {
    pattern <-
      append(
        pattern,
        expression_pattern(
          mergeset = classifier_data,
          ids = id,
          tumor = tumor,
          control = control
        )
      )
  }
  pattern_com <- data.frame(ids, pattern)

  colnames(pattern_com) <- c("Gene", "Over/Under")

  return(pattern_com)
}




ctrl_vec <- function(mergeset,
                     diagnosis) {
  t1 <- colnames(mergeset)
  t1 <- cbind(t1, diagnosis)
  t1 <- as.data.frame(t1)
  ctrl_m <- t1[with(t1, which(diagnosis == 0)), ]
  control <- as.vector(ctrl_m$t)
  return(control)
}

tum_vec <- function(mergeset,
                    diagnosis) {
  t1 <- colnames(mergeset)
  t1 <- cbind(t1, diagnosis)
  t1 <- as.data.frame(t1)
  tum_m <- t1[with(t1, which(diagnosis == 1)), ]
  tumor <- as.vector(tum_m$t)
  return(tumor)
}



expression_pattern <- function(mergeset, ids, tumor, control) {
  # selecting control-samples and computing mean
  ctrl <- mean(mergeset[ids, control])
  # selecting tumor-samples and computing mean
  tum <- mean(mergeset[ids, tumor])

  if (ctrl < tum) {
    pattern <- "Over"

  } else if (ctrl > tum) {
    # Over-/Underexpression is checked by comparing group-means
    pattern <- "Under" # vector is filled with respective symbols

  } else {
    pattern <- NA
  }
  return(pattern) # pattern vector is returned
}



#' @title prognostic_classifier
#'
#' @description Helper function to perform prognostic classification.
#'
#' @param validationstudiesinfo A list that contains specifications on
#'   the study that contains the validation information.
#' @param pattern_com A data.frame. The output of the function
#'   `generate_expression_pattern()`.
#'
#' @inheritParams get_survival_time
#'
#' @export
prognostic_classifier <- function(pattern_com,
                                  idtype,
                                  validationstudiesinfo,
                                  datadir) {

  targetcol <- "target"
  controlname <- "Control"
  targetname <- "Target"

  outlist <- list()

  for (st in names(validationstudiesinfo)) {

    # setd use_raw, if not provided with function arguments
    use_raw <- ifelse(
      is.null(validationstudiesinfo[[st]]$use_rawdata),
      FALSE,
      TRUE
    )

    eset <- tryCatch(
      expr = {
        eset <- eval(parse(text = st),
                     envir = 1L)
        cat(paste0("\nLoaded ", st, " from .Globalenv...\n"))
        eset
      }, error = function(e) {
        eset <- sigident.preproc::geo_load_eset(
          name = st,
          datadir = datadir,
          targetcolname = validationstudiesinfo[[st]]$targetcolname,
          targetcol = targetcol,
          targetname = targetname,
          controlname = controlname,
          targetlevelname = validationstudiesinfo[[st]]$targetlevelname,
          controllevelname = validationstudiesinfo[[st]]$controllevelname,
          use_rawdata = use_raw,
          setid = validationstudiesinfo[[st]]$setid
        )
        cat(paste0("\nLoaded ",
                   st,
                   " from URL\n"))
        eset
      }, finally = function(f) {
        return(eset)
      }
    )

    #% diagnosis <- create_diagnosis(vector = eset[[targetcol]],
    #%                               targetname = targetname,
    #%                               controlname = controlname)

    # generate vectors 'control' and 'tumor' for function expressionPattern
    #% classifier_data <- Biobase::exprs(eset)
    #% control <- ctrl_vec(classifier_data, diagnosis)
    #% tumor <- tum_vec(classifier_data, diagnosis)

    # filter only data of tumor samples
    if (!is.null(validationstudiesinfo[[st]]$targetlevelname)) {
      eset_targets <- eset[, which(eset[[targetcol]] == targetname)]
    } else {
      # use the whole dataset
      eset_targets <- eset
    }

    timestatus <- extract_time_status(
      eset_targets = eset_targets,
      timecol = validationstudiesinfo[[st]]$timecol,
      statuscol = validationstudiesinfo[[st]]$status$statuscol,
      statuslevels = validationstudiesinfo[[st]]$status$levels
    )


    time <- timestatus$time
    status <- timestatus$status
    risk_table <- data.frame(time, status)

    expr <- sigident.preproc::geo_create_expressionset(
      eset = eset_targets,
      idtype = idtype
    )

    # classification and Kaplan-Meier estimator
    sig <- sig_analysis(expr = expr, pattern_com = pattern_com)
    risk_table_sig <- classify(sig)

    # prepare classification and required data for Kaplan-Meier estimator
    groups <- risk_table_sig[, "RiskGroup"]
    risk_table <- cbind(risk_table, groups)
    colnames(risk_table) <- c("time", "status", "groups")
    rownames(risk_table) <- colnames(expr)

    kap <- fit_kaplan_estimator(risktable = risk_table)

    # prepare classification and required data for Kaplan-Meier estimator
    outlist[[st]] <- list(kaplan_estimator = kap,
                          risktable = risk_table)

  }
  return(outlist)
}

fit_kaplan_estimator <- function(risktable) {

  # create new data
  new_df <- with(risktable,
                 data.frame(
                   groups = c(0, 1),
                   survival_time = rep(mean(time, na.rm = TRUE), 2),
                   ph.ecog = c(1, 1)
                 ))

  # fit proportional hazards regression model
  fit <- survival::survfit(
    formula = survival::coxph(
      survival::Surv(time, status) ~ groups,
      data = risktable
    ),
    newdata = new_df # argument, passed to survfit.coxph
  )
  # return model (this is redundant; however, if you provide this
  # res_cox object to the formula argument in survfit, the ggsurvplot-
  # function will fail)
  res_cox <- survival::coxph(
    survival::Surv(time, status) ~ groups,
    data = risktable
  )
  return(list(res_cox = res_cox,
              fit = fit))
}


classify <- function(sigtable) {
  table <- sigtable
  rsums <- rowSums(table)
  vector <- c()
  for (i in seq_len(length(rsums))) {
    if (rsums[i] > 0.5 * length(colnames(table))) {
      a <- 1
    }
    else{
      a <- 0
    }
    vector <- append(vector, a)
  }
  new_column <- data.frame(vector)
  colnames(new_column)[1] <- "RiskGroup"
  table <- cbind(table, new_column)
  return(table)
}


extract_time_status <- function(eset_targets,
                                timecol,
                                statuscol,
                                statuslevels) {
  #% time <- eset_targets[, get("timecol")]
  # does not work with expression sets
  time <- eval(parse(text = paste0("eset_targets$", timecol)))
  time <- as.character(time)
  time[which(time == "")] <- "time: not available"
  time <- unlist(lapply(time, strsplit, ": ")) # string is split

  # first part of string ("overall survival") is discarded
  time <- time[c(FALSE, TRUE)]
  # results in a vector of survival times (in sample order)
  time <- as.numeric(time)

  # coding if patient survived or not
  status <- eval(parse(text = paste0("eset_targets$", statuscol)))

  # 1 = no event/alive; 2 = event/death
  if (is.na(statuslevels$na)) {
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
  return(list(time = time, status = status))
}



over <- function(avrg, level) {
  if (avrg < level) {
    a <- 1
  }
  else{
    a <- 0
  }
  return(a)
}

# defining functions to check for Over-/Underexpression
under <- function(avrg, level) {
  # if expression pattern is matching the signature,
  # 1 is returned, otherwise 0
  if (avrg > level) {
    a <- 1
  }
  else{
    a <- 0
  }
  return(a)
}

sig_analysis <- function(expr, pattern_com) {
  # Input is an eset, consisting of tumor-samples only and the
  # signature to be tested

  sigframe <- data.frame(c(seq_len(length(colnames(expr)))))

  for (i in seq_len(length(pattern_com[, "Gene"]))) {
    id <- as.character(pattern_com[i, "Gene"])
    avrg <- mean(expr[which(rownames(expr) == id), ])
    vector <- c()

    for (j in seq_len(length(colnames(expr)))) {
      level <- expr[id, j]

      if (pattern_com[i, "Over/Under"] == "Over") {
        a <- over(avrg, level)
      }
      else if (pattern_com[i, "Over/Under"] == "Under") {
        a <- under(avrg, level)
      }
      else{
        print("Wrong formatting of Signature")
      }
      vector <- append(vector, a)
    }

    new_column <- data.frame(vector)
    colnames(new_column)[1] <- id
    sigframe <- cbind(sigframe, new_column)
    rownames(sigframe) <- colnames(expr)
  }
  sigframe <- sigframe[-1]
  return(sigframe)
}
