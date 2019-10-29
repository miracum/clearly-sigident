deg_limma <- function(mergeset,
                       design) {

  fit <- limma_fitting(mergeset, design)

  t <- limma::topTable(
    fit,
    coef = 2,
    number = Inf,
    p.value = 0.05,
    lfc = 2
  )
  genes <- rownames(t)
  return(genes)
}

#' @title identify_degs
#'
#' @description Helper function to identify DEGs based on the limma package
#'
#' @param q_value A numeric value for the q-value (false discovery rate)
#'   (default=0.01).
#'
#' @inheritParams sigidentDEG
#'
#' @export
identify_degs <- function(mergeset,
                          design,
                          q_value = 0.01) {

  stopifnot(
    q_value > 0 | q_value <= 0.05,
    is.numeric(q_value)
  )

  fit <- limma_fitting(mergeset, design)

  t <- limma::topTable(
    fit,
    coef = 2,
    number = Inf,
    p.value = q_value,
    lfc = 2,
    adjust.method = "BH"
  )
  genes <- rownames(t)
  return(genes)
}

limma_fitting <- function(mergeset,
                      design) {

  fit <- limma::eBayes(limma::lmFit(mergeset, design))

  return(fit)
}


#' @title limma_top_table
#'
#' @description Helper function to get DEG results
#'
#' @inheritParams sigidentDEG
#' @inheritParams identify_degs
#'
#' @export
limma_top_table <- function(mergeset,
                           design,
                           q_value) {

  fit <- limma_fitting(mergeset, design)

  t <- limma::topTable(
    fit,
    coef = 2,
    number = Inf,
    p.value = q_value,
    lfc = 2,
    adjust.method = "BH"
  )

  t[, 2:4] <- NULL
  t[, 3] <- NULL
  t <- cbind(rownames(t), t)
  colnames(t) <- c("Probe ID", "logFC", "adj.q_value")

  return(t)
}

# #' @title qSelection_
# #'
# #' @description Helper function to select q_values
# #'
# #' @inheritParams sigidentDEG
# #'
# #' @export
#% qSelection_ <- function(sample_metadata,
#%                         study_metadata,
#%                         deg.q.selection = NULL) {
#%
#%   if (is.null(deg.q.selection)) {
#%     discovery <- discovery_func(
#%       sample_metadata = sample_metadata,
#%       study_metadata = study_metadata
#%     )
#%     deg_q <- 1/length(
#%       sample_metadata[sample_metadata$study %in% discovery, "sample"]
#%     )
#%   } else {
#%     deg_q <- as.numeric(deg.q.selection)
#%   }
#%   return(deg_q)
#% }
