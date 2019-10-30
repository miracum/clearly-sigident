#' @title extract_go_terms
#'
#' @description Helper function to extract GO terms
#'
#' @param gene A character vector containing the Entrez-IDs.
#' @param fdr The false discovery rate passed to `limma::goana`.
#'
#' @inheritParams sigidentEnrichment
#'
#' @export
extract_go_terms <- function(gene,
                             species,
                             fdr = NULL) {

  fit <- limma::topGO(limma::goana(
    de = gene,
    species = species,
    fdr = fdr
  ))
  return(fit)
}


#' @title extract_kegg_terms
#'
#' @description Helper function to extract KEGG terms
#'
#' @inheritParams extract_go_terms
#'
#' @export
extract_kegg_terms <- function(gene,
                               species) {

  fit <- limma::topKEGG(limma::kegga(de = gene,
                                     species = species))
  return(fit)
}


#' @title go_diff_reg
#'
#' @description Helper function to fitting linear models
#'
#' @param entrezids A character vector, containing entrez IDs.
#'   To be used only if 'idtype'= "affy" (default = NULL).
#'
#' @inheritParams identify_degs
#' @inheritParams batch_correction
#' @inheritParams sigidentDEG
#'
#' @export
go_diff_reg <- function(mergeset,
                        design,
                        idtype,
                        entrezids = NULL) {

  stopifnot(idtype %in% c("entrez", "affy"))

  if (idtype == "affy") {
    stopifnot(!is.null(entrezids))
    # map rownames to entrez-ids
    rownames(mergeset) <- entrezids
  }

  # run limma analysis
  fit <- limma_fitting(mergeset, design)
  return(fit)
}


#' @title go_enrichment_analysis
#'
#' @description Helper function to perform enrichment analysis
#'
#' @param org_db org_db
#' @param fitlm An MArrayLM object, returned by `go_diff_reg()`.
#' @param pathwayid KEGG pathway ID (like hsa04110 for cell cycle)
#'
#' @inheritParams sigidentEnrichment
#' @inheritParams extract_go_terms
#'
#' @export
go_enrichment_analysis <- function(gene,
                                   org_db,
                                   organism,
                                   fitlm,
                                   pathwayid,
                                   species,
                                   plotdir = NULL) {
  if (is.null(plotdir)) {
    plotdir <- "./plots/"
    if (!dir.exists("./plots/")) {
      dir.create("./plots/")
    }
  } else {
    plotdir <- clean_path_name(plotdir)
  }

  # TODO is this always like this?
  ego <- clusterProfiler::enrichGO(
    gene = gene,
    OrgDb = org_db,
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  kk <- clusterProfiler::enrichKEGG(
    gene = gene,
    organism = organism,
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01
  )

  # create matrix with Entrez-IDs and logFC from limma
  tt <- limma::topTable(
    fit = fitlm,
    coef = 2,
    number = Inf,
    p.value = 0.05,
    lfc = 2
  )
  # gene_fc <- as.data.frame(cbind(ID = rownames(tt),
  # logFC = as.numeric(tt[,"logFC"])))
  gene_fc <- tt[, 1:2]
  gene_fc <- gene_fc[order(gene_fc$ID), ]
  # removing empty Entrez-IDs
  gene_fc <- gene_fc[which(gene_fc$ID != ""), ]
  # removing Entrez-ID replicates
  gene_fc <- gene_fc[!duplicated(gene_fc$ID), ]
  rownames(gene_fc) <- gene_fc$ID
  gene_fc[, "ID"] <- NULL

  # pathview
  # https://github.com/egeulgen/pathfindR/issues/10
  utils::data("gene.idtype.bods", package = "pathview")
  # workaround to set correct workingdir for pathview
  oldwd <- getwd()
  setwd(paste0(oldwd, "/", plotdir))
  p_out1 <- pathview::pathview(gene.data = gene_fc,
                               pathway.id = pathwayid,
                               species = organism)
  setwd(oldwd)

  return(list(
    go = ego,
    kegg = kk,
    gene_fc = gene_fc,
    pathview = p_out1
  ))
}
