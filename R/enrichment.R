#' @title extractGOterms_
#'
#' @description Helper function to extract GO terms
#'
#' @param gene A character vector containing the Entrez-IDs.
#' @param FDR The false discovery rate passed to `limma::goana`.
#'
#' @inheritParams sigidentEnrichment
#'
#' @export
extractGOterms_ <- function(gene, species, FDR = NULL){
  fit <- limma::topGO(
    limma::goana(de = gene, species = species, FDR = FDR)
  )
  return(fit)
}


#' @title extractKEGGterms_
#'
#' @description Helper function to extract KEGG terms
#'
#' @inheritParams extractGOterms_
#'
#' @export
extractKEGGterms_ <- function(gene, species){
  fit <- limma::topKEGG(
    limma::kegga(de = gene, species = species)
  )
  return(fit)
}


#' @title goDiffReg_
#'
#' @description Helper function to fitting linear models
#'
#' @param entrezids A character vector, containing entrez IDs. To be used only if 'idtype'= "affy" (default = NULL).
#'
#' @inheritParams identifyDEGs_
#' @inheritParams batchCorrection_
#' @inheritParams sigidentDEG
#'
#' @export
goDiffReg_ <- function(mergeset, design, idtype, entrezids = NULL){
  stopifnot(
    idtype %in% c("entrez", "affy")
  )
  if (idtype == "affy"){
    stopifnot(!is.null(entrezids))
    # map rownames to entrez-ids
    rownames(mergeset) <- entrezids
  }
  # run limma analysis
  fit <- fitLimma_(mergeset, design)
  return(fit)
}


#' @title goEnrichmentAnalysis_
#'
#' @description Helper function to perform enrichment analysis
#'
#' @param OrgDB OrgDB
#' @param fitlm An MArrayLM object, returned by `goDiffReg_()`.
#' @param pathwayid KEGG pathway ID (like hsa04110 for cell cycle)
#'
#' @inheritParams sigidentEnrichment
#' @inheritParams extractGOterms_
#'
#' @export
goEnrichmentAnalysis_ <- function(gene, OrgDB, organism, fitlm, pathwayid, species, plotdir = NULL){

  if (is.null(plotdir)){
    plotdir <- "./plots/"
    if (!dir.exists("./plots/")){
      dir.create("./plots/")
    }
  } else {
    plotdir <- cleanPathName_(plotdir)
  }

  # TODO is this always like this?
  ego <- clusterProfiler::enrichGO(gene = gene,
                                   OrgDb = OrgDB,
                                   keyType = "ENTREZID",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.01,
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)

  kk <- clusterProfiler::enrichKEGG(gene = gene,
                                    organism = organism,
                                    keyType = "kegg",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.01)

  # create matrix with Entrez-IDs and logFC from limma
  tt <- limma::topTable(fit = fitlm,
                        coef = 2,
                        number = Inf,
                        p.value = 0.05,
                        lfc = 2)
  #geneFC <- as.data.frame(cbind(ID = rownames(tt), logFC = as.numeric(tt[,"logFC"])))
  geneFC <- tt[,1:2]
  geneFC <- geneFC[order(geneFC$ID),]
  # removing empty Entrez-IDs
  geneFC <- geneFC[which(geneFC$ID != ""),]
  # removing Entrez-ID replicates
  geneFC <- geneFC[!duplicated(geneFC$ID),]
  rownames(geneFC) <- geneFC$ID
  geneFC$ID <- NULL

  # pathview
  # https://github.com/egeulgen/pathfindR/issues/10
  require(pathfindR)
  # workaround to set correct workingdir for pathview
  oldwd <- getwd()
  setwd(paste0(oldwd, "/", plotdir))
  p.out1 <- pathview::pathview(gene.data = geneFC,
                               pathway.id = pathwayid,
                               species = organism)
  setwd(oldwd)

  return(list(go = ego, kegg = kk, geneFC = geneFC, pathview = p.out1))
}
