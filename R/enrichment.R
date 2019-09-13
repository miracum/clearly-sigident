#' @title extractGOterms_
#'
#' @description Helper function to extract GO terms
#'
#' @param entrez A character vector containing the Entrez-IDs.
#' @param FDR The false discovery rate passed to `limma::goana`.
#'
#' @inheritParams sigidentMicroarray
#'
#' @export
extractGOterms_ <- function(entrez, species, FDR = NULL){
  return(limma::topGO(
    limma::goana(de = entrez, species = species, FDR = FDR)
  ))
}


#' @title extractKEGGterms_
#'
#' @description Helper function to extract KEGG terms
#'
#' @inheritParams extractGOterms_
#'
#' @export
extractKEGGterms_ <- function(entrez, species){
  return(limma::topKEGG(
    limma::kegga(de = entrez, species = species)
  ))
}


#' @title goDiffReg_
#'
#' @description Helper function to ... (TODO what do we do here?)
#'
#' @inheritParams identifyDEGs_
#'
#' @export
goDiffReg_ <- function(mergeset, design){
  # run limma analysis
  return(fitLimma_(mergeset, design))
}


#' @title goEnrichmentAnalysis_
#'
#' @description Helper function to perform enrichment analysis
#'
#' @param fitlm A object, returned by `goDiffReg_()`.
#'
#' @inheritParams sigidentMicroarray
#' @inheritParams extractGOterms_
#'
#' @export
goEnrichmentAnalysis_ <- function(entrez, OrgDB, organism, fitlm, pathwayid, species, plotdir = NULL){

  if (is.null(plotdir)){
    plotdir <- "./plots/"
    if (!dir.exists("./plots/")){
      dir.create("./plots/")
    }
  } else {
    plotdir <- cleanPathName_(plotdir)
  }

  # TODO is this always like this?
  ego <- clusterProfiler::enrichGO(gene = entrez,
                                   OrgDb = OrgDB,
                                   keyType = "ENTREZID",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.01,
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)

  kk <- clusterProfiler::enrichKEGG(gene = entrez,
                                    organism = organism,
                                    keyType = "kegg",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.01)

  # create matrix with Entrez-IDs and logFC from limma
  tt <- limma::topTable(fitlm,
                        coef=2,
                        number=Inf,
                        p.value=0.05,
                        lfc=2)
  geneFC <- as.data.frame(cbind(ID = rownames(tt), logFC = as.numeric(tt[,"logFC"])))
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
                               species = species)
  setwd(oldwd)

  return(list(go = ego, kegg = kk, geneFC = geneFC, pathview = p.out1))
}
