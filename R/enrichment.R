enrichmentDataSelection_ <- function(mergeset, genes){
  # TODO ist das immer Spalte 12?
  Entrez <- Biobase::fData(mergeset)[12][genes,]
  # remove empty characters and replicates in EntrezIDs
  DEGs.Entrez <- Entrez[Entrez != ""]
  return(unique(DEGs.Entrez))
}


extractGOterms_ <- function(entrez, species, FDR = NULL){
  return(limma::topGO(
    limma::goana(de = entrez, species = species, FDR = FDR)
  ))
}


extractKEGGterms_ <- function(entrez, species){
  return(limma::topKEGG(
    limma::kegga(de = entrez, species = species)
  ))
}

goDiffReg_ <- function(mergeset, design){
  # run limma analysis
  return(fitLimma_(mergeset, design))
}

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
