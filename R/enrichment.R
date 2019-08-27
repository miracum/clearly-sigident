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

goDiffReg_ <- function(mergeset, BatchRemovedExprs, design){
  edata <- BatchRemovedExprs
  # TODO is this always like this?
  rownames(edata) <- mergeset@featureData@data$ENTREZ_GENE_ID

  # run limma analysis
  return(fitLimma_(edata, design))
}

goEnrichmentAnalysis_ <- function(entrez, OrgDB, organism, fitlm, pathwayid, species){
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
  geneFC <- tt[,1:2]
  geneFC <- geneFC[order(geneFC$ID),]
  # removing empty Entrez-IDs
  #geneFC <- geneFC[27:701,] # old code
  geneFC <- geneFC[which(geneFC$ID != ""),]
  # removing Entrez-ID replicates
  geneFC <- geneFC[!duplicated(geneFC$ID),]
  rownames(geneFC) <- geneFC$ID
  geneFC$ID <- NULL

  # pathview
  p.out1 <- pathview::pathview(gene.data = geneFC, pathway.id = pathwayid, species = species, out.suffix = "image1")

  return(list(go = ego, kegg = kk, geneFC = geneFC, pathview = p.out1))
}
