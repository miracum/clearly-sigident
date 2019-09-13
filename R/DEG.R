DEG.limma_ <- function(mergeset, design){
  fit <- fitLimma_(mergeset, design)
  t <- limma::topTable(fit, coef=2,number=Inf,p.value=0.05,lfc=2)
  genes <- rownames(t)
  return(genes)
}

identify.DEGs_ <- function(mergeset, design, qValue){
  fit <- fitLimma_(mergeset, design)
  t <- limma::topTable(fit, coef=2,number=Inf,p.value=qValue, lfc=2, adjust.method = "BH")
  genes <- rownames(t)
  return(genes)
}

fitLimma_ <- function(mergeset, design){
  return(limma::eBayes(limma::lmFit(mergeset, design)))
}

limmaTopTable_ <- function(mergeset, design, qValue){
  fit <- fitLimma_(mergeset, design)
  t <- limma::topTable(fit,
                       coef = 2,
                       number = Inf,
                       p.value = qValue,
                       lfc = 2,
                       adjust.method = "BH")
  t[,2:4] <- NULL
  t[,3] <- NULL
  t <- cbind(rownames(t),t)
  colnames(t) <- c("Probe ID", "logFC", "adj.qValue")
  return(t)
}

exportDEGannotations_ <- function(mergeset, genes){
  ids <- genes
  # TODO is this always so?
  sym <- Biobase::fData(mergeset)["Gene Symbol"][ids,]
  tit <- Biobase::fData(mergeset)["Gene Title"][ids,]
  gbACC <- Biobase::fData(mergeset)["GB_ACC"][ids,]
  Entrez <- Biobase::fData(mergeset)["ENTREZ_GENE_ID"][ids,]
  DEGsInfo <- data.table::data.table(cbind(ids,
                                           sym,
                                           tit,
                                           gbACC,
                                           Entrez))
  colnames(DEGsInfo) <- c("probe_ID","gene_symbol","gene_title","genebank_accession","entrez_id")
  return(DEGsInfo)
}


qSelection_ <- function(sampleMetadata, deg.q.selection = NULL){
  if (is.null(deg.q.selection)){
    deg_q <- 1/length(sampleMetadata$sample)
  } else {
    deg_q <- as.numeric(deg.q.selection)
  }
  return(deg_q)
}
