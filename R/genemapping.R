geneMap_ <- function(eset, genes){
  symbols <- list()
  for(i in 1:length(genes))
  {
    exrow <- which(rownames(Biobase::exprs(eset))==genes[i]) # which row contains DEG of interest
    # TODO is this always like this?
    symbol <- Biobase::fData(eset)$"Gene Symbol"[exrow] # getting the respective Gene symbol from fData
    symbol <- as.character(symbol)
    symbols <- c(symbols,symbol)
  }
  symbols <- as.character(symbols) # what to do with multiple Gene Symbols/additional annotation??
  return(symbols)
}


geneMapping_ <- function(mergeset, genes){
  syms <- geneMap_(mergeset,genes)

  # TODO is this always like this?
  genelist <- unlist(lapply(syms, strsplit, " /// "))
  genelist <- as.factor(genelist)
  genelist <- unique(genelist)
  genelist <- as.character(genelist)
  return(genelist)
}
