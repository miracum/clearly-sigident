loadEset_ <- function(name, datadir, targetcolname, targetcol, targetname, controlname, targetlevelname, controllevelname, setid = 1){
  # original GEO data
  eset <- GEOquery::getGEO(name, destdir = datadir)[[setid]]
  # rename targetcol
  colnames(Biobase::pData(eset))[which(colnames(Biobase::pData(eset)) == targetcolname)] <- targetcol

  # rename levels of targetcol
  if (!is.null(targetlevelname)){
    levelnames <- c(targetname, controlname)
    names(levelnames) <- c(targetlevelname,
                           controllevelname)
    eset[[targetcol]] <- plyr::revalue(eset[[targetcol]], levelnames)
  }
  return(eset)
}

createExpressionSet_ <- function(eset, idtype){
  expr <- Biobase::exprs(eset)
  expr <- idType_(expr = expr, eset = eset, idtype = idtype)
  return(expr)
}


idType_ <- function(expr, eset, idtype){
  stopifnot(
    idtype %in% c("entrez", "affy")
  )
  if (idtype == "entrez"){
    rownames(expr) <- as.character(eset@featureData@data$ENTREZ_GENE_ID)
    # remove empty characters and replicates in EntrezIDs
    expr <- expr[rownames(expr)!="",]
    expr <- expr[!duplicated(rownames(expr)),]
  } else if (idtype == "affy") {
    rownames(expr) <- as.character(eset@featureData@data$ID)
  }
  return(expr)
}
