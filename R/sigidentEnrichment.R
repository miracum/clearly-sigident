#' @title Perform Enrichtment Analysis in Gene Expression Datasets Derived from MicroArrays
#'
#' @description One function to perform gene enrichment analysis.
#'
#' @param mergedset A large merged Expression Data set. The output of the funtion `merge_()`.
#' @param species A character string indicating the sample's species. Currently supported: "Hs".
#' @param OrgDB A character string indicating the OrgDb. Currently supported: "org.Hs.eg.db".
#' @param organism A character string indicating the organism. Currently supported: "hsa".
#' @param pathwayid A character string indicating the pathway to show in the enrichment analysis. Currently supported: "hsa04110".'
#'
#' @inheritParams sigidentDEG
#'
#'
#' @export
sigidentEnrichment <- function(mergedset,
                               mergeset,
                               idtype,
                               design,
                               species,
                               OrgDB,
                               organism,
                               pathwayid,
                               plotdir = "./plots/",
                               csvdir = "./tables/"){

  stopifnot(
    class(mergedset) == "ExpressionSet",
    is.character(plotdir),
    is.character(csvdir),
    is.character(species),
    is.character(organism),
    is.character(OrgDB),
    is.character(pathwayid),
    idtype %in% c("entrez", "affy")
  )

  # create internal list for storage
  rv <- list()

  # store species, orgdb and orgamism
  rv$species <- species
  rv$orgdb <- OrgDB
  rv$organism <- organism
  rv$pathwayid <- pathwayid

  # store other arguments
  rv$idtype = idtype

  # store dirs
  rv$plotdir <- cleanPathName_(plotdir)
  rv$csvdir <- cleanPathName_(csvdir)

  # create output directories
  dir.create(rv$plotdir)
  dir.create(rv$csvdir)

  # add mergedset to list
  rv$mergedset <- mergedset
  # add mergeset to list
  rv$mergeset <- mergeset

  # gene enrichment
  rv$deg_entrez <- unique(mergedset@featureData@data$ENTREZ_GENE_ID)
  rv$deg_entrez <- rv$deg_entrez[rv$deg_entrez != ""]

  # test for over-representation of gene ontology terms
  rv$enr_topgo <- extractGOterms_(entrez = rv$deg_entrez,
                                  species = rv$species)
  data.table::fwrite(rv$enr_topgo, paste0(rv$csvdir, "Top_GO.csv"))

  # test for over-representation of KEGG pathways
  rv$enr_topkegg <- extractKEGGterms_(entrez = rv$deg_entrez,
                                      species = rv$species)
  data.table::fwrite(rv$enr_topkegg, paste0(rv$csvdir, "Top_KEGG.csv"))

  # take differential regulation between two groups (design) into account
  rv$enr_fitlm <- goDiffReg_(mergeset = rv$mergeset,
                             idtype = rv$idtype,
                             design = rv$design,
                             entrezids = rv$mergedset@featureData@data$ENTREZ_GENE_ID)

  # test for over-representation of gene ontology terms
  rv$enr_fitlm_topgo <- extractGOterms_(entrez = rv$enr_fitlm,
                                        species = rv$species,
                                        FDR = 0.01)
  data.table::fwrite(rv$enr_fitlm_topgo, paste0(rv$csvdir, "Top_GO_fitlm.csv"))


  # test for over-representation of KEGG pathways
  rv$enr_fitlm_topkegg <- extractKEGGterms_(entrez = rv$enr_fitlm,
                                            species = rv$species)
  data.table::fwrite(rv$enr_fitlm_topkegg, paste0(rv$csvdir, "Top_KEGG_fitlm.csv"))

  # perform enrichment analysis
  rv$enr_analysis <- goEnrichmentAnalysis_(entrez = rv$deg_entrez,
                                           OrgDB = rv$orgdb,
                                           organism = rv$organism,
                                           fitlm = rv$enr_fitlm,
                                           pathwayid = rv$pathwayid,
                                           species = rv$organism,
                                           plotdir = rv$plotdir)

  # plotting enrichmentanalysis
  createEnrichtedBarplot_(enrichmentobj = rv$enr_analysis$go,
                          type = "GO",
                          filename = paste0(rv$plotdir, "Enriched_GO.png"))

  createEnrichtedBarplot_(enrichmentobj = rv$enr_analysis$kegg,
                          type = "KEGG",
                          filename = paste0(rv$plotdir, "Enriched_KEGG.png"))
}
