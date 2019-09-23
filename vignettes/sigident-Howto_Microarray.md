---
title: "sigident - Howto Microarray"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{sigident-Howto_Microarray}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::knitr}
editor_options: 
  chunk_output_type: console
---




```r
library(sigident)
library(knitr)
```

# Prerequisites and Installation 

## Install sigident package


```r
options('repos' = 'https://ftp.fau.de/cran/')
install.packages("devtools")
devtools::install_git("https://gitlab.miracum.org/clearly/sigident.git")
```

# Data Preparation 

In order to use this R package and its functions, you need to prepare a merged gene expression dataset. 
This example uses the GEO lung cancer studies "GSE18842", "GSE19804" and "GSE19188", which contain in total 367 samples (197 tumor; 170 non-tumor) and 54,675 transcripts gained by using Affymetrix GeneChip Human Genome U133 Plus 2.0 Array (platform GPL570).

## Download of GEO datasets and normalization

To use microarray datasets from the GEO database the studies can be downloaded directly executing R code. Therefore, initializing the directory and creating the required folders is needed before downloading GEO studies. The objects `targetcol` should also be entitled "target". To run this package analyzing DEGs and identifying diagnostic and prognostic gene signatures for the seperation of two groups (or subtypes, healthy vs cancer, etc) we here named `controlname` and `targetname` exemplary as "control" and "Lung Cancer".  


```r
# initialize filePath:
filePath <- tempdir()
maindir <- "./geodata/"
datadir <- paste0(maindir, "data/")
dir.create(maindir)
dir.create(datadir)

# initialize more infos on the study
targetcol <- "target" # should be named "target"
controlname <- "Control"
targetname <- "Lung Cancer"
```

Insert the GEO accesion numbers of the desired GEO studies into the `getGEO()` function. This will download the datasets containing expressin data, pheno and feature data as well as the annotation of the probe sets and store it as .txt.gz files in the "data" folder. 


```r
GSE18842 <- GEOquery::getGEO("GSE18842", destdir = datadir)
GSE19804 <- GEOquery::getGEO("GSE19804", destdir = datadir)

# extracting expressionSets
eset1 <- GSE18842[[1]]
eset2 <- GSE19804[[1]]
```

Morover, the approach enables to download raw data as CEL files, uncompress it and conduct GCRMA normalization. For this, further commands are necessary to subsequently import the CEL files into the R environment. Accordingly `esetC` represents a normalized epressionSet though only containing expression data without pheno and feature data. These data can be added later. (Don't be confused by our example. We could have downloaded the study GSE19188 directly by applying `getGEO()` instead of downlaoding the raw data and adding pheno and feature data afterwards. We just wanted to demonstrating you how to use raw data in form of CEL files)


```r
# download raw data, uncompress them and execute GCRMA normalization
# actually not necessary for GSE19188 but to showcase feasibility of dealing with .CEL files
GEOquery::getGEOSuppFiles("GSE19188", baseDir = filePath)
utils::untar(paste0(filePath, "/GSE19188/GSE19188_RAW.tar"), exdir=maindir)
cels <- list.files(maindir, pattern = "[gz]")
tryCatch({
  sapply(paste(maindir, cels, sep="/"), GEOquery::gunzip)
  # at this point, first the textfile "phenodata.txt must be generated
  ## For downloading raw data use the function: getGEOSuppFiles(""). After downloading and uncompressing the raw file (*_RAW.tar) we need to capture
  ## the experimental information. Therefore a text file is created in the terminal window or rather in the command line.
  system("bash -c 'ls ./geodata/*.CEL > ./geodata/phenodata.txt'")
  ## The command: ls data/*.CEL > data/phenodata.txt puts all .CEl files in one text file.
  ## After this step the normalization can be started.
}, error = function(e){
  print(e)
})

cels <- list.files(maindir, pattern = "CEL$")

celfiles <- paste0(maindir, cels)
# use justGCRMA now due to memory issues with gcrma --> justGCRMA is much more memory efficient
# https://www.rdocumentation.org/packages/gcrma/versions/2.44.0/topics/justGCRMA
esetC <- gcrma::justGCRMA(filenames = celfiles, fast = TRUE)
# call garbage collection here
gc()

# old code
#celfiles <- affy::ReadAffy(filenames = gsub("\\.gz", "", cels), celfile.path = "./geodata")
#esetC <- gcrma::gcrma(celfiles, fast = TRUE)

# esetC = now is a normalized expressionSet, but without f- and pData
# the f- and pData are transfered manually from dataset GSE19188

GSE19188 <- GEOquery::getGEO("GSE19188", destdir = datadir)

eset3 <- GSE19188[[1]]
pData3 <- Biobase::pData(eset3)
fData3 <- Biobase::fData(eset3)
Biobase::pData(esetC) <- pData3
Biobase::fData(esetC) <- fData3
colnames(Biobase::exprs(esetC)) <- colnames(Biobase::exprs(eset3))
Biobase::annotation(esetC) <- Biobase::annotation(eset3)
eset3 <- esetC
```

## Formatting the phenoData 

The phenoData of the expresssionSets need to be formatted that the phenoData columns overlap perfectly. This means that all expressionSets should contain identical columns in the same order. Not consistent columns need to get removed. For a simpler procedure we recommend only retaining phenoData columns required for the downstream analysis (e.g. targetcol and survival data). 


```r
# create version b, that original expressionSets persist
eset1b <- eset1
eset2b <- eset2
eset3b <- eset3

# show number of samples and phenoData columns
dim(Biobase::pData(eset1b))
#> [1] 91 33
dim(Biobase::pData(eset2b))
#> [1] 120  37
dim(Biobase::pData(eset3b))
#> [1] 156  44
# depict phenoData column names
colnames(Biobase::pData(eset1b))
#>  [1] "title"                   "geo_accession"          
#>  [3] "status"                  "submission_date"        
#>  [5] "last_update_date"        "type"                   
#>  [7] "channel_count"           "source_name_ch1"        
#>  [9] "organism_ch1"            "characteristics_ch1"    
#> [11] "characteristics_ch1.1"   "molecule_ch1"           
#> [13] "extract_protocol_ch1"    "label_ch1"              
#> [15] "label_protocol_ch1"      "taxid_ch1"              
#> [17] "hyb_protocol"            "scan_protocol"          
#> [19] "description"             "data_processing"        
#> [21] "platform_id"             "contact_name"           
#> [23] "contact_email"           "contact_institute"      
#> [25] "contact_address"         "contact_city"           
#> [27] "contact_zip/postal_code" "contact_country"        
#> [29] "supplementary_file"      "data_row_count"         
#> [31] "relation"                "sample type:ch1"        
#> [33] "tissue:ch1"
colnames(Biobase::pData(eset2b))
#>  [1] "title"                   "geo_accession"          
#>  [3] "status"                  "submission_date"        
#>  [5] "last_update_date"        "type"                   
#>  [7] "channel_count"           "source_name_ch1"        
#>  [9] "organism_ch1"            "characteristics_ch1"    
#> [11] "characteristics_ch1.1"   "characteristics_ch1.2"  
#> [13] "characteristics_ch1.3"   "molecule_ch1"           
#> [15] "extract_protocol_ch1"    "label_ch1"              
#> [17] "label_protocol_ch1"      "taxid_ch1"              
#> [19] "hyb_protocol"            "scan_protocol"          
#> [21] "description"             "data_processing"        
#> [23] "platform_id"             "contact_name"           
#> [25] "contact_email"           "contact_department"     
#> [27] "contact_institute"       "contact_address"        
#> [29] "contact_city"            "contact_zip/postal_code"
#> [31] "contact_country"         "supplementary_file"     
#> [33] "data_row_count"          "age:ch1"                
#> [35] "gender:ch1"              "stage:ch1"              
#> [37] "tissue:ch1"
colnames(Biobase::pData(eset3b))
#>  [1] "title"                   "geo_accession"          
#>  [3] "status"                  "submission_date"        
#>  [5] "last_update_date"        "type"                   
#>  [7] "channel_count"           "source_name_ch1"        
#>  [9] "organism_ch1"            "characteristics_ch1"    
#> [11] "characteristics_ch1.1"   "characteristics_ch1.2"  
#> [13] "characteristics_ch1.3"   "characteristics_ch1.4"  
#> [15] "molecule_ch1"            "extract_protocol_ch1"   
#> [17] "label_ch1"               "label_protocol_ch1"     
#> [19] "taxid_ch1"               "hyb_protocol"           
#> [21] "scan_protocol"           "description"            
#> [23] "data_processing"         "platform_id"            
#> [25] "contact_name"            "contact_email"          
#> [27] "contact_phone"           "contact_fax"            
#> [29] "contact_department"      "contact_institute"      
#> [31] "contact_address"         "contact_city"           
#> [33] "contact_zip/postal_code" "contact_country"        
#> [35] "contact_web_link"        "supplementary_file"     
#> [37] "data_row_count"          "relation"               
#> [39] "relation.1"              "cell type:ch1"          
#> [41] "gender:ch1"              "overall survival:ch1"   
#> [43] "status:ch1"              "tissue type:ch1"


# remove unnecessary phenoData
# rename desired column names and if necessary the column entries (as we did for 'Tumor')

# eset1b
eset1b$relation <- NULL
eset1b$characteristics_ch1 <- NULL
eset1b$`sample type:ch1` <- NULL
eset1b$`tissue:ch1` <- NULL
colnames(Biobase::pData(eset1b))[8] = targetcol
levels(eset1b[[targetcol]]) = c(controlname, targetname)
dim(eset1b)
#> Features  Samples 
#>    54675       91


# eset2b
eset2b$characteristics_ch1.2 <- NULL
eset2b$characteristics_ch1.3 <- NULL
eset2b$contact_department <- NULL
eset2b$characteristics_ch1 <- NULL
eset2b$`age:ch1` <- NULL
eset2b$`gender:ch1` <- NULL
eset2b$`stage:ch1` <- NULL
eset2b$`tissue:ch1` <- NULL
colnames(Biobase::pData(eset2b))[8] = targetcol
levels(eset2b[[targetcol]]) = c(controlname, targetname)
dim(eset2b)
#> Features  Samples 
#>    54675      120


# eset3b
eset3b$characteristics_ch1.2 <- NULL
eset3b$characteristics_ch1.3 <- NULL
eset3b$characteristics_ch1.4 <- NULL
eset3b$contact_phone <- NULL
eset3b$contact_fax <- NULL
eset3b$contact_department <- NULL
eset3b$contact_web_link <- NULL
eset3b$relation <- NULL
eset3b$source_name_ch1 <- NULL
colnames(Biobase::pData(eset3b))[9] = targetcol
levels(eset3b[[targetcol]]) = c(controlname, targetname)
pNew <- Biobase::pData(eset3b)[,c(1:7,9,8,10:29)] # reorder columns
Biobase::pData(eset3b) = pNew
dim(eset3b)
#> Features  Samples 
#>    54675      156

# now, all phenoData columns are consistent
cbind(colnames(Biobase::pData(eset1b)), colnames(Biobase::pData(eset2b)), 
      colnames(Biobase::pData(eset3b)))
#>       [,1]                      [,2]                     
#>  [1,] "title"                   "title"                  
#>  [2,] "geo_accession"           "geo_accession"          
#>  [3,] "status"                  "status"                 
#>  [4,] "submission_date"         "submission_date"        
#>  [5,] "last_update_date"        "last_update_date"       
#>  [6,] "type"                    "type"                   
#>  [7,] "channel_count"           "channel_count"          
#>  [8,] "target"                  "target"                 
#>  [9,] "organism_ch1"            "organism_ch1"           
#> [10,] "characteristics_ch1.1"   "characteristics_ch1.1"  
#> [11,] "molecule_ch1"            "molecule_ch1"           
#> [12,] "extract_protocol_ch1"    "extract_protocol_ch1"   
#> [13,] "label_ch1"               "label_ch1"              
#> [14,] "label_protocol_ch1"      "label_protocol_ch1"     
#> [15,] "taxid_ch1"               "taxid_ch1"              
#> [16,] "hyb_protocol"            "hyb_protocol"           
#> [17,] "scan_protocol"           "scan_protocol"          
#> [18,] "description"             "description"            
#> [19,] "data_processing"         "data_processing"        
#> [20,] "platform_id"             "platform_id"            
#> [21,] "contact_name"            "contact_name"           
#> [22,] "contact_email"           "contact_email"          
#> [23,] "contact_institute"       "contact_institute"      
#> [24,] "contact_address"         "contact_address"        
#> [25,] "contact_city"            "contact_city"           
#> [26,] "contact_zip/postal_code" "contact_zip/postal_code"
#> [27,] "contact_country"         "contact_country"        
#> [28,] "supplementary_file"      "supplementary_file"     
#> [29,] "data_row_count"          "data_row_count"         
#>       [,3]                     
#>  [1,] "title"                  
#>  [2,] "geo_accession"          
#>  [3,] "status"                 
#>  [4,] "submission_date"        
#>  [5,] "last_update_date"       
#>  [6,] "type"                   
#>  [7,] "channel_count"          
#>  [8,] "target"                 
#>  [9,] "organism_ch1"           
#> [10,] "characteristics_ch1.1"  
#> [11,] "molecule_ch1"           
#> [12,] "extract_protocol_ch1"   
#> [13,] "label_ch1"              
#> [14,] "label_protocol_ch1"     
#> [15,] "taxid_ch1"              
#> [16,] "hyb_protocol"           
#> [17,] "scan_protocol"          
#> [18,] "description"            
#> [19,] "data_processing"        
#> [20,] "platform_id"            
#> [21,] "contact_name"           
#> [22,] "contact_email"          
#> [23,] "contact_institute"      
#> [24,] "contact_address"        
#> [25,] "contact_city"           
#> [26,] "contact_zip/postal_code"
#> [27,] "contact_country"        
#> [28,] "supplementary_file"     
#> [29,] "data_row_count"
```


# Preparations for utilizing the sigident package 

These variables need to be defined for the package functions to work. One could use them also as arguments directly in the respective function. However, we think it is more clearly to define them here at the beginning and to refer at each function to the respective variable.  



```r
plotdir <- "./plots/"
dir.create(plotdir)
#> Warning in dir.create(plotdir): './plots' already exists
csvdir <- "./csv/"
dir.create(csvdir)
#> Warning in dir.create(csvdir): './csv' already exists
species <- "Hs"
OrgDb <- "org.Hs.eg.db"
organism <- "hsa"
pathwayid <- "hsa04110"
seed <- 111
split <- 0.8
```


## Extract Metadata 

TODO Description here.

### Extract sample metadata 


```r
# workaround
GSE18842_meta <- data.table::data.table(
  cbind(
    study = "GSE18842",
    sample = eset1b@phenoData@data[,"geo_accession"],
    target = as.character(eset1b@phenoData@data[,targetcol])
  )
)

GSE19804_meta <- data.table::data.table(
  cbind(
    study = "GSE19804",
    sample = eset2b@phenoData@data[,"geo_accession"],
    target = as.character(eset2b@phenoData@data[,targetcol])
  )
)

GSE19188_meta <- data.table::data.table(
  cbind(
    study = "GSE19188",
    sample = eset3b@phenoData@data[,"geo_accession"],
    target = as.character(eset3b@phenoData@data[,targetcol])
  )
)

sampleMetadata <- rbind(GSE18842_meta,
                        GSE19804_meta,
                        GSE19188_meta)
```

### Define study metadata 


```r
studyMetadata <- data.frame(
  cbind(study = c("GSE18842", "GSE19804", "GSE19188", "GSE30219"),
        discovery = as.logical(c(1, 1, 1, 0)),
        validation = as.logical(c(0, 0, 0, 1))
  ),
  stringsAsFactors = FALSE
)
```

## Merging and Batch Correction

TODO description here

First, a boxplot is created with the included samples on the x-axis and the standardised expression values on the y-axis. `mergedset` results as output from above described merging and represents a matrix containing the genes (Entrez ID) as rows and the samples as columns.


```r
#-----------------------------------------#
##### 4. Merging and Batch Correction #####
#-----------------------------------------#


esets=c(eset1b,eset2b,eset3b)
mergedset <- sigident::merge_(esets)
# merging 3 esets, resulting in a mergeset containing 367 samples and 54675 genes
```



```r
# visualize log2 transformed expression values of the merged data set
filename <- paste0(plotdir, "boxplot_merged_data.jpeg")
jpeg(filename, width = 2000, height = 1000, res = 150, units = "px")
boxplot(mergedset@assayData$exprs, main = "Merged data before batch correction", 
        xlab = "Samples", ylab ="Expression value")
dev.off()
#> png 
#>   2
```


```r
knitr::include_graphics(filename)
```

<img src="./plots/boxplot_merged_data.jpeg" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" width="80%" />

### Discovering Batch Effects 


```r
# generate vector with column diagnosis
table(mergedset@phenoData@data[[targetcol]])
#> 
#>     Control Lung Cancer 
#>         170         197

dd <- sigident::createDiagnosisDesign_(sampleMetadata = sampleMetadata,
                                       studyMetadata = studyMetadata,
                                       controlname = controlname,
                                       targetname = targetname,
                                       targetcol = targetcol)
diagnosis <- dd$diagnosis
length(diagnosis)
#> [1] 367
table(diagnosis)
#> diagnosis
#>   0   1 
#> 170 197

design <- dd$design

# determine batches with number of samples
batch <- sigident::createBatch_(studyMetadata = studyMetadata,
                                sampleMetadata = sampleMetadata)
# check batch assignment
table(batch)
#> batch
#>   1   2   3 
#>  91 120 156
length(batch)
#> [1] 367

# Batch Effect Correction 
# removing batch effects computing linear models, take variance between the diagnosis into consideration
mergeset <- sigident::createCombat_(mergedset = mergedset, batch = batch, design = design)
#> Standardizing Data across genes
dim(mergeset)
#> [1] 21879   367
```

`mergeset` results as output from above described merging and represents a matrix containing the genes (Entrez ID) as rows and the samples as columns.


```r
DF <- mergedset@assayData$exprs
gPCA_before <- sigident::batchDetection_(mergeset = DF,
                                          batch = batch)
filename <- paste0(plotdir, "PCplot_before.png")
sigident::createBatchPlot_(correction_obj = gPCA_before,
                           filename = filename,
                           time = "before")
```


```r
knitr::include_graphics(filename)
```

<img src="./plots/PCplot_before.png" title="plot of chunk unnamed-chunk-15" alt="plot of chunk unnamed-chunk-15" width="80%" />


```r
filename <- paste0(plotdir, "import_boxplot.png")
sigident::createImportBoxplot_(mergeset, filename)
```


```r
knitr::include_graphics(filename)
```

<img src="./plots/import_boxplot.png" title="plot of chunk unnamed-chunk-17" alt="plot of chunk unnamed-chunk-17" width="80%" />


```r
# remove unneeded objects
rm("DF", "GSE18842_meta", "GSE19188_meta", "GSE19804_meta")
```


### Visualize batch effects  

TODO Description here.


```r
gPCA_after <- sigident::batchDetection_(mergeset = mergeset,
                                        batch = batch)
filename <- paste0(plotdir, "PCplot_after.png")
sigident::createBatchPlot_(correction_obj = gPCA_after,
                           filename = filename,
                           time = "after")
```


```r
knitr::include_graphics(filename)
```

<img src="./plots/PCplot_after.png" title="plot of chunk unnamed-chunk-20" alt="plot of chunk unnamed-chunk-20" width="80%" />

# DEG Analysis 

TODO Description here.


```r
q_selection <- NULL # use default setting
deg_q <- sigident::qSelection_(sampleMetadata = sampleMetadata, 
                               studyMetadata = studyMetadata,
                               deg.q.selection = q_selection)

genes <- sigident::identifyDEGs_(mergeset = mergeset,
                                 design = design,
                                 qValue = deg_q)

# heatmap creation
filename <- paste0(plotdir, "DEG_heatmap.png")
# create colors for map
ht_colors <- sigident::colorHeatmap_(sampleMetadata = sampleMetadata,
                                     studyMetadata = studyMetadata,
                                     targetcol = targetcol,
                                     controlname = controlname) # cancer = red
sigident::createDEGheatmap_(mergeset = mergeset,
                            genes = genes,
                            patientcolors = ht_colors,
                            filename = filename)
```


```r
knitr::include_graphics(filename)
```

<img src="./plots/DEG_heatmap.png" title="plot of chunk unnamed-chunk-22" alt="plot of chunk unnamed-chunk-22" width="80%" />

TODO description here


```r
deg_info <- sigident::exportDEGannotations_(mergedset = mergedset)
data.table::fwrite(deg_info, paste0(csvdir, "DEG_info.csv"))
```


```r
dim(deg_info)
#> [1] 54675     5
knitr::kable(head(deg_info))
```



|probe_ID  |gene_symbol      |gene_title                                                    |genebank_accession |entrez_id          |
|:---------|:----------------|:-------------------------------------------------------------|:------------------|:------------------|
|1007_s_at |DDR1 /// MIR4640 |discoidin domain receptor tyrosine kinase 1 /// microRNA 4640 |U48705             |780 /// 100616237  |
|1053_at   |RFC2             |replication factor C (activator 1) 2, 40kDa                   |M87338             |5982               |
|117_at    |HSPA6            |heat shock 70kDa protein 6 (HSP70B')                          |X51757             |3310               |
|121_at    |PAX8             |paired box 8                                                  |X69699             |7849               |
|1255_g_at |GUCA1A           |guanylate cyclase activator 1A (retina)                       |L36861             |2978               |
|1294_at   |MIR5193 /// UBA7 |microRNA 5193 /// ubiquitin-like modifier activating enzyme 7 |L13852             |7318 /// 100847079 |


```r
deg_results <- sigident::limmaTopTable_(mergeset = mergeset,
                                        design = design,
                                        qValue = deg_q)
data.table::fwrite(deg_results, paste0(csvdir, "DEG_results.csv"))
```


```r
dim(deg_results)
#> [1] 370   3
knitr::kable(head(deg_results))
```



|                   |Probe ID           |     logFC| adj.qValue|
|:------------------|:------------------|---------:|----------:|
|177                |177                | -5.357857|          0|
|2823               |2823               | -3.946057|          0|
|7134               |7134               | -4.090512|          0|
|11170              |11170              | -3.420523|          0|
|2277 /// 100532742 |2277 /// 100532742 | -4.132929|          0|
|100505495          |100505495          | -2.455868|          0|

# Gene enrichment 

TODO Description here.


```r
deg_entrez <- unique(genes)
```


## Test for over-representation 

TODO Description here.


```r
enr_topgo <- sigident::extractGOterms_(entrez = deg_entrez,
                                       species = species)
```

```r
dim(enr_topgo)
#> [1] 20  5
knitr::kable(head(enr_topgo))
```



|           |Term                                   |Ont |    N|  DE| P.DE|
|:----------|:--------------------------------------|:---|----:|---:|----:|
|GO:0005576 |extracellular region                   |CC  | 4295| 141|    0|
|GO:0050896 |response to stimulus                   |BP  | 8836| 219|    0|
|GO:0140014 |mitotic nuclear division               |BP  |  237|  28|    0|
|GO:0070887 |cellular response to chemical stimulus |BP  | 3101| 110|    0|
|GO:0000819 |sister chromatid segregation           |BP  |  157|  23|    0|
|GO:0008283 |cell proliferation                     |BP  | 1969|  82|    0|

TODO Description here.


```r
enr_topkegg <- sigident::extractKEGGterms_(entrez = deg_entrez,
                                           species = species)
```

```r
dim(enr_topkegg)
#> [1] 20  4
knitr::kable(head(enr_topkegg))
```



|              |Pathway                                                       |   N| DE|      P.DE|
|:-------------|:-------------------------------------------------------------|---:|--:|---------:|
|path:hsa04110 |Cell cycle                                                    | 124| 13| 0.0000009|
|path:hsa05144 |Malaria                                                       |  49|  8| 0.0000043|
|path:hsa04512 |ECM-receptor interaction                                      |  88|  8| 0.0003218|
|path:hsa05418 |Fluid shear stress and atherosclerosis                        | 139| 10| 0.0003999|
|path:hsa04061 |Viral protein interaction with cytokine and cytokine receptor | 100|  8| 0.0007634|
|path:hsa04610 |Complement and coagulation cascades                           |  79|  7| 0.0008918|

TODO Description here.


```r
enr_fitlm <- sigident::goDiffReg_(mergeset = mergeset,
                                  design = design)

enr_fitlm_topgo <- sigident::extractGOterms_(entrez = enr_fitlm,
                                             species = species,
                                             FDR = 0.01)
data.table::fwrite(enr_fitlm_topgo, paste0(csvdir, "Top_GO_fitlm.csv"))
```

```r
dim(enr_fitlm_topgo)
#> [1] 20  7
knitr::kable(head(enr_fitlm_topgo))
```



|           |Term                                     |Ont |     N|   Up| Down| P.Up|    P.Down|
|:----------|:----------------------------------------|:---|-----:|----:|----:|----:|---------:|
|GO:0070013 |intracellular organelle lumen            |CC  |  4844| 2168| 1015|    0| 0.6619777|
|GO:0031974 |membrane-enclosed lumen                  |CC  |  4844| 2168| 1015|    0| 0.6619777|
|GO:0043233 |organelle lumen                          |CC  |  4844| 2168| 1015|    0| 0.6619777|
|GO:0043231 |intracellular membrane-bounded organelle |CC  |  9294| 3683| 2057|    0| 0.0008912|
|GO:0044446 |intracellular organelle part             |CC  |  8295| 3318| 1881|    0| 0.0000055|
|GO:0005622 |intracellular                            |CC  | 12445| 4626| 2931|    0| 0.0000000|


```r
enr_fitlm_topkegg <- sigident::extractKEGGterms_(entrez = enr_fitlm,
                                                 species = species)
data.table::fwrite(enr_fitlm_topkegg, paste0(csvdir, "Top_KEGG_fitlm.csv"))
```

```r
dim(enr_fitlm_topkegg)
#> [1] 20  6
knitr::kable(head(enr_fitlm_topkegg))
```



|              |Pathway                     |   N| Up| Down|      P.Up|    P.Down|
|:-------------|:---------------------------|---:|--:|----:|---------:|---------:|
|path:hsa03013 |RNA transport               | 143| 99|   16| 0.0000000| 0.9999231|
|path:hsa05332 |Graft-versus-host disease   |  36|  0|   29| 1.0000000| 0.0000000|
|path:hsa04380 |Osteoclast differentiation  | 126| 20|   65| 0.9999999| 0.0000000|
|path:hsa04062 |Chemokine signaling pathway | 181| 43|   84| 0.9999275| 0.0000000|
|path:hsa04640 |Hematopoietic cell lineage  |  92| 12|   52| 0.9999999| 0.0000000|
|path:hsa05133 |Pertussis                   |  69| 11|   42| 0.9999586| 0.0000000|

## GO Analysis

TODO Description here.


```r
enr_analysis <- sigident::goEnrichmentAnalysis_(entrez = deg_entrez,
                                                OrgDB = OrgDb,
                                                organism = organism,
                                                fitlm = enr_fitlm,
                                                pathwayid = pathwayid,
                                                species = organism,
                                                plotdir = plotdir)
```

```r
filename <- paste0(plotdir, "/", organism, "04110.png")
```

```r
knitr::include_graphics(filename)
```

<img src="./plots//hsa04110.png" title="plot of chunk unnamed-chunk-38" alt="plot of chunk unnamed-chunk-38" width="80%" />

```r
filename <- paste0(plotdir, "/", organism, "04110.pathview.png")
```

```r
knitr::include_graphics(filename)
```

<img src="./plots//hsa04110.pathview.png" title="plot of chunk unnamed-chunk-40" alt="plot of chunk unnamed-chunk-40" width="80%" />

TODO Description here.


```r
filename <- paste0(plotdir, "Enriched_GO.png")
sigident::createEnrichtedBarplot_(enrichmentobj = enr_analysis$go,
                                  type = "GO",
                                  filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/Enriched_GO.png" title="plot of chunk unnamed-chunk-42" alt="plot of chunk unnamed-chunk-42" width="80%" />

TODO Description here.


```r
filename <- paste0(plotdir, "Enriched_KEGG.png")
sigident::createEnrichtedBarplot_(enrichmentobj = enr_analysis$kegg,
                                  type = "KEGG",
                                  filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/Enriched_KEGG.png" title="plot of chunk unnamed-chunk-44" alt="plot of chunk unnamed-chunk-44" width="80%" />


# Identification of diagnostic signatures 

TODO Description here.

## Split data into training and test dataset

TODO Description here.


```r
training_list <- sigident::createTrainingTest_(diagnosis = diagnosis,
                                               mergeset = mergeset,
                                               split = split,
                                               seed = seed)
```

## Lasso regression

TODO Description here.


```r
diagnostic_lasso <- sigident::signature_(traininglist = training_list,
                                         type = "lasso",
                                         nfolds = 10,
                                         seed = seed)
filename <- paste0(plotdir, "CV_lasso.png")
sigident::createCVPlot_(cv_obj = diagnostic_lasso$fitCV,
                        filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/CV_lasso.png" title="plot of chunk unnamed-chunk-47" alt="plot of chunk unnamed-chunk-47" width="80%" />


```r
filename <- paste0(plotdir, "ROC_Lasso.min.png")
sigident::createROCplot_(roc = diagnostic_lasso$roc.min,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_Lasso.min.png" title="plot of chunk unnamed-chunk-49" alt="plot of chunk unnamed-chunk-49" width="80%" />


```r
filename <- paste0(plotdir, "ROC_Lasso.1se.png")
sigident::createROCplot_(roc = diagnostic_lasso$roc.1se,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_Lasso.1se.png" title="plot of chunk unnamed-chunk-51" alt="plot of chunk unnamed-chunk-51" width="80%" />


## Elastic net regression 

TODO Description here.


```r
diagnostic_elasticnet <- sigident::signature_(traininglist = training_list,
                                              type = "elastic",
                                              alpha = 0.9,
                                              nfolds = 10,
                                              seed = seed)
filename <- paste0(plotdir, "CV_elasticNet.png")
sigident::createCVPlot_(cv_obj = diagnostic_elasticnet$fitCV,
                        filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/CV_elasticNet.png" title="plot of chunk unnamed-chunk-53" alt="plot of chunk unnamed-chunk-53" width="80%" />


```r
filename <- paste0(plotdir, "ROC_elasticNet.min.png")
sigident::createROCplot_(roc = diagnostic_elasticnet$roc.min,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_elasticNet.min.png" title="plot of chunk unnamed-chunk-55" alt="plot of chunk unnamed-chunk-55" width="80%" />


```r
filename <- paste0(plotdir, "ROC_elasticNet.1se.png")
sigident::createROCplot_(roc = diagnostic_elasticnet$roc.1se,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_elasticNet.1se.png" title="plot of chunk unnamed-chunk-57" alt="plot of chunk unnamed-chunk-57" width="80%" />


## Gridsearch to find alpha and lambda

TODO Description here.


```r
diagnostic_glmGrid <- sigident::signature_(traininglist = training_list,
                                           type = "grid",
                                           nfolds = 10,
                                           seed = seed)
# plot model of gridsearch
filename <- paste0(plotdir, "Gridsearch_model.png")
sigident::createGridModelPlot_(model = diagnostic_glmGrid$caret.train,
                               filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/Gridsearch_model.png" title="plot of chunk unnamed-chunk-59" alt="plot of chunk unnamed-chunk-59" width="80%" />


```r
# plot variable importance of gridsearch
filename <- paste0(plotdir, "Gridsearch_variable_importance.png")
sigident::createGridVarImpPlot_(model = diagnostic_glmGrid$caret.train,
                                filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/Gridsearch_variable_importance.png" title="plot of chunk unnamed-chunk-61" alt="plot of chunk unnamed-chunk-61" width="80%" />


```r
# create roc plot
filename <- paste0(plotdir, "ROC_elasticNet.grid.png")
sigident::createROCplot_(roc = diagnostic_glmGrid$roc.elasticNet,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_elasticNet.grid.png" title="plot of chunk unnamed-chunk-63" alt="plot of chunk unnamed-chunk-63" width="80%" />

# Analysis of Diagnostic Models 

## Comparing Model Performance  

These models can be included into one big list and the results can then be printed using the function `compareDiagnosticModels()`.


```r
diagnosticModels <- list(
  lasso.min = list(model = diagnostic_lasso$lambda.min,
                   confmat = diagnostic_lasso$confmat.min,
                   auc = as.numeric(diagnostic_lasso$roc.min$auc)),
  lasso.1se = list(model = diagnostic_lasso$lambda.1se,
                   confmat = diagnostic_lasso$confmat.1se,
                   auc = as.numeric(diagnostic_lasso$roc.1se$auc)),
  elastic.min = list(model = diagnostic_elasticnet$lambda.min,
                     confmat = diagnostic_elasticnet$confmat.min,
                     auc = as.numeric(diagnostic_elasticnet$roc.min$auc)),
  elastic.1se = list(model = diagnostic_elasticnet$lambda.1se,
                     confmat = diagnostic_elasticnet$confmat.1se,
                     auc = as.numeric(diagnostic_elasticnet$roc.1se$auc)),
  elastic.grid = list(model = diagnostic_glmGrid$elasticNet.auto,
                      confmat = diagnostic_glmGrid$confmat.elasticNet,
                      auc = as.numeric(diagnostic_glmGrid$roc.elasticNet$auc))
  )
```


```r
knitr::kable(
  sigident::compareDiagnosticModels(diagnosticModels)
)
```



|Model        |   AUC|
|:------------|-----:|
|lasso.min    | 0.983|
|lasso.1se    | 0.983|
|elastic.min  | 0.982|
|elastic.1se  | 0.983|
|elastic.grid | 0.975|

## Export selected Entrez-IDs


```r
# map selected variables and export 
signature_genes <- sigident::geneMapSig_(mergeset = mergeset, model = diagnostic_lasso$lambda.min)
data.table::fwrite(signature_genes, paste0(csvdir, "signature_genes_model1.csv"))
```


```r
dim(signature_genes)
#> [1] 53  1
head(signature_genes)
#>            Entrez_ID
#> 1             123591
#> 2             151306
#> 3             158158
#> 4             284656
#> 5 8974 /// 101927705
#> 6              64399
```


# Identification of Prognostic Signature

## Create a list that contains specifications of the study/studies that contain(s) survival time information 


```r
discoverystudies.w.timedata <- list("GSE19188" = list(timecol = "characteristics_ch1.2",
                                                      statuscol = "characteristics_ch1.3",
                                                      targetcolname = "characteristics_ch1",
                                                      controllevelname = "tissue type: healthy",
                                                      targetlevelname = "tissue type: tumor"))
```

## Identification of survival correlated genes 

### Extract survival data 


```r
#---------------------------------------------------#
##### 8. Identification of Prognostic Signature #####
#---------------------------------------------------#

### 8.A Identification of survival correlated genes
survivalData <- sigident::getSurvivalTime_(studyMetadata = studyMetadata,
                                           sampleMetadata = sampleMetadata,
                                           discoverystudies.w.timedata = discoverystudies.w.timedata,
                                           targetname = targetname,
                                           controlname = controlname,
                                           targetcol = targetcol,
                                           datadir = datadir)
survTable <- survivalData[["GSE19188"]]$survTable
entrezIDs <- survivalData[["GSE19188"]]$entrezIDs
```


```r
dim(survTable)
#> [1] 156  84
head(survTable)
#>           time status 780 /// 100616237 7318 /// 100847079        7067
#> GSM475656 12.5     NA       0.278536376         -0.3665863  0.47354494
#> GSM475657   NA      1      -0.005101812          0.1405779  0.43960675
#> GSM475658   NA      1      -0.619296868          0.3662965 -0.52670499
#> GSM475659   NA      1      -0.400919703          0.6306993 -0.08487694
#> GSM475660   NA      1      -0.025467569          0.4217360  0.09975505
#> GSM475661   NA      1      -1.277200911         -1.5330581  0.07878199
#>                 6352        1548      23170       10406        5594
#> GSM475656  0.8875675  0.38299987  0.2955272  0.91245794 -1.16121413
#> GSM475657  0.9791955  0.16722715 -0.3997257 -0.19852185 -0.44646059
#> GSM475658  0.2855274 -0.19247590 -0.4227795 -0.22619978 -0.07487306
#> GSM475659  0.7294725 -0.18182166 -0.5488857 -0.06133111  0.10894611
#> GSM475660 -0.1806954  0.05652838 -0.2739466  0.13451377 -0.22083194
#> GSM475661 -4.4714688 -0.20820236  1.1235085 -0.22022012 -0.36546548
#>                128153      163154      54899       57617      113235
#> GSM475656 -0.07582036  0.12431423 -1.1389247  0.53926500  0.07835584
#> GSM475657  0.73072404 -0.27225240 -0.4493079 -0.17457677 -0.41100526
#> GSM475658 -0.11863885 -0.36833212  0.7679284  0.08809509 -0.42789991
#> GSM475659 -0.28841898 -0.17162174 -0.5296897  0.19254496  0.15441037
#> GSM475660  0.13282074 -0.01726770 -0.2708726  0.44122320 -0.31260598
#> GSM475661  0.26276471  0.03603483  2.2488019 -0.04348902 -0.28182835
#>                91937 79844 /// 653082        172     148113       54965
#> GSM475656 -0.1617181        0.1318657 -0.2801627  0.5408657  0.97308481
#> GSM475657 -0.3877970        2.1742008  0.6111931 -0.8620225 -0.27847977
#> GSM475658 -0.2950219        1.1420088  0.1214691 -0.7928894 -0.64063761
#> GSM475659 -0.2509978        1.8861580  0.6650443 -0.7168694 -0.08739238
#> GSM475660 -0.2523125       -0.2037104 -0.2957881 -0.5780904 -0.55568117
#> GSM475661 -0.1408117       -2.2607172  0.7984810 -0.2761466  0.95754971
#>                 91252      221264      113277      84920      91624
#> GSM475656  0.10406473  0.61309769 -0.27860441  0.1174461 -1.3514100
#> GSM475657 -0.25599605 -0.15101972  0.47419215 -0.3135123  0.9981771
#> GSM475658 -0.27187265 -0.05174585 -0.70553514 -0.5018055  1.1186026
#> GSM475659 -0.58801060 -0.27440977  0.01214489 -0.3266313  1.5886667
#> GSM475660  0.08653026 -0.30907107 -0.31270477 -0.2501872  0.1419729
#> GSM475661  0.00526674 -0.40916447  0.03347595 -0.2346556 -1.9271205
#>                 4238     170575       85478      159091     220136
#> GSM475656 -0.2078267 -0.9112034 -0.15163468  0.12178865 -0.7252915
#> GSM475657 -0.3294961  0.7381955  0.83172874  0.24266971  2.8881303
#> GSM475658 -0.2357807  0.8773266 -0.01261386  0.08480011 -0.5941382
#> GSM475659 -0.4746650  1.2085270 -0.16921275  0.14216124  0.1375772
#> GSM475660 -0.5325492  0.4102301  0.88760538 -0.19504977  2.2262432
#> GSM475661  0.2801247 -1.5761351 -0.17531030 -0.06676253 -0.6752415
#>                  5930      92806       11078        3233        5150
#> GSM475656 -0.16110583  0.2431688 -0.65124772  0.56759661 -0.02645938
#> GSM475657 -0.06302431 -0.2619877  0.82603106 -0.16573512  0.54508755
#> GSM475658  0.29415973 -0.2954114 -0.48354213 -0.12874636 -0.69667993
#> GSM475659  0.29858607  0.3345876  1.03017049 -0.36890757  0.39476103
#> GSM475660 -1.17882497  0.4414852 -0.45822590  0.15854543  0.18927089
#> GSM475661  0.96915173  1.4120668 -0.04659303  0.03343629  0.31892804
#>                29883     260429      255057      114609     124540
#> GSM475656 -0.3841586  0.1625087  0.42836254 -0.35342172 -0.1196380
#> GSM475657 -0.3885528 -0.1149259  0.05422687  0.37183347 -0.9238617
#> GSM475658 -0.1486220 -0.0856261 -0.00042600  0.13887136 -0.9781208
#> GSM475659 -0.2093516 -0.1394833 -0.02706574  0.19682419 -0.8358295
#> GSM475660 -0.9767422  0.4303705  0.36854862 -0.07562991 -0.2147314
#> GSM475661  2.8255018 -0.4330280  0.31337882 -0.21008887  1.5669405
#>                85477      132321       84449       157506       135295
#> GSM475656  2.0469556 -1.28852588 -0.10839580 -0.078252838 -0.137168565
#> GSM475657 -0.4765486  0.45659061  0.53593648  0.040095106 -0.326381726
#> GSM475658 -0.4150521  0.28437154  0.29818371 -0.002294732  0.272287357
#> GSM475659 -0.4384465  0.45064700  0.54803838  0.142776263 -0.008534535
#> GSM475660 -0.6432211  0.04404273  0.02614389 -0.130727883 -0.233841553
#> GSM475661 -1.1557350 -0.41218374  0.18249916 -0.057319801  1.063468311
#>                 202309     203111      150350       81629      140870
#> GSM475656 -0.637652966  0.5243408  0.13300716  0.27115411 -0.09117015
#> GSM475657 -0.111551519  0.3378671 -0.05389970  0.05574549  0.51269486
#> GSM475658  0.009126382 -0.2029729  0.03736802 -0.40694137 -0.32229488
#> GSM475659  0.411190558 -0.3834744 -0.01008170 -0.48586082  0.04661639
#> GSM475660 -0.181443306 -0.2896606  0.18044051  0.23407509 -0.12340252
#> GSM475661 -0.855063105 -0.5154584 -0.16857919  0.20740689  0.13874551
#>               160364         2972      123591 571 /// 100379661
#> GSM475656 -1.3158803  0.224653987  0.04450449        0.10701042
#> GSM475657 -0.2223639 -0.042325894 -0.20482490        0.09389595
#> GSM475658  2.5103626 -0.156351653 -0.18044753       -0.02005741
#> GSM475659  0.6968696 -0.148395252 -0.16658369       -0.24895364
#> GSM475660  1.3723640 -0.007567775 -0.19924807        0.02025438
#> GSM475661 -2.4132640 -0.009225840 -0.15952311        0.05768440
#>                126206      165530 245909 /// 503841      259240     121441
#> GSM475656  0.09760723  0.02604634        0.10865766  0.06613000 -0.4525503
#> GSM475657 -0.13054675 -0.04842479        0.03704578 -0.04696868  0.0806653
#> GSM475658 -0.23724318 -0.35664951       -0.72306703 -0.14210613 -0.2063048
#> GSM475659 -0.23964383  0.09686435        0.07284487 -0.12954904 -0.3802460
#> GSM475660  0.42173484  0.23982987        0.20812146  0.26165222 -0.3935905
#> GSM475661  0.12728337  0.07415260        0.49254374  0.17126916  0.9658228
#>                254173      125972      220979        2117       317719
#> GSM475656  0.28960049 -0.25707171  0.62106444 -0.47101319 -0.014389032
#> GSM475657 -0.02496260 -0.10924410 -0.08845968  0.12716405  0.023269085
#> GSM475658 -0.19333057 -0.18773429  0.68069336  0.42587910 -0.005139628
#> GSM475659 -0.02966311 -0.08009281  0.52800933 -0.64405488 -0.190897384
#> GSM475660  0.64369195 -0.04114249  0.10938982 -0.27645374 -0.106289287
#> GSM475661  0.16597716 -0.15212603 -0.15332420 -0.07025383  0.086453224
#>                220992      162387 64072 /// 100653137       84465
#> GSM475656  0.38467515 -0.09202981          0.03046847  2.44610616
#> GSM475657 -0.06545603 -0.22127891         -0.06839991  0.46251777
#> GSM475658  0.18012662 -0.81065784         -0.37587378 -0.04861362
#> GSM475659 -0.08222380 -0.16606038         -0.33936931  0.57853291
#> GSM475660 -0.24850125 -0.02731787          0.01067467 -0.28217290
#> GSM475661  1.30315497  0.09767578          0.11140833 -0.18640771
#>                 11318       80712 147199 /// 653486      285126
#> GSM475656 -0.13368814  0.12754646       0.009682731  0.16993225
#> GSM475657 -0.31987640 -0.09137923      -0.286299152 -0.29759720
#> GSM475658  1.06945924 -0.06434080      -0.155012869 -0.07046859
#> GSM475659  0.04024538  0.04079875       0.112325121 -0.07768369
#> GSM475660  0.01116855  0.07382120      -0.146128999  0.36579427
#> GSM475661  0.04528832  0.11639248      -0.100699510  0.04117625
#>                126248      158471       125997        92949        85509
#> GSM475656 -0.06117258  1.55819034 -0.300846986  0.017812505  0.371721193
#> GSM475657  0.16679397 -0.07261160 -0.434666422  0.011226717  0.126828186
#> GSM475658 -0.01910797 -0.11760205 -0.119045969 -0.009253004 -0.424766600
#> GSM475659 -0.02114144  0.21192769 -0.006409731 -0.075286382 -0.121452997
#> GSM475660  0.19095790 -0.04059053 -0.185981481 -0.175078145 -0.003727694
#> GSM475661  0.10019611 -0.51028690  0.190430798  2.177766365  0.178767803
#>                89778      118421      259234        83451      23527
#> GSM475656 -0.3373551 -0.11676088  0.56101992 -0.314147138  0.2128014
#> GSM475657  2.0518806 -0.07544516 -0.08450644  0.136568284 -0.6620996
#> GSM475658 -0.2639671  0.39272939 -0.09101155 -0.290484843 -0.2358734
#> GSM475659 -0.3418193 -0.21847370 -0.47412441 -0.208416239 -0.7324841
#> GSM475660  0.2374486 -0.04384302  0.13974709 -0.004050877 -1.0119304
#> GSM475661 -0.3363428  0.12474152 -0.10206668  0.447004465 -0.8630060
#>                  2593       3664
#> GSM475656 -0.12241182  0.4714302
#> GSM475657 -0.39589762 -0.9062970
#> GSM475658 -0.60780261 -1.3907893
#> GSM475659 -0.53227797 -1.6478505
#> GSM475660 -0.05546723 -0.8218098
#> GSM475661 -0.33322556 -2.4100944
```


```r
print(entrezIDs)
#>  [1] "780 /// 100616237"   "7318 /// 100847079"  "7067"               
#>  [4] "6352"                "1548"                "23170"              
#>  [7] "10406"               "5594"                "128153"             
#> [10] "163154"              "54899"               "57617"              
#> [13] "113235"              "91937"               "79844 /// 653082"   
#> [16] "172"                 "148113"              "54965"              
#> [19] "91252"               "221264"              "113277"             
#> [22] "84920"               "91624"               "4238"               
#> [25] "170575"              "85478"               "159091"             
#> [28] "220136"              "5930"                "92806"              
#> [31] "11078"               "3233"                "5150"               
#> [34] "29883"               "260429"              "255057"             
#> [37] "114609"              "124540"              "85477"              
#> [40] "132321"              "84449"               "157506"             
#> [43] "135295"              "202309"              "203111"             
#> [46] "150350"              "81629"               "140870"             
#> [49] "160364"              "2972"                "123591"             
#> [52] "571 /// 100379661"   "126206"              "165530"             
#> [55] "245909 /// 503841"   "259240"              "121441"             
#> [58] "254173"              "125972"              "220979"             
#> [61] "2117"                "317719"              "220992"             
#> [64] "162387"              "64072 /// 100653137" "84465"              
#> [67] "11318"               "80712"               "147199 /// 653486"  
#> [70] "285126"              "126248"              "158471"             
#> [73] "125997"              "92949"               "85509"              
#> [76] "89778"               "118421"              "259234"             
#> [79] "83451"               "23527"               "2593"               
#> [82] "3664"
```

### Compute univariate cox regression and determine significance of each gene through separate univariate Cox regressions 


```r
surv_correlated <- sigident::univCox_(survTable = survTable,
                                      entrezIDs = entrezIDs)
# export table with survival correlated genes
data.table::fwrite(surv_correlated, paste0(csvdir, "survival_correlated_genes.csv"))
```


```r
dim(surv_correlated)
#> [1] 8 5
head(surv_correlated)
#>        Entrez_ID  beta HR (95% CI for HR) wald.test p.value
#> 7067        7067     2       7.4 (1.6-34)       6.5  0.0110
#> 6352        6352 -0.31   0.73 (0.57-0.94)         6  0.0140
#> 221264    221264   1.6       4.9 (1.6-15)       7.6  0.0058
#> 113277    113277  -1.1   0.34 (0.16-0.73)       7.6  0.0060
#> 29883      29883  0.74      2.1 (1.2-3.8)       6.2  0.0130
#> 255057    255057   3.8       45 (4.6-440)        11  0.0011
```


## Prognostic classifier and Kaplan-Meier estimator

classification of patients into high and low risk groups regarding to expression profile of signature genes

due to selection bias, only GSE18842 and GSE19804 are used for calculating classifier


```r
classifier_studies <- c("GSE18842", "GSE19804")
```


### Evaluate expression pattern (over-/underrepresentation)

generate vectors 'control' and 'tumor' for function expressionPattern_


```r
exprPattern <- sigident::generateExpressionPattern_(classifier_studies = classifier_studies,
                                                    sigCov = surv_correlated,
                                                    mergeset = mergeset,
                                                    studyMetadata = studyMetadata,
                                                    sampleMetadata = sampleMetadata,
                                                    controlname = controlname,
                                                    targetname = targetname,
                                                    targetcol = targetcol)
```


```r
dim(exprPattern)
#> [1] 8 2
head(exprPattern)
#>     Gene Over/Under
#> 1   7067      Under
#> 2   6352      Under
#> 3 221264       Over
#> 4 113277       Over
#> 5  29883      Under
#> 6 255057       Over
```

### Apply the prognostic classifier on validation data set

#### Create a list that contains specifications of the study that contains the validation information


```r
validationstudiesinfo <- list("GSE30219" = list(timecol = "characteristics_ch1.9",
                                                statuscol = "characteristics_ch1.8",
                                                targetcolname = "source_name_ch1",
                                                controllevelname = "Non Tumoral Lung",
                                                targetlevelname = "Lung Tumour"))
```


```r
pC <- sigident::prognosticClassifier_(PatternCom = exprPattern,
                                      validationstudiesinfo = validationstudiesinfo,
                                      datadir = datadir,
                                      controlname = controlname,
                                      targetname = targetname,
                                      targetcol = targetcol)
fit <- pC[["GSE30219"]]$kaplan.estimator$fit
RiskTable <- pC[["GSE30219"]]$risktable
```


```r
dim(RiskTable)
#> [1] 307   3
head(RiskTable)
#>           time status Groups
#> GSM748053  121      1      1
#> GSM748054   21      2      0
#> GSM748055   86      2      1
#> GSM748056   10      2      1
#> GSM748057   81      2      0
#> GSM748058   68      1      1
```

#### View proportional hazards regression models


```r
pC[["GSE30219"]]$kaplan.estimator$res.cox
#> Call:
#> survival::coxph(formula = survival::Surv(time, status) ~ Groups, 
#>     data = RiskTable)
#> 
#>           coef exp(coef) se(coef)      z    p
#> Groups -0.1243    0.8831   0.1574 -0.789 0.43
#> 
#> Likelihood ratio test=0.63  on 1 df, p=0.4264
#> n= 278, number of events= 188 
#>    (29 observations deleted due to missingness)
```


```r
fit
#> Call: survfit(formula = res.cox, newdata = new_df)
#> 
#>     n events median 0.95LCL 0.95UCL
#> 1 278    188     48      28      70
#> 2 278    188     57      33     103
```


#### Plot prognostic Kaplan-Meier plot


```r
# create roc plot
filename <- paste0(plotdir, "Prognostic_Kaplan-Meier_Plot.png")
sigident::createSurvPlot_(fit = fit,
                          RiskTable = RiskTable,
                          filename = filename)
```


```r
knitr::include_graphics(filename)
```

<img src="./plots/Prognostic_Kaplan-Meier_Plot.png" title="plot of chunk unnamed-chunk-83" alt="plot of chunk unnamed-chunk-83" width="80%" />

# An alternative workflow: two generic functions for signature identification 

In order to simplify the whole abovedescribed workflow, we wrapped all these functions into one big function.

## Preprocessing: 

The preprocessing needs to be conform with the above described approach.

## Run `sigidentDiagnostic`-Function 



```r
results <- sigident::sigidentDiagnostic(mergeset = mergeset,
                                        controlname = "Control",
                                        targetname = "Lung Cancer",
                                        studyMetadata = studyMetadata,
                                        sampleMetadata = sampleMetadata,
                                        species = "Hs",
                                        OrgDB = "org.Hs.eg.db",
                                        organism = "hsa",
                                        pathwayid = "hsa04110",
                                        deg.q.selection = NULL,
                                        seed = 111,
                                        nfolds = 10,
                                        split = 0.8,
                                        plotdir = "./plots/",
                                        csvdir = "./tables/",
                                        targetcol = "target")
```


```r
diagnosticModels <- results$diagnosticModels
sigUtils <- results$utils
```


```r
knitr::kable(
  sigident::compareDiagnosticModels(diagnosticModels)
)
```

## Run `sigidentPrognostic`-Function 


```r
progn_results <- sigident::sigidentPrognostic(mergeset = mergeset,
                                              studyMetadata = studyMetadata,
                                              sampleMetadata = sampleMetadata,
                                              discoverystudies.w.timedata = discoverystudies.w.timedata,
                                              classifier_studies = classifier_studies,
                                              validationstudiesinfo = validationstudiesinfo,
                                              controlname = "Control",
                                              targetname = "Lung Cancer",
                                              datadir = datadir,
                                              plotdir = "./plots/",
                                              csvdir = "./tables/",
                                              targetcol = "target")
```


```r
# create roc plot
filename <- paste0(plotdir, "Prognostic_Kaplan-Meier_Plot.png")
sigident::createSurvPlot_(fit = progn_results$GSE19188$GSE30219$fit,
                          RiskTable = progn_results$GSE19188$GSE30219$risktable,
                          filename = filename)
```


```r
knitr::include_graphics(filename)
```

