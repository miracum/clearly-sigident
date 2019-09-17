---
title: "sigident - Howto Microarray"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{sigident-Howto_Microarray}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
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

TODO description


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

#------------------------------------------------------#
###### 1. Download of GEO datasets & Normalization #####
#------------------------------------------------------#

# insert the desired GEO dataset into the getGEO() function
GSE18842 <- GEOquery::getGEO("GSE18842", destdir = datadir)
GSE19804 <- GEOquery::getGEO("GSE19804", destdir = datadir)

# extract expressionSets
eset1 <- GSE18842[[1]]
eset2 <- GSE19804[[1]]

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

TODO description here.


```r
#-------------------------------------#
##### 2. Formatting the phenoData #####
#-------------------------------------#

# the phenoData of the expresssionSets must be formatted in such a way as that the phenoData columns overlap perfectly
# all expressionSets should have the same columns in the same order
# not consistent columns must get removed
# for simpler handling only phenoData columns of interest should be retained
# we here retained all consistent phenoData

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
colnames(pData(eset2b))[8] = targetcol
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
cbind(colnames(Biobase::pData(eset1b)), colnames(Biobase::pData(eset2b)), colnames(Biobase::pData(eset3b)))
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

## Checking normalization


```r
#-----------------------------------#
##### 3. Checking normalization #####
#-----------------------------------#

graphics::boxplot(Biobase::exprs(eset1b))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

```r
graphics::boxplot(Biobase::exprs(eset2b))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png)

```r
graphics::boxplot(Biobase::exprs(eset3b))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-3.png)

```r

#overview of how many tumor and cancer samples are in esets
table(Biobase::pData(eset1b)[,targetcol])
#> 
#>     Control Lung Cancer 
#>          45          46
table(Biobase::pData(eset2b)[,targetcol])
#> 
#>     Control Lung Cancer 
#>          60          60
table(Biobase::pData(eset3b)[,targetcol])
#> 
#>     Control Lung Cancer 
#>          65          91
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
  )
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

<img src="./plots/boxplot_merged_data.jpeg" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" width="80%" />

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

<img src="./plots/PCplot_before.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" width="80%" />


```r
filename <- paste0(plotdir, "import_boxplot.png")
sigident::createImportBoxplot_(mergeset, filename)
```


```r
knitr::include_graphics(filename)
```

<img src="./plots/import_boxplot.png" title="plot of chunk unnamed-chunk-15" alt="plot of chunk unnamed-chunk-15" width="80%" />


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

<img src="./plots/PCplot_after.png" title="plot of chunk unnamed-chunk-18" alt="plot of chunk unnamed-chunk-18" width="80%" />

# DEG Analysis 

TODO Description here.


```r
q_selection <- NULL # use default setting
deg_q <- sigident::qSelection_(sampleMetadata = sampleMetadata,
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

<img src="./plots/DEG_heatmap.png" title="plot of chunk unnamed-chunk-20" alt="plot of chunk unnamed-chunk-20" width="80%" />

TODO description here


```r
deg_info <- sigident::exportDEGannotations_(mergedset = mergedset)
data.table::fwrite(deg_info, paste0(csvdir, "DEG_info.csv"))
```


```r
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
|path:hsa04512 |ECM-receptor interaction                                      |  88|  8| 0.0003220|
|path:hsa05418 |Fluid shear stress and atherosclerosis                        | 139| 10| 0.0004003|
|path:hsa04061 |Viral protein interaction with cytokine and cytokine receptor | 100|  8| 0.0007640|
|path:hsa04610 |Complement and coagulation cascades                           |  79|  7| 0.0008924|

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

<img src="./plots//hsa04110.png" title="plot of chunk unnamed-chunk-36" alt="plot of chunk unnamed-chunk-36" width="80%" />

```r
filename <- paste0(plotdir, "/", organism, "04110.pathview.png")
```

```r
knitr::include_graphics(filename)
```

<img src="./plots//hsa04110.pathview.png" title="plot of chunk unnamed-chunk-38" alt="plot of chunk unnamed-chunk-38" width="80%" />

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

<img src="./plots/Enriched_GO.png" title="plot of chunk unnamed-chunk-40" alt="plot of chunk unnamed-chunk-40" width="80%" />

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

<img src="./plots/Enriched_KEGG.png" title="plot of chunk unnamed-chunk-42" alt="plot of chunk unnamed-chunk-42" width="80%" />


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

<img src="./plots/CV_lasso.png" title="plot of chunk unnamed-chunk-45" alt="plot of chunk unnamed-chunk-45" width="80%" />


```r
filename <- paste0(plotdir, "ROC_Lasso.min.png")
sigident::createROCplot_(roc = diagnostic_lasso$roc.min,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_Lasso.min.png" title="plot of chunk unnamed-chunk-47" alt="plot of chunk unnamed-chunk-47" width="80%" />


```r
filename <- paste0(plotdir, "ROC_Lasso.1se.png")
sigident::createROCplot_(roc = diagnostic_lasso$roc.1se,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_Lasso.1se.png" title="plot of chunk unnamed-chunk-49" alt="plot of chunk unnamed-chunk-49" width="80%" />


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

<img src="./plots/CV_elasticNet.png" title="plot of chunk unnamed-chunk-51" alt="plot of chunk unnamed-chunk-51" width="80%" />


```r
filename <- paste0(plotdir, "ROC_elasticNet.min.png")
sigident::createROCplot_(roc = diagnostic_elasticnet$roc.min,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_elasticNet.min.png" title="plot of chunk unnamed-chunk-53" alt="plot of chunk unnamed-chunk-53" width="80%" />


```r
filename <- paste0(plotdir, "ROC_elasticNet.1se.png")
sigident::createROCplot_(roc = diagnostic_elasticnet$roc.1se,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_elasticNet.1se.png" title="plot of chunk unnamed-chunk-55" alt="plot of chunk unnamed-chunk-55" width="80%" />


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

<img src="./plots/Gridsearch_model.png" title="plot of chunk unnamed-chunk-57" alt="plot of chunk unnamed-chunk-57" width="80%" />


```r
# plot variable importance of gridsearch
filename <- paste0(plotdir, "Gridsearch_variable_importance.png")
sigident::createGridVarImpPlot_(model = diagnostic_glmGrid$caret.train,
                                filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/Gridsearch_variable_importance.png" title="plot of chunk unnamed-chunk-59" alt="plot of chunk unnamed-chunk-59" width="80%" />


```r
# create roc plot
filename <- paste0(plotdir, "ROC_elasticNet.grid.png")
sigident::createROCplot_(roc = diagnostic_glmGrid$roc.elasticNet,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_elasticNet.grid.png" title="plot of chunk unnamed-chunk-61" alt="plot of chunk unnamed-chunk-61" width="80%" />




# An alternative workflow: one generic function 

In order to simplify the whole abovedescribed workflow, we wrapped all these functions into one big function.

## Preprocessing: 

The preprocessing needs to be conform with the above described approach.

## Run `sigidentMicroarray`-Function 



```r
diagnosticModels <- sigident::sigidentMicroarray(mergeset = mergeset,
                                                 controlname = "Control",
                                                 targetname = "Lung Cancer",
                                                 studyMetadata = studyMetadata,
                                                 sampleMetadata = sampleMetadata,
                                                 species = "Hs",
                                                 OrgDb = "org.Hs.eg.db",
                                                 organism = "hsa",
                                                 pathwayid = "hsa04110",
                                                 deg.q.selection = NULL,
                                                 seed = 111,
                                                 nfolds = 10,
                                                 split = 0.8,
                                                 plotdir = "./plots/",
                                                 csvdir = "./tables/",
                                                 targetcol = "target"))
```

