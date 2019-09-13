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
library(mergeGEO)
```

# Prerequisites and Installation 

## Install mergeGEO package


```r
devtools::install_git("https://gitlab.miracum.org/clearly/mergegeo.git")
```

## Install sigident package


```r
devtools::install_git("https://gitlab.miracum.org/clearly/sigident.git")
```

# Data Preparation 

In order to use this R package and its functions, you need to prepare a merged gene expression dataset. 
This example uses the GEO lung cancer studies "GSE18842", "GSE19804" and "GSE19188", which contain in total 367 samples (197 tumor; 170 non-tumor) and 54,675 transcripts gained by using Affymetrix GeneChip Human Genome U133 Plus 2.0 Array (platform GPL570).

We therefore created the R package 'mergeGEO' that organizes a method to merge GEO datasets, originally described and published by Hughey & Butte (2015) in one R function.
For a detailled background, please refer to the original publication ([https://doi.org/10.1093/nar/gkv229](https://doi.org/10.1093/nar/gkv229)) and its supplemental R files (licensed under GPLv3, [https://doi.org/10.5281/zenodo.16006](https://doi.org/10.5281/zenodo.16006)). 

To get a detailled background on how to use 'mergeGEO', please visit its package vignette: ([link](link))

## Download of GEO datasets using mergeGEO 

During the utilization of mergeGEO, the GEO datasets are downloaded, normalized and batch corrected. In order to provide a working example, we included exemplary files containing the required study metadata and the sample metadata.


```r
studymetadata = "lungcancer_study_metadata.csv"
samplemetadata = "lungcancer_sample_metadata.csv"
studyname <- "lungcancer"
denovo <- F
metadatadir <- "./metadata/"
dir.create(metadatadir)

# copy prepared files to metadatadir
file.copy(from=system.file("./demofiles/lungcancer_study_metadata.csv",
                           package = "mergeGEO"),
          to="./metadata/lungcancer_study_metadata.csv")
file.copy(from=system.file("./demofiles/lungcancer_sample_metadata.csv",
                           package = "mergeGEO"),
          to="./metadata/lungcancer_sample_metadata.csv")

mergeset <- mergeGEO::mergeGEO(studymetadata = studymetadata,
                               samplemetadata = samplemetadata,
                               studyname = studyname,
                               denovo = denovo,
                               metadatadir = metadatadir)
```


## Loading Metadata 

TODO Description here.


```r
studyMetadata <- mergeGEO::readStudyMetadata_(studymetadataFilename = paste0(metadatadir, studymetadata))
sampleMetadata <- mergeGEO::readSampleMetadata_(samplemetadataFilename = paste0(metadatadir, samplemetadata),
                                             studyMetadata = studyMetadata)
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
targetcol <- "target"
controlname <- "Control"
targetname <- "Lung Cancer"
species <- "Hs"
OrgDb <- "org.Hs.eg.db"
organism <- "hsa"
pathwayid <- "hsa04110"
seed <- 111
split <- 0.8
```

First, a boxplot is created with the included samples on the x-axis and the standardised expression values on the y-axis. `mergeset` results as output from the function `mergeGEO()` and represents a matrix containing the genes (Entrez ID) as rows and the samples as columns.


```r
filename <- paste0(plotdir, "import_boxplot.png")
sigident::createImportBoxplot_(mergeset, filename)
```


```r
knitr::include_graphics(filename)
```

<img src="./plots/import_boxplot.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" width="80%" />


# Create experiment desciption and visualize batch effects  

TODO Description here.


```r
dd <- sigident::createDiagnosisDesign_(sampleMetadata = sampleMetadata,
                                       studyMetadata = studyMetadata,
                                       controlname = controlname,
                                       targetname = targetname,
                                       targetcol = targetcol)
diagnosis <- dd$diagnosis
design <- dd$design

batch <- sigident::createBatch_(studyMetadata = studyMetadata,
                                sampleMetadata = sampleMetadata)

gPCA_after <- sigident::batchCorrection_(mergeset = mergeset,
                                         batch = batch)
filename <- paste0(plotdir, "PCplot_after.png")
sigident::createBatchPlot_(correction_obj = gPCA_after,
                           filename = filename,
                           time = "after")
```


```r
knitr::include_graphics(filename)
```

<img src="./plots/PCplot_after.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" width="80%" />

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

<img src="./plots/DEG_heatmap.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" width="80%" />

<!-- TODO start here 13.9. Lorenz -->

<!-- ```{r} -->
<!-- deg_info <- exportDEGannotations_(mergeset = mergeset, -->
<!--                                   genes = genes) -->
<!-- data.table::fwrite(deg_info, paste0(csvdir, "DEG_info.csv")) -->
<!-- ``` -->

<!-- ```{r} -->
<!-- knitr::kable(head(deg_info)) -->
<!-- ``` -->

<!-- ```{r} -->
<!-- deg_results <- limmaTopTable_(BatchRemovedExprs = combat, -->
<!--                               design = design, -->
<!--                               qValue = deg_q) -->
<!-- data.table::fwrite(deg_results, paste0(csvdir, "DEG_results.csv")) -->
<!-- ``` -->

<!-- ```{r} -->
<!-- knitr::kable(head(deg_results)) -->
<!-- ``` -->

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



|           |Term                             |Ont |    N|  DE| P.DE|
|:----------|:--------------------------------|:---|----:|---:|----:|
|GO:0005576 |extracellular region             |CC  | 4295| 153|    0|
|GO:0044421 |extracellular region part        |CC  | 3312| 128|    0|
|GO:0050896 |response to stimulus             |BP  | 8836| 236|    0|
|GO:0008283 |cell proliferation               |BP  | 1969|  89|    0|
|GO:0042127 |regulation of cell proliferation |BP  | 1664|  80|    0|
|GO:0031012 |extracellular matrix             |CC  |  468|  39|    0|

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



|              |Pathway                                |   N| DE|      P.DE|
|:-------------|:--------------------------------------|---:|--:|---------:|
|path:hsa05144 |Malaria                                |  49|  9| 0.0000010|
|path:hsa04110 |Cell cycle                             | 124| 11| 0.0000874|
|path:hsa04657 |IL-17 signaling pathway                |  93|  9| 0.0002005|
|path:hsa04974 |Protein digestion and absorption       |  95|  9| 0.0002358|
|path:hsa05418 |Fluid shear stress and atherosclerosis | 139| 11| 0.0002424|
|path:hsa04610 |Complement and coagulation cascades    |  79|  8| 0.0003330|

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



|           |Term                          |Ont |     N|   Up| Down|      P.Up| P.Down|
|:----------|:-----------------------------|:---|-----:|----:|----:|---------:|------:|
|GO:0002376 |immune system process         |BP  |  2678|  893| 1062| 1.0000000|      0|
|GO:0006955 |immune response               |BP  |  1821|  547|  787| 1.0000000|      0|
|GO:0050896 |response to stimulus          |BP  |  7683| 3013| 2466| 1.0000000|      0|
|GO:0001775 |cell activation               |BP  |  1275|  370|  576| 1.0000000|      0|
|GO:0051716 |cellular response to stimulus |BP  |  6236| 2458| 2043| 0.9999997|      0|
|GO:0065007 |biological regulation         |BP  | 10319| 4199| 3113| 0.9999517|      0|


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



|              |Pathway                                 |   N| Up| Down|      P.Up| P.Down|
|:-------------|:---------------------------------------|---:|--:|----:|---------:|------:|
|path:hsa04380 |Osteoclast differentiation              | 126| 25|   75| 1.0000000|      0|
|path:hsa04010 |MAPK signaling pathway                  | 283| 97|  133| 0.9999954|      0|
|path:hsa04640 |Hematopoietic cell lineage              |  92| 16|   56| 1.0000000|      0|
|path:hsa05332 |Graft-versus-host disease               |  36|  1|   29| 1.0000000|      0|
|path:hsa05166 |Human T-cell leukemia virus 1 infection | 214| 79|  103| 0.9989278|      0|
|path:hsa04611 |Platelet activation                     | 119| 26|   66| 1.0000000|      0|

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

<img src="./plots//hsa04110.png" title="plot of chunk unnamed-chunk-24" alt="plot of chunk unnamed-chunk-24" width="80%" />

```r
filename <- paste0(plotdir, "/", organism, "04110.pathview.png")
```

```r
knitr::include_graphics(filename)
```

<img src="./plots//hsa04110.pathview.png" title="plot of chunk unnamed-chunk-26" alt="plot of chunk unnamed-chunk-26" width="80%" />

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

<img src="./plots/Enriched_GO.png" title="plot of chunk unnamed-chunk-28" alt="plot of chunk unnamed-chunk-28" width="80%" />

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

<img src="./plots/Enriched_KEGG.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" width="80%" />


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

<img src="./plots/CV_lasso.png" title="plot of chunk unnamed-chunk-33" alt="plot of chunk unnamed-chunk-33" width="80%" />


```r
filename <- paste0(plotdir, "ROC_Lasso.min.png")
sigident::createROCplot_(roc = diagnostic_lasso$roc.min,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_Lasso.min.png" title="plot of chunk unnamed-chunk-35" alt="plot of chunk unnamed-chunk-35" width="80%" />


```r
filename <- paste0(plotdir, "ROC_Lasso.1se.png")
sigident::createROCplot_(roc = diagnostic_lasso$roc.1se,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_Lasso.1se.png" title="plot of chunk unnamed-chunk-37" alt="plot of chunk unnamed-chunk-37" width="80%" />


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

<img src="./plots/CV_elasticNet.png" title="plot of chunk unnamed-chunk-39" alt="plot of chunk unnamed-chunk-39" width="80%" />


```r
filename <- paste0(plotdir, "ROC_elasticNet.min.png")
sigident::createROCplot_(roc = diagnostic_elasticnet$roc.min,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_elasticNet.min.png" title="plot of chunk unnamed-chunk-41" alt="plot of chunk unnamed-chunk-41" width="80%" />


```r
filename <- paste0(plotdir, "ROC_elasticNet.1se.png")
sigident::createROCplot_(roc = diagnostic_elasticnet$roc.1se,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_elasticNet.1se.png" title="plot of chunk unnamed-chunk-43" alt="plot of chunk unnamed-chunk-43" width="80%" />


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

<img src="./plots/Gridsearch_model.png" title="plot of chunk unnamed-chunk-45" alt="plot of chunk unnamed-chunk-45" width="80%" />


```r
# plot variable importance of gridsearch
filename <- paste0(plotdir, "Gridsearch_variable_importance.png")
sigident::createGridVarImpPlot_(model = diagnostic_glmGrid$caret.train,
                                filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/Gridsearch_variable_importance.png" title="plot of chunk unnamed-chunk-47" alt="plot of chunk unnamed-chunk-47" width="80%" />


```r
# create roc plot
filename <- paste0(plotdir, "ROC_elasticNet.grid.png")
sigident::createROCplot_(roc = diagnostic_glmGrid$roc.elasticNet,
                         filename = filename)
```

```r
knitr::include_graphics(filename)
```

<img src="./plots/ROC_elasticNet.grid.png" title="plot of chunk unnamed-chunk-49" alt="plot of chunk unnamed-chunk-49" width="80%" />




# An alternative workflow: one generic function 

In order to simplify the whole abovedescribed workflow, we wrapped all these functions into one big function.

## Preprocessing: 

TODO Description here.

### Merge data sets  

TODO Description here.


```r
studymetadata = "lungcancer_study_metadata.csv"
samplemetadata = "lungcancer_sample_metadata.csv"
studyname <- "lungcancer"
denovo <- TRUE
metadatadir <- "./metadata/"
dir.create(metadatadir)

# copy prepared files to metadatadir
file.copy(from=system.file("./demofiles/lungcancer_study_metadata.csv",
                           package = "mergeGEO"),
          to="./metadata/lungcancer_study_metadata.csv")
file.copy(from=system.file("./demofiles/lungcancer_sample_metadata.csv",
                           package = "mergeGEO"),
          to="./metadata/lungcancer_sample_metadata.csv")

mergeset <- mergeGEO::mergeGEO(studymetadata = studymetadata,
                               samplemetadata = samplemetadata,
                               studyname = studyname,
                               denovo = denovo,
                               metadatadir = metadatadir)
```

### Read metadata 

TODO Description here.


```r
studyMetadata <- mergeGEO::readStudyMetadata_(studymetadataFilename = paste0(metadatadir,
                                                                             studymetadata))
sampleMetadata <- mergeGEO::readSampleMetadata_(samplemetadataFilename = paste0(metadatadir,
                                                                                samplemetadata),
                                                studyMetadata = studyMetadata)
```

## Run `sigidentMicroarray`-Function 

TODO Description here.


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

