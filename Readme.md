# sigident (!!! currently under development !!!)

This is the repository of the R package 'sigident'. It provides core functionalities to identify diagnostic and prognostic signatures from gene expression datasets.

Currently implemented features are:

- merging of microarray datasets
- DEG analysis
- functional analysis (gene enrichment)
- Identification of diagnostic signatures, using
  + Lasso regression
  + Elastic net regression
  + glmnet (grid search for best alpha and lambda)
- Identification of prognostic signatures

## Installation

You can install the development version of *sigident* with:

``` r
options('repos' = 'https://ftp.fau.de/cran/')
install.packages("devtools")
devtools::install_git("https://gitlab.miracum.org/clearly/sigident.git")
```

## Example

Please view the package's vignette to see a detailled description how to perform merging, signature and functional analyses of microarray data.

