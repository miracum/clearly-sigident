# sigident (!!! under development !!!)

<!-- badges: start -->
[![R CMD Check via {tic}](https://github.com/miracum/clearly-sigident/workflows/R%20CMD%20Check%20via%20{tic}/badge.svg?branch=master)](https://github.com/miracum/clearly-sigident/actions)
[![linting](https://github.com/miracum/clearly-sigident/workflows/lint/badge.svg?branch=master)](https://github.com/miracum/clearly-sigident/actions)
[![test-coverage](https://github.com/miracum/clearly-sigident/workflows/test-coverage/badge.svg?branch=master)](https://github.com/miracum/clearly-sigident/actions)
[![codecov](https://codecov.io/gh/miracum/clearly-sigident/branch/master/graph/badge.svg)](https://codecov.io/gh/miracum/clearly-sigident)
[![pipeline status](https://gitlab.miracum.org/clearly/sigident/badges/master/pipeline.svg)](https://gitlab.miracum.org/clearly/sigident/commits/master)
[![coverage report](https://gitlab.miracum.org/clearly/sigident/badges/master/coverage.svg)](https://gitlab.miracum.org/clearly/sigident/commits/master)
<!-- badges: end -->

This is the repository of the R package 'sigident'. It provides core functionalities to identify diagnostic and prognostic signatures from gene expression datasets.

Currently implemented features are:

- merging of microarray datasets (via the R package [`sigident.preproc`](https://github.com/miracum/clearly-sigident.preproc.git))
- DEG analysis and functional analysis (via the R package [`sigident.func`](https://github.com/miracum/clearly-sigident.func.git))
- Identification and validation of diagnostic signatures, using
  + Lasso regression
  + Elastic net regression
  + glmnet (grid search for best alpha and lambda)
- Identification of prognostic signatures

# Installation

You can install *sigident* with the following commands in R:

```r
install.packages("devtools")
devtools::install_github("miracum/clearly-sigident")
```

The version of the package, which was used for the publication [A Toolbox for Functional Analysis and the Systematic Identification of Diagnostic and Prognostic Gene Expression Signatures Combining Meta-Analysis and Machine Learning](https://www.mdpi.com/2072-6694/11/10/1606) can anytime be reproduced using the version tag *v0.0.2* during the installation process:

```r
devtools::install_github("miracum/clearly-sigident", ref = "v0.0.2")
```

# Example

Please view the [package's vignette](vignettes/) to see a detailled description how to prepare datasets in order to be suitable for usage with the `sigident` package and to learn, how to perform merging, signature and functional analyses of microarray data.

Since the building the package vignette takes rather long (~ 40 min.), we provide the already built vignettes in [this repository](https://github.com/miracum/clearly-sigident_vignettes). 

# Notice 

The *sigident* package is under active development and not on CRAN yet - this means, that from time to time, the API can break, due to extending and modifying its functionality. It can also happen, that previoulsy included functions and/or function arguments are no longer supported. 
However, a detailed package vignette will be provided alongside with every major change in order to describe the currently supported workflow.

# Citation  

To cite the 'sigident' package in publications, please use: 

```
@Article{,
    title = {A Toolbox for Functional Analysis and the Systematic Identification of Diagnostic and Prognostic Gene Expression Signatures Combining Meta-Analysis and Machine Learning},
    volume = {11},
    doi = {10.3390/cancers11101606},
    pages = {14},
    number = {1606},
    journal = {Cancers},
    author = {Johannes Vey and Lorenz A. Kapsner and Maximilian Fuchs and Philipp Unberath and Giulia Veronesi and Meik Kunz},
    year = {2019},
}
```
and

```
@Manual{,
    title = {sigident: Signature Analyses in Genomic Expression Sets},
    author = {Lorenz A. Kapsner and Johannes Vey and Meik Kunz and Andreas Pittroff},
    year = {2019},
    note = {R package version 0.0.1.9000},
    url = {https://gitlab.miracum.org/clearly/sigident},
}
```

# More Infos:

- about CLEARLY: [https://www.transcanfp7.eu/index.php/abstract/clearly.html](https://www.transcanfp7.eu/index.php/abstract/clearly.html)
- about MIRACUM: [https://www.miracum.org/](https://www.miracum.org/)
- about the Medical Informatics Initiative: [https://www.medizininformatik-initiative.de/index.php/de](https://www.medizininformatik-initiative.de/index.php/de)
