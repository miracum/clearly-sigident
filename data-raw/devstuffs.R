packagename <- "sigident"

# remove existing description object
unlink("DESCRIPTION")
# Create a new description object
my_desc <- desc::description$new("!new")
# Set your package name
my_desc$set("Package", packagename)
# Set author names
my_desc$set_authors(c(
  person("Lorenz A.", "Kapsner", email = "lorenz.kapsner@uk-erlangen.de", role = c('cre', 'aut'),
         comment = c(ORCID = "0000-0003-1866-860X")),
  person("Johannes", "Vey", role = c('aut'),
         comment = c(ORCID = "0000-0002-2610-9667")),
  person("Meik", "Kunz", role = c('ctb')),
  person("Andreas", "Pittroff", role = c('ctb')),
  person("Tobias", "Leis", role = c('ctb')),
  person("Julian", "Henke", role = c('ctb')),
  person("Nicholas", "Dickel", role = c('ctb'))
))
# Remove some author fields
my_desc$del("Maintainer")
# Vignette Builder
my_desc$set("VignetteBuilder" = "knitr")
# Set the version
my_desc$set_version("0.1.2")
# The title of your package
my_desc$set(Title = "Signature Analyses in Genomic Expression Sets")
# The description of your package
my_desc$set(Description = "Identify diagnostic and prognostic signatures from gene expression datasets.")
# The description of your package
my_desc$set("Date" = as.character(Sys.Date()))
# The urls
my_desc$set("URL", "https://github.com/miracum/clearly-sigident")
my_desc$set("BugReports", "https://github.com/miracum/clearly-sigident/issues")
# License
my_desc$set("License", "GPL-3")

# BioConductor stuff
my_desc$set("biocViews" = "")


# Save everyting
my_desc$write(file = "DESCRIPTION")

# License
#usethis::use_gpl3_license(name="Universitätsklinikum Erlangen")


# add Imports and Depends
# Listing a package in either Depends or Imports ensures that it’s installed when needed
# Imports just loads the package, Depends attaches it
# Loading will load code, data and any DLLs; register S3 and S4 methods; and run the .onLoad() function.
##      After loading, the package is available in memory, but because it’s not in the search path,
##      you won’t be able to access its components without using ::.
##      Confusingly, :: will also load a package automatically if it isn’t already loaded.
##      It’s rare to load a package explicitly, but you can do so with requireNamespace() or loadNamespace().
# Attaching puts the package in the search path. You can’t attach a package without first loading it,
##      so both library() or require() load then attach the package.
##      You can see the currently attached packages with search().

# Depends

# Imports (CRAN packages)
usethis::use_package("data.table", type="Imports")
usethis::use_package("grDevices", type="Imports")
usethis::use_package("magrittr", type="Imports")
usethis::use_package("stats", type="Imports")
usethis::use_package("graphics", type="Imports")
usethis::use_package("caret", type="Imports")
usethis::use_package("glmnet", type="Imports")
usethis::use_package("pROC", type="Imports")
usethis::use_package("doParallel", type="Imports")
usethis::use_package("survival", type="Imports")
usethis::use_package("plyr", type="Imports")
usethis::use_package("survminer", type="Imports")
usethis::use_package("randomForest", type="Imports")
usethis::use_package("gbm", type="Imports")
usethis::use_package("sigident.preproc", type="Imports")
usethis::use_package("sigident.func", type="Imports")

# Bioconductor
# https://github.com/r-lib/devtools/issues/700

# Suggests
usethis::use_package("e1071", type="Suggests")
usethis::use_package("testthat", type = "Suggests")
usethis::use_package("devtools", type = "Suggests")
usethis::use_package("rmarkdown", type = "Suggests")
usethis::use_package("qpdf", type = "Suggests")
usethis::use_package("knitr", type = "Suggests")
usethis::use_package("lintr", type = "Suggests")


# Development package
preproc_tag <- "v0.0.5"
func_tag <- "v0.0.4"

devtools::install_github("miracum/clearly-sigident.preproc", ref = preproc_tag, upgrade = "always")
devtools::install_github("miracum/clearly-sigident.func", ref = func_tag, upgrade = "always")

# https://cran.r-project.org/web/packages/devtools/vignettes/dependencies.html
desc::desc_set_remotes(c(
  paste0(
    "github::miracum/clearly-sigident.preproc@", preproc_tag),
  # sigident.func is only required for the vignettes (therefore a "suggests"-package)
  paste0(
    "github::miracum/clearly-sigident.func@", func_tag)
),
file = usethis::proj_get())

# buildignore and gitignore
usethis::use_build_ignore("docker")
usethis::use_build_ignore("vignettes/geodata")
usethis::use_build_ignore("vignettes/plots")
usethis::use_build_ignore("vignettes/csv")
usethis::use_build_ignore("vignettes/GSE19188")
usethis::use_build_ignore("tests/testthat/testdata")
usethis::use_build_ignore("tests/testthat/plots")
usethis::use_build_ignore(".RData")
usethis::use_build_ignore("data")
usethis::use_build_ignore("docker")
usethis::use_build_ignore("Glio_main.rds")
usethis::use_build_ignore("metadatadir")
usethis::use_build_ignore("vignettes/metadata")
usethis::use_build_ignore("vignettes/preview.dir")
usethis::use_build_ignore("vignettes/*.csv")
usethis::use_build_ignore("vignettes/figure")
usethis::use_build_ignore("vignettes/tables")
usethis::use_build_ignore("geodata")
usethis::use_build_ignore("csv")
usethis::use_build_ignore("plots")
usethis::use_build_ignore("Readme.md")
usethis::use_build_ignore("infos.R")
usethis::use_build_ignore(".gitlab-ci.yml")
usethis::use_build_ignore("ci")
usethis::use_build_ignore(".vscode")
usethis::use_build_ignore(".lintr")
usethis::use_build_ignore("doc")
usethis::use_build_ignore("Meta")


usethis::use_git_ignore("/*")
usethis::use_git_ignore("/*/")
usethis::use_git_ignore("*.log")
usethis::use_git_ignore("!/.gitignore")
usethis::use_git_ignore("!/.gitlab-ci.yml")
usethis::use_git_ignore("!/data-raw/")
usethis::use_git_ignore("!/DESCRIPTION")
usethis::use_git_ignore("!/inst/")
usethis::use_git_ignore("!/LICENSE.md")
usethis::use_git_ignore("!/man/")
usethis::use_git_ignore("!NAMESPACE")
usethis::use_git_ignore("!/R/")
usethis::use_git_ignore("!/ci/")
usethis::use_git_ignore("!/README.md")
usethis::use_git_ignore("!/tests/")
usethis::use_git_ignore("/.Rhistory")
usethis::use_git_ignore("/.Rproj*")
usethis::use_git_ignore("/.RData")
usethis::use_git_ignore("!/vignettes/")
usethis::use_git_ignore("/vignettes/*")
usethis::use_git_ignore("!/vignettes/*.Rmd")
usethis::use_git_ignore("tests/testthat/testdata")
usethis::use_git_ignore("tests/testthat/plots")
usethis::use_git_ignore("geodata")
usethis::use_git_ignore("csv")
usethis::use_git_ignore("plots")
usethis::use_git_ignore("!/*.Rproj")
usethis::use_git_ignore(".Rproj*")
usethis::use_git_ignore("*.rds")
usethis::use_git_ignore("/.vscode")
usethis::use_git_ignore("!/.lintr")
