context("lints")

#% prefix <- "./"

prefix <- "../../00_pkg_src/sigident/"

test_that(
  desc = "test lints",
  code = {
    lintlist <- list(
      "R" = list(
        "batchcorrection.R" = NULL,
        "compare_diagnostic_models.R" = NULL,
        "data_handling.R" = NULL,
        "deg.R" = NULL,
        "diagnostic_sig.R" = NULL,
        "enrichment.R" = NULL,
        "otherfunctions.R" = NULL,
        "plotting.R" = NULL,
        "prognostic_sig.R" = NULL
      ),
      "tests/testthat" = list(
        "test-lints.R" = NULL
      )
    )
    for (directory in names(lintlist)) {
      print(directory)
      for (fname in names(lintlist[[directory]])) {
        print(fname)
        #% print(list.files(prefix))

        # skip on cran
        #% skip_on_cran()

        lintr::expect_lint(
          file = paste0(
            prefix,
            directory,
            "/",
            fname
          ),
          checks = lintlist[[directory]][[fname]]
        )
      }
    }
  }
)
