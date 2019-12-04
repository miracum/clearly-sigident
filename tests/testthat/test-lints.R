context("lints")

if (dir.exists("../../00_pkg_src")) {
  prefix <- "../../00_pkg_src/sigident/"
} else if (dir.exists("../../R")) {
  prefix <- "../../"
} else if (dir.exists("./R")) {
  prefix <- "./"
}

test_that(
  desc = "test lints",
  code = {
    lintlist <- list(
      "R" = list(
        "compare_diagnostic_models.R" = NULL,
        "DEG.R" = NULL,
        "diagnostic_sig.R" = NULL,
        "enrichment.R" = NULL,
        "otherfunctions.R" = NULL,
        "plotting.R" = NULL,
        "prognostic_sig.R" = NULL,
        "sigidentDEG.R" = "snake_case",
        "sigidentDiagnostic.R" = "snake_case",
        "sigidentEnrichment.R" = "snake_case",
        "sigidentPrognostic.R" = "snake_case",
        "utils.R" = NULL
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

        # skip on covr
        skip_on_covr()

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
