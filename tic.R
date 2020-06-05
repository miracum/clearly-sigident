# installs dependencies, runs R CMD check, runs covr::codecov()
do_package_checks(
  args = "--no-vignettes",
  build_args = "--no-build-vignettes"
)
