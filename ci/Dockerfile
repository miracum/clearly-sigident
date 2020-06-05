FROM rocker/verse:4.0.0

# build args
ARG PKG_NAME

# set cran repo
RUN echo "options('repos' = 'https://cloud.r-project.org/')" >> /usr/local/lib/R/etc/Rprofile.site

# install system dependencies
RUN apt-get install -y --no-install-recommends \
    libjpeg-dev # dependency of qpdf

ADD ${PKG_NAME} ${PKG_NAME}

# install package dependencies
RUN R -e "devtools::install_git(url = 'https://gitlab.miracum.org/clearly/sigident.preproc', ref = 'development', upgrade = 'always')"
RUN R -e "devtools::install_dev_deps(pkg = '${PKG_NAME}', upgrade = 'always')"

# clean up
RUN rm -rf ${PKG_NAME}
