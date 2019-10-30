FROM rocker/tidyverse:3.6.1

# set cran repo
RUN echo "options('repos' = 'https://ftp.fau.de/cran/')" >> /usr/local/lib/R/etc/Rprofile.site

# install system dependencies
RUN apt-get install -y --no-install-recommends \
    libjpeg-dev # dependency of qpdf

ARG PKG_NAME=sigident
ADD ${PKG_NAME} ${PKG_NAME}

# install package dependencies
RUN R -e "devtools::install_deps(pkg = '${PKG_NAME}', upgrade = 'always')"
RUN R -e "devtools::install_dev_deps(pkg = '${PKG_NAME}', upgrade = 'always')"

# clean up
RUN rm -rf ${PKG_NAME}

CMD ["R"]
