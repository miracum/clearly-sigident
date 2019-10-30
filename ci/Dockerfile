FROM rocker/tidyverse:3.6.1

RUN echo "options('repos' = 'https://ftp.fau.de/cran/')" >> /usr/local/lib/R/etc/Rprofile.site
RUN apt-get install libjpeg-dev -y --no-install-recommends # dependency of qpdf

ARG PKG_NAME=sigident
RUN ls -al && cd ../ && ls -al
ADD ../${PKG_NAME} ${PKG_NAME}

RUN R -e "devtools::install_deps(pkg = '${PKG_NAME}', upgrade = 'always')"
RUN R -e "devtools::install_dev_deps(pkg = '${PKG_NAME}', upgrade = 'always')"

CMD ["R"]
