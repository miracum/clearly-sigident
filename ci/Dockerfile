FROM rocker/tidyverse:3.6.1

RUN echo "options('repos' = 'https://ftp.fau.de/cran/')" >> /usr/local/lib/R/etc/Rprofile.site
RUN apt-get install libjpeg-dev -y --no-install-recommends # dependency of qpdf

RUN R -e "devtools::install_deps(pkg = '.', upgrade = 'always')"
RUN R -e "devtools::install_dev_deps(pkg = '.', upgrade = 'always')"

CMD ["R"]
