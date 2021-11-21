FROM rocker/r-base

USER root

RUN apt-get update --allow-releaseinfo-change && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    librsvg2-dev \
    libavfilter-dev \
    libwebp-dev \
    cargo \
    libpoppler-cpp-dev \
    libmagick++-dev

RUN Rscript -e 'install.packages("remotes"); remotes::install_github("icbi-lab/immunedeconv");'
RUN Rscript -e 'install.packages(c("tidyverse", "dplyr"), repos = "http://cran.rstudio.com/");'

COPY ./immunedeconv.R /app/immunedeconv.R

CMD /bin/bash