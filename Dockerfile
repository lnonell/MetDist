#################################################################
# Dockerfile
# Author:           Lara Nonell
# Version:          1
# Date              210419
# Software:         R
# Description:      Packages to test methylation distributions 
# Tags:             None, for the moment
# Base Image:       R
#################################################################

FROM r-base:3.5.1

RUN apt-get update && apt-get install -y \
    r-cran-xml \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    ghostscript
ENV PATH=pkg-config:$PATH

RUN install2.r --error --deps TRUE \
    doParallel \
	foreach \
	fitdistrplus \
	VGAM \
	ZOIP \
	aod \
	simplexreg \
	betareg \
	ZOIP \
	gamlss \
	quantreg \
    && rm -rf /tmp/downloaded_packages/

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("Biobase");biocLite("limma");biocLite("biomaRt");biocLite("RnBeads");biocLite("RnBeads.hg38")'

##That's all for the moment
