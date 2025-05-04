#
# Create mpem/global image
#

# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:4.4.3

## create directories
RUN mkdir -p /Code

## copy global files
COPY ./Runfiles-Docker/install_packages.r Code/.
COPY ./Code/*.r Code/.

## install R-packages
RUN Rscript /Code/install_packages.r
