FROM rocker/rstudio:3.6.1

COPY install.R /home/rstudio/

RUN if [ -f /home/rstudio/install.R ]; then R --quiet -f /home/rstudio/install.R; fi
