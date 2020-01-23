FROM rocker/r-base

COPY install.R /home/rstudio/
RUN chown -R  rstudio /home/rstudio/

RUN if [ -f /home/rstudio/install.R ]; then R --quiet -f /home/rstudio/install.R; fi
