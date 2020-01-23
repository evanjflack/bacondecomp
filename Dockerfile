FROM rocker/r-base

RUN R -e "install.packages('knitr', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('rmarkdown', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('testthat', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('ggplot2', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('covr', repos = 'http://cran.us.r-project.org')"
