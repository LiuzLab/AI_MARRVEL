FROM rocker/tidyverse

# Install R libs
RUN R -e "install.packages('ontologyIndex',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ontologySimilarity',dependencies=TRUE, repos='http://cran.rstudio.com/')"
