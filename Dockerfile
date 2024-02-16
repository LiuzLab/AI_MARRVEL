# List of dependencies installed
#	bcftools v1.9-1-deb_cv1
#	ensemblorg/ensembl-vep:release_104.3
#	python 2.7
#	python 3.8.13
#	R 4.2.1

# dc7a542b7e1f
#FROM ubuntu:18.04
FROM ensemblorg/ensembl-vep:release_104.3
USER root
ENV DEBIAN_FRONTEND noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    curl \
    git \
    build-essential \
    libncurses-dev \
    zlib1g-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    tabix \
    python2.7 \
    python3.8 \
    python3.8-dev \
    python3.8-distutils \
    python3-apt

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python3.8 get-pip.py

# Install python 3.8 dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip3 install --upgrade pip
RUN pip3 install -r /opt/requirements.txt
RUN pip3 install bgzip


# Install R
RUN apt-get update
RUN apt install -y --no-install-recommends software-properties-common dirmngr
# Add the keys
RUN apt install wget
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 || \
    apt-key adv --keyserver ha.pool.sks-keyservers.net --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 || \
    apt-key adv --keyserver pgp.mit.edu --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 || \
    apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 || \
    apt-key adv --keyserver keyserver.pgp.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt install -y r-base r-base-core r-recommended r-base-dev

# Install R libs
RUN R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ontologyIndex',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ontologySimilarity',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"



# Install bcftools
COPY bcftools-1.9.tar.bz2 /opt/bcftools-1.9.tar.bz2
RUN tar -xf /opt/bcftools-1.9.tar.bz2 -C /opt/ && \
  rm /opt/bcftools-1.9.tar.bz2 && \
  cd /opt/bcftools-1.9 && \
  ./configure && \
  make && \
  make install && \
  rm -rf /opt/bcftools-1.9


# Copy the pipeline into Docker image
COPY run /run/
RUN chmod +x /run/proc.sh



