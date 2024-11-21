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

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
RUN mv bcftools-1.20.tar.bz2 /opt/bcftools-1.20.tar.bz2
RUN tar -xf /opt/bcftools-1.20.tar.bz2 -C /opt/ && \
  rm /opt/bcftools-1.20.tar.bz2 && \
  cd /opt/bcftools-1.20 && \
  ./configure && \
  make && \
  make install && \
  rm -rf /opt/bcftools-1.20

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
RUN mv bedtools.static.binary /run/bedtools
RUN chmod a+x /run/bedtools



