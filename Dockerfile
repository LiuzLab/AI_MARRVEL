FROM python:3.8

# Install dependencies
RUN apt update && apt install -y \
    curl \
    git \
    build-essential \
    libncurses-dev \
    zlib1g-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    tabix \
    jq

# Install python 3.8 dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip3 install --upgrade pip
RUN pip3 install -r /opt/requirements.txt

# Install bcftools
RUN curl -OL https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
RUN mv bcftools-1.20.tar.bz2 /opt/bcftools-1.20.tar.bz2
RUN tar -xf /opt/bcftools-1.20.tar.bz2 -C /opt/ && \
  rm /opt/bcftools-1.20.tar.bz2 && \
  cd /opt/bcftools-1.20 && \
  ./configure && \
  make && \
  make install && \
  rm -rf /opt/bcftools-1.20

# Install bedtools
RUN curl -OL https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
RUN mv bedtools.static.binary /run/bedtools
RUN chmod a+x /run/bedtools
