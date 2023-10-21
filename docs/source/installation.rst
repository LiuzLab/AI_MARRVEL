.. _install:

*************
Installation
*************

AIM is constructed with R and Python in different environment.

To free users from setting up all the environment, we provide the software as docker image.


AIM-Lite
=============
AIM-Lite takes VCF and HPO as input to generate predictions. 

It is the samllest function unit of AI-MARRVEL software.

.. note::

   An more powerful batch-enabled AIM docker image is under development.

**Pull AIM-Lite Docker from Docker Hub**
.. code-block:: bash
    
    docker pull chaozhongliu/aim-lite:latest


**Build Local AIM-Lite Docker**

If you want to keep a stable local version of AIM, follow the instruction below.

.. code-block:: bash
    
    # Clone repositories
    git clone https://github.com/LiuzLab/AI_MARRVEL.git
    cd AI_MARRVEL

    # Download necessary software
    #   bcftools
    wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
    #   bedtools
    wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
    mv bedtools.static.binary run/bedtools
    chmod a+x run/bedtools

    # Build docker image (takes some time)
    docker build -t aim-lite .


Install Required Data Dependencies
===================================
AIM requires multiple database sources to annotate variants, we have compiled all for users to download.

.. warning::

   HGMD-related database is included but not provided in the public download link, for license reasons.



