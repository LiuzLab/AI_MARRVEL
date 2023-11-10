.. _install:

*************
Installation
*************

To streamline the user experience and eliminate the need for environment setup, AIM is delivered as a Docker image


AIM-Lite
=============
AIM-Lite processes VCF and HPO inputs to output predictions, representing the most compact functional component of the AI-MARRVEL suite.

.. note::

   An more powerful batch-enabled AIM docker image is under development.

**Pull AIM-Lite Docker from Docker Hub**

.. code-block:: bash
    
    docker pull chaozhongliu/aim-lite:latest


**Alternatively, Build Local AIM-Lite Docker**

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
AIM utilizes various databases for variant annotation, all of which have been compiled and made available for user download.

We use Globus for data accessing `here <https://app.globus.org/file-manager?origin_id=6810458e-b702-423f-9f0c-070c1691482d&origin_path=%2F>`_

.. warning::

   Due to licensing restrictions, the HGMD-related database is not included in the public download. 
   However, AIM operates effectively without this data.

   To prepare HGMD database, see below.


Prepare HGMD-related data
===================================

**HGMD databse VCF**

Copy your HGMD ``bgzip`` compressed VCF file and the index file (\*.tbi) into your database folder under ``vep/hg19(hg38)/``, and rename it as ``HGMD_Pro_2022.2_hg19(hg38).vcf.gz(.tbi)``.


**HGMD Phenotype Information**

Get the ``hgmd_phenbase-*.dump.sql`` from HGMD.

Run `HGMD_phenbase.sh <https://github.com/LiuzLab/AI_MARRVEL/blob/main/utils/HGMD_phenbase.sh>`_ to get ``HGMD_phen.tsv``. Move the file to your database folder under ``omim_annotate/hg19(hg38)/``
