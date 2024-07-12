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

    # Build docker image (takes some time)
    docker build -t aim-lite .


Install Required Data Dependencies
=================================== 
AIM utilizes various databases for variant annotation, all of which have been compiled and are available for download. We use AWS S3 for data access, and the data can be downloaded by following these steps:

1. **Install the AWS CLI**:
   Follow the instructions provided in the `AWS CLI Installation Guide <https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html>`_


2. **Navigate to Your Desired Directory**:
   Change to the directory where you want your data dependencies downloaded. For example, in Ubuntu, use:

   .. code-block:: bash

      cd <desired/folder/path>

3. **Use the following command to sync the S3 bucket to your local directory**:

   .. code-block:: bash

      aws s3 sync s3://aim-data-dependencies-public . --no-sign-request


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
