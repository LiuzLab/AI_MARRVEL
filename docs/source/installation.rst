.. _install:

*************
Installation
*************


Installation Java and nextflow
=============================

It is highly recommend to install java using `SDKMAN <https://sdkman.io/install/>``

   .. code-block:: bash
   curl -s https://get.sdkman.io | bash # install sdkman
   sdk install java 17.0.10-tem         # install java 17
   java -version                        # confirm the correct version of java is installed

After then install nextflow with the command-line

   .. code-block:: bash
   curl -s https://get.nextflow.io | bash # donwload nextflow
   chmod +x nextflow                      # make nextflow excutable
   sudo mv nextflow /usr/local/bin        # move to an executable path (in $PATH)
   nextflow info                          # confirm the installation is done.

## Running AIM

Use following command-line to run AIM.

   .. code-block:: bash
   nextflow run Liuzlab/AI_MARRVEL -r nextflow_conversion \
                  --ref_dir <PATH_TO_REFERENCE_DIRECTORY> \
                  --ref_ver [Reference genome: hg19/hg38] \
                  --input_vcf <PATH_TO_INPUT_VCF_FILE> \
                  --input_hpo <PATH_TO_INPUT_HPO_FILE> \
                  --outdir <PATH_TO_OUTPUT_DIRECTORY> \
                  --run_id [Sample Id] 



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

   Due to licensing restrictions, the HGMD-related database is not included in the public download. To prepare HGMD database, see below.


Prepare HGMD-related data
===================================

**HGMD databse VCF**

Copy your HGMD ``bgzip`` compressed VCF file and the index file (\*.tbi) into your database folder under ``vep/hg19(hg38)/``, and rename it as ``HGMD_Pro_2022.2_hg19(hg38).vcf.gz(.tbi)``.


**HGMD Phenotype Information**

Get the ``hgmd_phenbase-*.dump.sql`` from HGMD.

Run `HGMD_phenbase.sh <https://github.com/LiuzLab/AI_MARRVEL/blob/main/utils/HGMD_phenbase.sh>`_ to get ``HGMD_phen.tsv``. Move the file to your database folder under ``omim_annotate/hg19(hg38)/``
