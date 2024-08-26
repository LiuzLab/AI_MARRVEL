.. _runningaim:

*************
Running AIM
*************

Quick Start
=============

AIM
=========

Install the software and data dependencies. Download our example `here <https://app.globus.org/file-manager?origin_id=bacf3c02-a7a5-4a4d-b706-43cf37f0445a&origin_path=%2F>`_.

.. code-block:: bash
    
   nextflow run Liuzlab/AI_MARRVEL -r nextflow_conversion \
                     --ref_dir <PATH_TO_REFERENCE_DIRECTORY> \
                     --ref_ver [Reference genome: hg19/hg38] \
                     --input_vcf <PATH_TO_INPUT_VCF_FILE> \
                     --input_hpo <PATH_TO_INPUT_HPO_FILE> \
                     --outdir <PATH_TO_OUTPUT_DIRECTORY> \
                     --run_id [Sample Id]     



Parameter Explanation
=======================

.. list-table:: AIM parameters
   :header-rows: 1
   :widths: 40 60
   
   *  -  Argument
      -  Explanation
   *  -  ``input_vcf``
      -  Full path to the sample VCF file, must be compressed (e.g., sample.vcf.gz).
   *  -  ``input_hpo``
      -  Full path to the HPO list, formatted as a txt file with one HPO term per row.
   *  -  ``ref_dir``
      -  Full path to the AIM pipeline dependencies directory.
   *  - ``out_dir``
      -  Full path to the output folder. The directory must already exist.
   *  -  ``run_id``
      -  Unique identifier for this run, used as the output file prefix.
   *  -  ``ref_ver``
      -  Specifies the reference genome version to use for variant analysis. Supported values are hg19 and hg38 (default: hg19).

For detailed information about each process, please refer to the documentation: https://ai-marrvel.readthedocs.io/en/latest/