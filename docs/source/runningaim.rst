.. _runningaim:

*************
Running AIM
*************

Quick Start
=============

**AIM**

Install the software and data dependencies. Download our example `here <https://app.globus.org/file-manager?origin_id=bacf3c02-a7a5-4a4d-b706-43cf37f0445a&origin_path=%2F>`_. 

.. code-block:: bash

    nextflow run /path/to/AI_MARRVEL/main.nf \
        -profile <debug/docker/singularity/> \
        --ref_dir /path/to/dependencies/ \
        --input_vcf /path/to/input.vcf \
        --input_hpo /path/to/input.hpos.txt \
        --outdir /path/to/output \
        --storedir /path/to/store \
        --run_id <Sample_ID> \
        --ref_ver <hg19/hg38>

Command Line Arguments
=======================

**General Options**
--------------------

*  ``--help``  
   [boolean, string] Show the help message for all top-level parameters.  
   When a parameter is given to ``--help``, the full help message of that parameter will be printed.

*  ``--help_full``  
   [boolean] Show the help message for all non-hidden parameters.

*  ``--show_hidden``  
   [boolean] Show all hidden parameters in the help message.  
   Must be used with ``--help`` or ``--help_full``.

**Input Options**
------------------

*  ``--input_vcf``  
   [string] Path to input VCF file.

*  ``--input_phenopacket``  
   [string] Path to input phenopacket JSON file.

*  ``--input_variant``  
   [string] Input variant in the format ``chrom-pos-ref-alt``, e.g., ``1-123456-A-T``.

*  ``--input_hpo``  
   [string] Path to input HPO file.

*  ``--input_ped``  
   [string] Path to the pedigree file for trio mode.

**Output Options**
-------------------

*  ``--outdir``  
   [string] Output directory.  
   *Default:* ``./out``

*  ``--storedir``  
   [string] Path to a shared directory for intermediate files.  
   Reusing this directory across runs skips regenerating common files and cuts runtime.  
   *Default:* ``./store``

**Reference and Run Options**
------------------------------

*  ``--ref_dir``  
   [string] Path to AIM pipeline dependencies directory.

*  ``--run_id``  
   [string] Unique identifier for this run.  
   *Default:* ``1``

*  ``--ref_ver``  
   [string] Reference genome version (accepted: ``hg19``, ``hg38``).  
   *Default:* ``hg19``

**Filtering Options**
----------------------

*  ``--bed_filter``  
   [string] Path to a BED file to perform analysis only for regions of interest.

*  ``--exome_filter``  
   [boolean] Enable exonic filter.  
   Additionally filters out variants outside of exonic regions.

*  ``--impact_filter``  
   [boolean] Enable low impact filter.  
   Additionally enables low impact transcripts.
