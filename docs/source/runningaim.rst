.. _runningaim:

*************
Running AIM
*************

Quick Start
=============

**AIM-Lite**

Install the software and data dependencies. Download our example `here <https://app.globus.org/file-manager?origin_id=bacf3c02-a7a5-4a4d-b706-43cf37f0445a&origin_path=%2F>`_. 

.. code-block:: bash
    
    docker run -u $(id -u):$(id -g) \
            -v <Path to VCF File>:/input/vcf.gz \
            -v <Path to HPO file>:/input/hpo.txt \
            -v <Path to downloaded database>:/run/data_dependencies \
            -v <Path to output folder>:/out \
        chaozhongliu/aim-lite /run/proc.sh [Sample ID] [Reference genome: hg19/hg38] [Memory Limit (G)]



Arguments Explanation
=======================

.. list-table:: AIM arguments
   :header-rows: 1
   :widths: 40 60
   
   *  -  Argument
      -  Explanation
   *  -  Path to VCF File
      -  Full path to the sample VCF. See Input Data Preparation for details 
   *  -  Path to HPO file
      -  Full path to HPO list, formatted as txt, each row is an HPO term.
   *  -  Path to downloaded database
      -  Full path to data dependencies downloaded
   *  -  Path to output folder
      -  Full path to output folder, must already exist
   *  -  Sample ID
      -  Can be any ID. It's used as output frefix 
   *  -  Reference genome
      -  Reference genome used to map reads and call variants. AIM support hg19 and hg38
   *  -  Memory Limit
      -  The maximum memory allocated to AIM. If the system has 48G RAM, 32 to 36 as recommended
   


