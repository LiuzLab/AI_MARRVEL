.. _quickstart:

*************
Quick Start
*************

AIM-Lite
=============

Install the software. And download our example here. 

.. code-block:: bash
    
    docker run -v <Path to VCF File>:/input/vcf.gz \
            -v <Path to HPO file>:/input/hpo.txt \
            -v <Path to downloaded database>:/run/data_dependencies \
            -v <Path to output folder>:/out \
        chaozhongliu/aim-lite /run/proc.sh [Sample ID] [Reference genome: hg19/hg38] [Memory Limit (G)]




