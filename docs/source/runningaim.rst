.. _runningaim:

*************
Running AIM
*************

Quick Start
=============

**AIM-Lite**

Install the software and data dependencies. Download our example here. 

.. code-block:: bash
    
    docker run -u $(id -u):$(id -g) \
            -v <Path to VCF File>:/input/vcf.gz \
            -v <Path to HPO file>:/input/hpo.txt \
            -v <Path to downloaded database>:/run/data_dependencies \
            -v <Path to output folder>:/out \
        chaozhongliu/aim-lite /run/proc.sh [Sample ID] [Reference genome: hg19/hg38] [Memory Limit (G)]



Arguments Explanation
=======================



