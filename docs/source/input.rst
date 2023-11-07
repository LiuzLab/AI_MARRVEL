.. _input:

**************************
Input Data Preparation
**************************

Input VCF file
===================
AIM accepts most of the variant calling VCF output. Please see below the checklist before running the software.

#. The VCF file **must** contain FILTER column with PASS flag
#. The VCF **should not** contain any unknown chromosomes (future version of AIM plans to handle this automatically)
#. Both hg19 and hg38 referene genomes are acceptable.
#. The VCF file **must** be compressed using ``bgzip``


Input HPO file
===================
AIM accepts phenotype information as a list of HPO (Human Phenotype Ontology) terms.

Here is an example file called ``hpos.txt``::

    HP:0000365
    HP:0000508
    HP:0001249
    HP:0001263
    HP:0001999

If you only have clinical notes in hand, please try our `AI-MARRVEL website <https://ai.marrvel.org/>`_ to convert it into HPO terms.


.. note::

    We are packaging this ``clinical notes to HPO terms`` step into AIM pipeline.

