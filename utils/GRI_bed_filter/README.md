# Genomic Region of Interest filter

\
\

##### Data files for generation:

1.  Gencode v46 gene annotation\
2.  HGMD 2022 database\
3.  Clinvar 8/6/2024 database\
4.  SpliceAI v1.3 prediction data\

##### BED file description:

1.  Version 1 --- gene only: BED \<--- protein_coding gene + HGMD mutations (50 flanking bp) + Cinvar mutations (50 flanking bp)\
    (Note: Promoter upstream 1kb is dropped as requested by Dr.Liu)\
2.  Version 2 --- exon only: BED \<--- protein_coding gene + HGMD mutations (50 flanking bp) + Cinvar mutations (50 flanking bp) + SpliceAI predicted positions (50 flanking bp)\
3.  Additionally, initial bed files which also contain a source column can be generated, that indicates where each position entry comes from (gencode, hgmd, clinvar, or spliceAI)

##### BED file format:

1.  chr \| start \| end\
2.  The 1st column contains prefix of "chr" for chromosomes, for example, "chr1, chr2, .... chr Y".