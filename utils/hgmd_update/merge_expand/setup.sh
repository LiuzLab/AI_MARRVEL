#!/bin/bash

#----------------------------- Setup -----------------------------

HGMD_hg19_DIR=HGMD_Pro_2024.3_hg19.vcf.gz
HGMD_hg38_DIR=HGMD_Pro_2024.3_hg38.vcf.gz

OUT_hg19_DIR=output/hg19
OUT_hg38_DIR=output/hg38

REF_hg19_DIR=references/hg19
REF_hg38_DIR=references/hg38
REF_aa=references/AA_abbreviations.txt

RSCRIPT_PATH=expandvar.R

# download reference genomes
wget -P $REF_hg19_DIR https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
wget -P $REF_hg38_DIR https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# decompress reference genomes
gunzip $REF_hg19_DIR/*.gz

# download AA abbreviations
# just Amino Acids 3 letters and 1 letter names, there is a copy at skyriver: /home/zhijiany/workdir/references/AA_abbreviations.txt

# setup conda environment
conda create -n hgmd r-base python=3.11 samtools -y
conda activate hgmd

# index reference genomes
conda run -n hgmd samtools faidx $REF_hg19_DIR/hg19.fa
conda run -n hgmd samtools faidx $REF_hg38_DIR/hg38.fa

# install transvar
conda run -n hgmd pip install transvar

# setup transvar annotations and reference genomes
#hg19
transvar config --download_anno --refversion hg19
transvar config -k reference -v $REF_hg19_DIR/hg19.fa --refversion hg19

#hg38
transvar config --download_anno --refversion hg38
transvar config -k reference -v $REF_hg38_DIR/hg38.fa --refversion hg38



#----------------------------- Run the R script -----------------------------

#Run the Rscript for hg19
conda run -n hgmd Rscript $RSCRIPT_PATH $HGMD_hg19_DIR $OUT_hg19_DIR $REF_aa hg19

#Run the R script for hg38
conda run -n hgmd Rscript $RSCRIPT_PATH $HGMD_hg38_DIR $OUT_hg38_DIR $REF_aa hg38

