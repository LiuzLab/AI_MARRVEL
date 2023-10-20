#!/bin/bash
#$1 submission id
#$2 ref genome
#$1 memory limitation


#make output subdirectories
echo $(date +"%T") > /out/whole_log.txt
mkdir -m777 /out/phrank
mkdir -m777 /out/input
mkdir -m777 /out/tier-test-false
mkdir -m777 /out/rami-test
mkdir -m777 /out/scores
mkdir -m777 /out/final_matrix
mkdir -m777 /out/final_matrix_expanded

REF_DIR=$2
CHUNK_RAM=$1

if [ "$REF_DIR" = "hg19" ]
then
    REF_GEN="GRCh37"
else
    REF_GEN="GRCh38"
fi


echo "bcf annotate"
#annotate with new chromosomes and preserve original coordinates in ID
bcftools annotate --rename-chrs /run/data_dependencies/bcf_annotate/chrmap.txt -x ID /input/vcf.gz -Oz -o /out/$1-annot.txt
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' /out/$1-annot.txt -Oz -o /out/$1-add-id.vcf.gz


echo "quality filters"
#run quality filters #removed 27th line
bcftools filter /out/$1-add-id.vcf.gz -i'FILTER == "PASS"' -Oz -o /out/$1.filt.vcf.gz


echo "remove mitochondrial"
gunzip /out/$1.filt.vcf.gz
grep -vE "chrM" /out/$1.filt.vcf > /out/$1.filt.rmMT.vcf


echo "remove unknown chromosomes"
###############
# To-Do
###############


#Phrank annotation
echo "phrank annotation"
zcat /input/vcf.gz | awk 'substr($0, 1, 1) != "#"' | cut -f1,2,4,5 | sed 's/\t/:/g' > /out/$1-var.txt
cat /out/$1-var.txt | sed 's/chr//g' | sort -u > /out/$1-var-filt.txt


#annotate with ENSG
echo "ensg annotates"
#/opt/conda/envs/marrvel-py/bin/python $4/phrank/src/location_to_gene.py /out/$1-var-filt.txt ${REF_DIR} | sed 's/:/\t/g' | sed 's/X\t/23\t/g' | sed 's/Y\t/24\t/g' | sed 's/MT\t/25\t/g' > /out/$1-ensmbl.txt
python2.7 /run/phrank/src/location_to_gene.py /out/$1-var-filt.txt ${REF_DIR} | sed 's/:/\t/g' | sed 's/X\t/23\t/g' | sed 's/Y\t/24\t/g' | sed 's/MT\t/25\t/g' > /out/$1-ensmbl.txt


echo "gene-sym"
cat /out/$1-ensmbl.txt | sort -k5,5 | join -1 5 -2 1 - /run/data_dependencies/phrank/hg19/ensembl_to_symbol.txt | sed 's/ /\t/g' | cut -f2- > /out/$1-genesym.txt
cat /out/$1-genesym.txt | cut -f5 | sort -u | join -t$'\t' -1 1 -2 2 - /run/data_dependencies/phrank/hg19/gene_to_symbol_sorted.txt | cut -f2 | sort -u > /out/$1-gene.txt


# #run phrank scoring
echo "phrank scoring"
python3.8 /run/phrank/src/run_phrank.py /out/$1-gene.txt /input/hpo.txt ${REF_DIR} > /out/phrank/$1.txt


# #annotate with OMIM and HPO path
echo "annotate with OMIM and HPO path"
Rscript /run/chaozhong-annot.R /input/hpo.txt $1 ${REF_DIR} > /dev/null


#Filter proband vcf with blacklist
echo "Filter proband vcf with blacklist"
bgzip /out/$1.filt.rmMT.vcf

tabix /out/$1.filt.rmMT.vcf.gz
mkdir -m777 /out/isec_tmp1
bcftools isec -p /out/isec_tmp1 -w 1 -Oz \
/out/$1.filt.rmMT.vcf.gz /run/data_dependencies/filter_vep/${REF_DIR}/gnomad.${REF_DIR}.blacklist.genomes.vcf.gz

tabix /out/isec_tmp1/0000.vcf.gz
mkdir -m777 /out/isec_tmp2
bcftools isec -p /out/isec_tmp2 -w 1 -Oz \
/out/$1.filt.rmMT.vcf.gz /run/data_dependencies/filter_vep/${REF_DIR}/gnomad.${REF_DIR}.blacklist.exomes.vcf.gz

tabix /out/isec_tmp2/0000.vcf.gz
mkdir -m777 /out/isec_tmp3
bcftools isec -p /out/isec_tmp3 -Ov \
/out/isec_tmp1/0000.vcf.gz /out/isec_tmp2/0000.vcf.gz

mv /out/isec_tmp3/0002.vcf /out/$1.filt.rmBL.vcf
rm -rf /out/isec_tmp1
rm -rf /out/isec_tmp2
rm -rf /out/isec_tmp3


#annotate with vep
echo "VEP annotation"

if [ "$REF_DIR" = "hg38" ]
then
    
    /opt/vep/src/ensembl-vep/vep \
        --dir_cache /run/data_dependencies/vep/hg38/  \
        --dir_plugins /run/data_dependencies/vep/hg38/Plugins/ \
        --fork 10 --everything --format vcf \
        --cache --offline --tab --force_overwrite \
        --species homo_sapiens --assembly GRCh38 \
        --custom /run/data_dependencies/vep/hg38/gnomad.genomes.GRCh38.v3.1.2.sites.vcf.bgz,gnomADg,vcf,exact,0,AF,AF_popmax,controls_nhomalt \
        --custom /run/data_dependencies/vep/hg38/clinvar_20220730.vcf.gz,clinvar,vcf,exact,0,CLNREVSTAT,CLNSIG,CLNSIGCONF \
        --custom /run/data_dependencies/vep/hg38/HGMD_Pro_2022.2_hg38.vcf.gz,hgmd,vcf,exact,0,CLASS,GENE,PHEN,RANKSCORE \
        --af_gnomad  --plugin REVEL,/run/data_dependencies/vep/hg38/new_tabbed_revel_grch38.tsv.gz,ALL \
        --plugin SpliceAI,snv=/run/data_dependencies/vep/hg38/spliceai_scores.masked.snv.hg38.vcf.gz,indel=/run/data_dependencies/vep/hg38/spliceai_scores.masked.indel.hg38.vcf.gz,cutoff=0.5 \
        --plugin CADD,/run/data_dependencies/vep/hg38/hg38_whole_genome_SNV.tsv.gz,ALL \
        --plugin dbNSFP,/run/data_dependencies/vep/hg38/dbNSFP4.1a_grch38.gz,ALL \
        --individual all --output_file /out/$1-vep.txt --input_file /out/$1.filt.rmBL.vcf
else

     /opt/vep/src/ensembl-vep/vep \
        --dir_cache /run/data_dependencies/vep/hg19/ \
        --dir_plugins /run/data_dependencies/vep/hg19/Plugins/ \
        --fork 10 --everything --format vcf \
        --cache --offline --tab --force_overwrite \
        --species homo_sapiens  --assembly GRCh37 \
        --custom /run/data_dependencies/vep/hg19/gnomad.genomes.r2.1.sites.grch37_noVEP.vcf.gz,gnomADg,vcf,exact,0,AF,AF_popmax,controls_nhomalt \
        --custom /run/data_dependencies/vep/hg19/clinvar_20220730.vcf.gz,clinvar,vcf,exact,0,CLNREVSTAT,CLNSIG,CLNSIGCONF \
        --custom /run/data_dependencies/vep/hg19/HGMD_Pro_2022.2_hg19.vcf.gz,hgmd,vcf,exact,0,CLASS,GENE,PHEN,RANKSCORE \
        --af_gnomad  --plugin REVEL,/run/data_dependencies/vep/hg19/new_tabbed_revel.tsv.gz,ALL \
        --plugin SpliceAI,snv=/run/data_dependencies/vep/hg19/spliceai_scores.masked.snv.hg19.vcf.gz,indel=/run/data_dependencies/vep/hg19/spliceai_scores.masked.indel.hg19.vcf.gz,cutoff=0.5 \
        --plugin CADD,/run/data_dependencies/vep/hg19/hg19_whole_genome_SNVs.tsv.gz,ALL \
        --plugin dbNSFP,/run/data_dependencies/vep/hg19/dbNSFP4.3a_grch37.gz,ALL \
        --individual all --output_file /out/$1-vep.txt --input_file /out/$1.filt.rmBL.vcf
fi

#annotate w/ marrvel flatfiles
AIM_FREE_RAM=$(free -g | awk 'NR==2{printf $7}')

if (( AIM_FREE_RAM > CHUNK_RAM )); then
    RAM_MAX=${CHUNK_RAM}
else
    RAM_MAX=${AIM_FREE_RAM}
fi

python3.8 /run/split_chunks.py $1 $RAM_MAX

while read -r INDEX LINEH LINEA LINEB
do
    sed -n -e "${LINEH}p" -e "${LINEA},${LINEB}p" /out/$1-vep.txt > /out/$1-vep-${INDEX}.txt

    python3.8 /run/annotation/main.py \
        -outPrefix /out/rami-test/r1 \
        -patientID $1 \
        -patientHPOsimiOMIM /out/$1-dx \
        -patientHPOsimiHGMD /out/$1-cz \
        -varFile /out/$1-vep-$INDEX.txt \
        -inFileType vepAnnotTab \
        -patientFileType one \
        -genomeRef ${REF_DIR} \
        -diseaseInh AD \
        -modules curate,conserve
    
    if [ $INDEX -gt 1 ]; then
        sed -n "2,\$p" /out/rami-test/$1_scores.csv > /out/rami-test/$1_scores_$INDEX.csv
        rm /out/rami-test/$1_scores.csv
        sed -n "2,\$p" /out/rami-test/r1_$1_scores.txt > /out/rami-test/r1_$1_scores_$INDEX.txt
        rm /out/rami-test/r1_$1_scores.txt
    else
        mv /out/rami-test/$1_scores.csv /out/rami-test/$1_scores_$INDEX.csv
        mv /out/rami-test/r1_$1_scores.txt /out/rami-test/r1_$1_scores_$INDEX.txt
    fi
done < /out/vep_split.txt


# Combine all $1_scores_$INDEX.csv
for INDEX in $(cut -d$'\t' -f1 /out/vep_split.txt)
do
    cat /out/rami-test/$1_scores_$INDEX.csv
done > /out/rami-test/$1_scores.csv

for INDEX in $(cut -d$'\t' -f1 /out/vep_split.txt)
do
    cat /out/rami-test/r1_$1_scores_$INDEX.txt
done > /out/rami-test/r1_$1_scores.txt

rm /out/$1-vep-*.txt
rm /out/rami-test/r1_$1_scores_*.txt
rm /out/rami-test/$1_scores_*.csv

Rscript /run/VarTierDiseaseDBFalse.R $1 ${REF_DIR}

python3.8 /run/generate_new_matrix_2.py $1 ${REF_DIR}

python3.8 /run/predict_new/run_final.py $1

python3.8 /run/merge_rm.py $1

python3.8 /run/extraModel/main.py -id $1


echo $(date +"%T") >> /out/whole_log.txt


######################
# Delete any output?
######################










