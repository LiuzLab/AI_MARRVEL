# #!/bin/bash
# #$1 submission id
# #$2 ref genome
# #$3 memory limitation

# #make output subdirectories
# echo $(date +"%T") > /out/whole_log.txt
# mkdir -m777 /out/phrank
# mkdir -m777 /out/input
# mkdir -m777 /out/tier-test-false
# mkdir -m777 /out/rami-test
# mkdir -m777 /out/scores
# mkdir -m777 /out/final_matrix
# mkdir -m777 /out/final_matrix_expanded

# REF_DIR=$2
# CHUNK_RAM=$3
# KEEP_INTERMEDIATE=${4:-True}

# if [ "$REF_DIR" = "hg19" ]
# then
#     REF_GEN="GRCh37"
# else
#     REF_GEN="GRCh38"
# fi

# #Check if hpo.txt is empty or has no valid HPO ID
# if [[ -z $(egrep 'HP:[0-9]{7}' /input/hpo.txt) ]] ; then
#     #Replace it with HP:0000001
#     echo "HP:0000001" > /input/hpo.txt
# fi


# cp /input/vcf.gz /out/vcf.gz  ### New Step: Copying file to output directory


# # Initial path of the VCF file
# INPUT_VCF_PATH="/out/vcf.gz"  ### Changed: Adjust path to new output location

# #Check if the file exists
# if [ -e "$INPUT_VCF_PATH" ]; then
#     #Use file command to determine the file type
#     INPUT_VCF_TYPE=$(file -b "$INPUT_VCF_PATH")

#     # Determine if the file is gzip or bgzip format
#     if echo "$INPUT_VCF_TYPE" | grep -q 'gzip compressed data'; then
#         # Check if it's bgzip format for tabix compatibility
#         if echo "$INPUT_VCF_TYPE" | grep -q 'BGZF'; then
#             echo "The file is in BGZF format, ready for tabix."
#             tabix -p vcf "$INPUT_VCF_PATH"  ### Adjusted Step: Only index if BGZF
#         else
#             echo "GZIP format detected, converting to BGZF."
#             gunzip -c "$INPUT_VCF_PATH" | bgzip > "${INPUT_VCF_PATH%.gz}"
#             mv "${INPUT_VCF_PATH%.gz}" "$INPUT_VCF_PATH"  ### Ensuring the correct .gz extension
#             tabix -p vcf "$INPUT_VCF_PATH"  ### Index the newly bgzipped file
#         fi
#     elif echo "$INPUT_VCF_TYPE" | grep -q 'ASCII text'; then
#         echo "Plain VCF file detected, compressing and indexing."
#         bgzip -c "$INPUT_VCF_PATH" > "${INPUT_VCF_PATH}.gz"
#         mv "${INPUT_VCF_PATH}.gz" "$INPUT_VCF_PATH"  ### Replace uncompressed with compressed
#         tabix -p vcf "$INPUT_VCF_PATH"  ### Index the newly bgzipped file
#     fi
# else
#     echo "The file $INPUT_VCF_PATH does not exist."
# fi



# echo "VCF pre-processing"
# #annotate with new chromosomes and preserve original coordinates in ID
# bcftools annotate --rename-chrs /run/data_dependencies/bcf_annotate/chrmap.txt -x ID /out/vcf.gz -Oz -o /out/$1-annot.txt
# bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' /out/$1-annot.txt -Oz -o /out/$1-add-id.vcf.gz


# echo "VCF quality filtering"
# #run quality filters #removed 27th line
# bcftools filter /out/$1-add-id.vcf.gz -i'FILTER == "PASS"' -Oz -o /out/$1.filt.vcf.gz

# #check number of variants left
# variant_count=$(bcftools view -H /out/$1.filt.vcf.gz | wc -l)
# if [ "$variant_count" -gt 0 ]; then
#     echo "Quality filtering completed successfully. Variants passing the filters: $variant_count"
# else
#     echo "The VCF file doesn't have FILTER annotation, or all variants filtered."
#     echo "Pipeline will proceed with unfiltered VCF file."
#     cp /out/$1-add-id.vcf.gz /out/$1.filt.vcf.gz
# fi

# echo "Remove mitochondrial and unknown chromosome variants"
# #gunzip /out/$1.filt.vcf.gz
# #grep -vE "chrM" /out/$1.filt.vcf > /out/$1.filt.rmMT.vcf
# tabix -p vcf /out/$1.filt.vcf.gz
# bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y /out/$1.filt.vcf.gz -o /out/$1.filt.rmMT.vcf

# #Phrank annotation
# echo "Phrank scoring"
# zcat /out/vcf.gz | awk 'substr($0, 1, 1) != "#"' | cut -f1,2,4,5 | sed 's/\t/:/g' > /out/$1-var.txt
# cat /out/$1-var.txt | sed 's/chr//g' | sort -u > /out/$1-var-filt.txt


# #annotate with ENSG
# #echo "ensg annotates"
# python2.7 /run/phrank/src/location_to_gene.py /out/$1-var-filt.txt ${REF_DIR} | sed 's/:/\t/g' | sed 's/X\t/23\t/g' | sed 's/Y\t/24\t/g' | sed 's/MT\t/25\t/g' > /out/$1-ensmbl.txt


# #echo "gene-sym"
# cat /out/$1-ensmbl.txt | sort -k5,5 | join -1 5 -2 1 - /run/data_dependencies/phrank/hg19/ensembl_to_symbol.txt | sed 's/ /\t/g' | cut -f2- > /out/$1-genesym.txt
# cat /out/$1-genesym.txt | cut -f5 | sort -u | join -t$'\t' -1 1 -2 2 - /run/data_dependencies/phrank/hg19/gene_to_symbol_sorted.txt | cut -f2 | sort -u > /out/$1-gene.txt


# # #run phrank scoring
# #echo "phrank scoring"
# python3.8 /run/phrank/src/run_phrank.py /out/$1-gene.txt /input/hpo.txt ${REF_DIR} > /out/phrank/$1.txt


# # #annotate with OMIM and HPO path
# echo "Calculating phenotype similarities"
# Rscript /run/phenoSim.R /input/hpo.txt $1 ${REF_DIR} > /dev/null


# #Filter proband vcf with blacklist
# echo "Filtering VCF with blacklist"
# bgzip /out/$1.filt.rmMT.vcf

# tabix /out/$1.filt.rmMT.vcf.gz
# mkdir -m777 /out/isec_tmp1
# bcftools isec -p /out/isec_tmp1 -w 1 -Oz \
# /out/$1.filt.rmMT.vcf.gz /run/data_dependencies/filter_vep/${REF_DIR}/gnomad.${REF_DIR}.blacklist.genomes.vcf.gz

# tabix /out/isec_tmp1/0000.vcf.gz
# mkdir -m777 /out/isec_tmp2
# bcftools isec -p /out/isec_tmp2 -w 1 -Oz \
# /out/$1.filt.rmMT.vcf.gz /run/data_dependencies/filter_vep/${REF_DIR}/gnomad.${REF_DIR}.blacklist.exomes.vcf.gz

# tabix /out/isec_tmp2/0000.vcf.gz
# mkdir -m777 /out/isec_tmp3
# bcftools isec -p /out/isec_tmp3 -Ov \
# /out/isec_tmp1/0000.vcf.gz /out/isec_tmp2/0000.vcf.gz

# mv /out/isec_tmp3/0002.vcf /out/$1.filt.rmBL.vcf
# rm -rf /out/isec_tmp1
# rm -rf /out/isec_tmp2
# rm -rf /out/isec_tmp3


# #annotate with vep
# echo "VEP annotation"

# if [ "$REF_DIR" = "hg38" ]
# then
    
#     /opt/vep/src/ensembl-vep/vep \
#         --dir_cache /run/data_dependencies/vep/hg38/  \
#         --dir_plugins /run/data_dependencies/vep/hg38/Plugins/ \
#         --fork 10 --everything --format vcf \
#         --cache --offline --tab --force_overwrite \
#         --species homo_sapiens --assembly GRCh38 \
#         --custom /run/data_dependencies/vep/hg38/gnomad.genomes.GRCh38.v3.1.2.sites.vcf.bgz,gnomADg,vcf,exact,0,AF,AF_popmax,controls_nhomalt \
#         --custom /run/data_dependencies/vep/hg38/clinvar_20220730.vcf.gz,clinvar,vcf,exact,0,CLNREVSTAT,CLNSIG,CLNSIGCONF \
#         --custom /run/data_dependencies/vep/hg38/HGMD_Pro_2022.2_hg38.vcf.gz,hgmd,vcf,exact,0,CLASS,GENE,PHEN,RANKSCORE \
#         --af_gnomad  --plugin REVEL,/run/data_dependencies/vep/hg38/new_tabbed_revel_grch38.tsv.gz,ALL \
#         --plugin SpliceAI,snv=/run/data_dependencies/vep/hg38/spliceai_scores.masked.snv.hg38.vcf.gz,indel=/run/data_dependencies/vep/hg38/spliceai_scores.masked.indel.hg38.vcf.gz,cutoff=0.5 \
#         --plugin CADD,/run/data_dependencies/vep/hg38/hg38_whole_genome_SNV.tsv.gz,ALL \
#         --plugin dbNSFP,/run/data_dependencies/vep/hg38/dbNSFP4.1a_grch38.gz,ALL \
#         --individual all --output_file /out/$1-vep.txt --input_file /out/$1.filt.rmBL.vcf
# else

#      /opt/vep/src/ensembl-vep/vep \
#         --dir_cache /run/data_dependencies/vep/hg19/ \
#         --dir_plugins /run/data_dependencies/vep/hg19/Plugins/ \
#         --fork 10 --everything --format vcf \
#         --cache --offline --tab --force_overwrite \
#         --species homo_sapiens  --assembly GRCh37 \
#         --custom /run/data_dependencies/vep/hg19/gnomad.genomes.r2.1.sites.grch37_noVEP.vcf.gz,gnomADg,vcf,exact,0,AF,AF_popmax,controls_nhomalt \
#         --custom /run/data_dependencies/vep/hg19/clinvar_20220730.vcf.gz,clinvar,vcf,exact,0,CLNREVSTAT,CLNSIG,CLNSIGCONF \
#         --custom /run/data_dependencies/vep/hg19/HGMD_Pro_2022.2_hg19.vcf.gz,hgmd,vcf,exact,0,CLASS,GENE,PHEN,RANKSCORE \
#         --af_gnomad  --plugin REVEL,/run/data_dependencies/vep/hg19/new_tabbed_revel.tsv.gz,ALL \
#         --plugin SpliceAI,snv=/run/data_dependencies/vep/hg19/spliceai_scores.masked.snv.hg19.vcf.gz,indel=/run/data_dependencies/vep/hg19/spliceai_scores.masked.indel.hg19.vcf.gz,cutoff=0.5 \
#         --plugin CADD,/run/data_dependencies/vep/hg19/hg19_whole_genome_SNVs.tsv.gz,ALL \
#         --plugin dbNSFP,/run/data_dependencies/vep/hg19/dbNSFP4.3a_grch37.gz,ALL \
#         --individual all --output_file /out/$1-vep.txt --input_file /out/$1.filt.rmBL.vcf
# fi

# echo "Feature engineering"
# #annotate w/ marrvel flatfiles
# AIM_FREE_RAM=$(free -g | awk 'NR==2{printf $7}')

# if (( AIM_FREE_RAM > CHUNK_RAM )); then
#     RAM_MAX=${CHUNK_RAM}
# else
#     RAM_MAX=${AIM_FREE_RAM}
# fi

# python3.8 /run/split_chunks.py $1 $RAM_MAX

# while read -r INDEX LINEH LINEA LINEB
# do
#     sed -n -e "${LINEH}p" -e "${LINEA},${LINEB}p" /out/$1-vep.txt > /out/$1-vep-${INDEX}.txt

#     python3.8 /run/annotation/main.py \
#         -outPrefix /out/rami-test/r1 \
#         -patientID $1 \
#         -patientHPOsimiOMIM /out/$1-dx \
#         -patientHPOsimiHGMD /out/$1-cz \
#         -varFile /out/$1-vep-$INDEX.txt \
#         -inFileType vepAnnotTab \
#         -patientFileType one \
#         -genomeRef ${REF_DIR} \
#         -diseaseInh AD \
#         -modules curate,conserve
    
#     if [ $INDEX -gt 1 ]; then
#         sed -n "2,\$p" /out/rami-test/$1_scores.csv > /out/rami-test/$1_scores_$INDEX.csv
#         rm /out/rami-test/$1_scores.csv
#         sed -n "2,\$p" /out/rami-test/r1_$1_scores.txt > /out/rami-test/r1_$1_scores_$INDEX.txt
#         rm /out/rami-test/r1_$1_scores.txt
#     else
#         mv /out/rami-test/$1_scores.csv /out/rami-test/$1_scores_$INDEX.csv
#         mv /out/rami-test/r1_$1_scores.txt /out/rami-test/r1_$1_scores_$INDEX.txt
#     fi
# done < /out/vep_split.txt


# # Combine all $1_scores_$INDEX.csv
# for INDEX in $(cut -d$'\t' -f1 /out/vep_split.txt)
# do
#     cat /out/rami-test/$1_scores_$INDEX.csv
# done > /out/rami-test/$1_scores.csv

# for INDEX in $(cut -d$'\t' -f1 /out/vep_split.txt)
# do
#     cat /out/rami-test/r1_$1_scores_$INDEX.txt
# done > /out/rami-test/r1_$1_scores.txt

# rm /out/$1-vep-*.txt
# rm /out/rami-test/r1_$1_scores_*.txt
# rm /out/rami-test/$1_scores_*.csv

# Rscript /run/VarTierDiseaseDBFalse.R $1 ${REF_DIR}

# python3.8 /run/generate_new_matrix_2.py $1 ${REF_DIR}

# echo "Making prediction"
# python3.8 /run/predict_new/run_final.py $1

# python3.8 /run/merge_rm.py $1

# python3.8 /run/extraModel/main.py -id $1


# echo $(date +"%T") >> /out/whole_log.txt

# # delete intermediate files
# if [ "$KEEP_INTERMEDIATE" = "False" ]; then
#     rm /out/$1-add-id.vcf.gz
#     rm /out/$1-annot.txt
#     rm /out/$1-cz
#     rm /out/$1-dx
#     rm /out/$1-ensmbl.txt
#     rm /out/$1.filt.rmBL.vcf
#     rm /out/$1.filt.rmMT.vcf.gz
#     rm /out/$1.filt.rmMT.vcf.gz.tbi
#     rm /out/$1.filt.vcf
#     rm /out/$1-genesym.txt
#     rm /out/$1-gene.txt
#     rm /out/$1.matrix.txt
#     rm /out/$1-var-filt.txt
#     rm /out/$1-var.txt
#     rm /out/$1-vep.txt
#     rm /out/$1-vep.txt*
#     rm /out/vep_split.txt
#     rm -r /out/final_matrix
#     rm -r /out/final_matrix_expanded
#     rm -r /out/input
#     rm -r /out/phrank
#     rm -r /out/rami-test
#     rm -r /out/scores
#     rm -r /out/tier-test-false
# fi

# mv /out/conf_4Model/*.csv /out/
# mv /out/conf_4Model/integrated/*.csv /out/
# rm -r /out/conf_4Model

