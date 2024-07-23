# Trio case prediction pipeline

#$1 submission id
#$2 ref genome
#$3 memory limitation
#$4 keep intermediate

#input joint vcf file mount as---> input/joint.vcf
#input hpo file mount as---> input/hpo.txt
#input ped file mount as---> input/trio.ped
#output folder mount as---> out/
#data dependency folder mount as---> /run/data_dependencies


# Make output subdirectories
echo $(date +"%T") > /out/whole_log_trio.txt
mkdir -m777 /out/inheritance

# Extract proband, maternal, and paternal id from ped file, sample IDs in ped should be the same as in the vcf file
lastrow=$(tail -n 1 "input/trio.ped")
SAMPLE_P=$(echo "$lastrow" | awk '{print $2}')
SAMPLE_F=$(echo "$lastrow" | awk '{print $3}')
SAMPLE_M=$(echo "$lastrow" | awk '{print $4}')

# GATK for De Novo Calling
#=======================================================
# gatk software
# dependence files
# get ref genome version
# father mother patient ID extraction
# python docker env envifinne
#=======================================================

echo "Normalize joint VCF file, split multiallelics rows"
bcftools norm --multiallelics -both \
              -Oz -o /out/${1}.tmp.vcf.gz \
              input/joint.vcf

tabix out/${1}.tmp.vcf.gz   #* index the vcf file


#? is this a redundant step? since we already normalized the vcf file
bcftools norm --rm-dup none \
              -Oz -o /out/${1}.joint.vcf.gz \
              /out/${1}.tmp.vcf.gz

tabix out/${1}.joint.vcf.gz  #* index the vcf file
rm -f out/${1}.tmp.vcf.gz
rm out/${1}.tmp.vcf.gz.tbi

#! It could be helpful if there are more than 3 samples in the vcf file, and it will only keep the targeted trio samples
zcat out/${1}.joint.vcf.gz | grep '^#CHROM' > out/vcfheaders.txt    #* extract header from vcf file


#* python3.8
python /trio_pipeline/rename_VCFheader.py ${SAMPLE_P} ${SAMPLE_F} ${SAMPLE_M} /out/vcfheaders.txt

bcftools reheader --samples /out/sample_ids.txt \   #sample_ids.txt from rename_VCFheader.py
                  -Oz \
                  -o /out/${1}.joint2.vcf.gz \
                  /out/${1}.joint.vcf.gz

mv -f out/${1}.joint2.vcf.gz out/${1}.joint.vcf.gz
rm out/${1}.joint.vcf.gz.tbi    #* update the index file
tabix out/${1}.joint.vcf.gz

# De Novo Calling using VariantAnnotator
echo "GATK De Novo Calling..."
    gatk VariantAnnotator \         #* GATK 4.2.6.0
    -A PossibleDeNovo \
    -R /run/data_dependencies/var_annotator/${2}/${2}.fa \
    --pedigree /input/trio.ped \
    -V /out/${1}.joint.vcf.gz \
    -O /out/inheritance/${1}.DeNovo.g.vcf

#rm $6/inheritance/$3.vcf.gz
#rm $6/inheritance/$3.vcf.gz.tbi


# *generate inheritance summary and txt files in the inheritance folder
echo "Check variants inheritance..."
python trio_pipeline/check_inheritance.GATK.py ${1} /out/inheritance/${1}.DeNovo.g.vcf ${2} ${SAMPLE_P} ${SAMPLE_F} ${SAMPLE_M}

# Extract patient VCF from joint VCF
#=======================================================
# TODOs
#=======================================================
echo "Extract proband VCF..."
bcftools view -c1 -Oz -o /out/${1}.vcf.gz -s ${SAMPLE_P} /out/${1}.joint.vcf.gz #* extract the proband variants from the joint vcf file

echo "Fix VCF formatting..."
(cat trio_pipeline/acceptable_header.txt | grep -v '#CHROM'; zcat out/${1}.vcf.gz | grep '#CHROM' | head -1; zcat out/${3}.vcf.gz | grep -v '#' | awk 'BEGIN {FS = "\t"} {$7 = "PASS"; $8 = "."; gsub(":.*", "", $9); gsub(":.*", "", $10); gsub("\\|", "/", $10); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}') | gzip > out/${1}_fixed.vcf.gz

# Run marrvel singleton pipeline
echo "Run AIM default pipeline..."
cp out/${1}_fixed.vcf.gz input/vcf.gz #* make input file available for singleton script
bash /run/proc.sh ${1} ${2} ${3} ${4} #* run the singleton pipeline for proband


# Prepare trio feature matrix
echo "Prepare trio feature matrix"
python /trio_pipeline/merge_inheritance.py \
        /out/rami-test/${1}_scores.csv \
        /out/final_matrix/${1}.csv \
        /out/inheritance/${1}.inheritance.txt \
        /out/inheritance/${1}.trio.mtx.csv

# Make prediction
#=======================================================
# change model in predict_trio
# change features in predict_trio
#=======================================================
echo "Run trio model to predict risk score"
#docker run  -v $6:/out -v $4:/run \
#-v /houston_30t/chaozhong/aimarrvel_pipeline/trio_pipeline/predict_trio:/run/predict_trio \
#lucianli123/marrvel-py python /run/predict_trio/run_final.py $3

# docker run -u $(id -u):$(id -g) -v $6:/out -v $4:/run \
# -v $4/predict_trio:/run/predict_trio \
# chaozhongliu/aimarrvel-trio-python python3.8 /run/predict_trio/run_final.py $3

python trio_pipeline/predict_trio/run_final.py ${1}

# docker run -u $(id -u):$(id -g) -v $6:/out -v $4:/run \
# -v $4/predict_trio:/run/predict_trio \
# chaozhongliu/aimarrvel-trio-python python3.8 /run/predict_trio/run_final_NDG.py $3

python  trio_pipeline/predict_trio/run_final_NDG.py ${1} 

# docker run -u $(id -u):$(id -g) \
# -v $6:/out -v $4:/run \
# -v $4:/run/inheritance \
# lucianli123/marrvel-py \
# python /run/inheritance/merge_rm_trio.py $3

python trio_pipeline/merged_rm_trio.py ${1}

