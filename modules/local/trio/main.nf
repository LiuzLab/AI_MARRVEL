process VCF_PRE_PROCESS_TRIO {
    container "broadinstitute/gatk"
    publishDir "${params.outdir}/${params.run_id}/vcf_trio", mode: 'copy'

    input:
    path vcf
    path ped
    path fasta
    path fasta_index
    path fasta_dict
    path chrmap_file

    output:
    path "${params.run_id}_fixed.vcf.gz", emit: "vcf"
    path "${params.run_id}.inheritance.txt", emit: "inheritance"

    script:
    """
    # Extract proband, maternal, and paternal id from ped file, sample IDs in ped should be the same as in the vcf file

    # 1) Make sure there are exactly 3 non-empty lines
    row_count=\$(awk 'NF' "$ped" | wc -l)
    if [ "\$row_count" -ne 3 ]; then
        echo "ERROR: Expected 3 rows in '$ped' (found \$row_count)" >&2
        exit 1
    fi

    # 2) Find the proband row (both parent fields != 0)
    #    PED format: 1.FAMID 2.INDIV 3.PID 4.MID 5.SEX 6.PHENOTYPE
    proband_line=\$(awk '\$3 != "0" && \$4 != "0" { print; exit }' "$ped")
    if [ -z "\$proband_line" ]; then
        echo "ERROR: No row with both paternal (field 3) and maternal (field 4) IDs found" >&2
        exit 1
    fi

    # 3) Parse out the IDs
    #    We only care about fields 2,3,4 here.
    read -r _ SAMPLE_P SAMPLE_F SAMPLE_M _ <<< "\$proband_line"

    echo "Normalize joint VCF file, split multiallelics rows"
    bcftools norm --multiallelics -both \
                -Oz -o ./${params.run_id}.tmp.vcf.gz \
                $vcf
    tabix -p vcf ./${params.run_id}.tmp.vcf.gz

    bcftools norm --rm-dup none \
                -Oz -o ./${params.run_id}.joint.vcf.gz \
                ./${params.run_id}.tmp.vcf.gz

    tabix -p vcf ./${params.run_id}.joint.vcf.gz
    zcat ./${params.run_id}.joint.vcf.gz | grep '^#CHROM' > ./vcfheaders.txt    #* extract header from vcf file

    # create ./sample_ids.txt
    trio_rename_VCFheader.py \
        \$SAMPLE_P \
        \$SAMPLE_F \
        \$SAMPLE_M \
        ./vcfheaders.txt

    bcftools reheader --samples ./sample_ids.txt \
                  -o ./${params.run_id}.joint2.vcf.gz \
                  ./${params.run_id}.joint.vcf.gz

    mv ./${params.run_id}.joint2.vcf.gz ./${params.run_id}.joint.vcf.gz
    tabix -f -p vcf ./${params.run_id}.joint.vcf.gz

    bcftools annotate --rename-chrs "$chrmap_file" -x ID ./${params.run_id}.joint.vcf.gz -Oz -o step1.vcf.gz
    tabix step1.vcf.gz
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y step1.vcf.gz -Oz -o step1_1.vcf.gz

    mv step1_1.vcf.gz ./${params.run_id}.joint.vcf.gz

    # Index the final output using tabix
    echo "Indexing the final output with tabix..."
    tabix -f -p vcf ./${params.run_id}.joint.vcf.gz

    # De Novo Calling using VariantAnnotator
    #* GATK 4.2.6.0
    echo "GATK De Novo Calling..."
    gatk VariantAnnotator \
        -A PossibleDeNovo \
        -R ${fasta} \
        --pedigree $ped \
        -V ./${params.run_id}.joint.vcf.gz \
        -O ./${params.run_id}.DeNovo.g.vcf

    echo "Check variants inheritance..."

    # Extract patient VCF from joint VCF
    #=======================================================
    # TODOs
    #=======================================================

    echo "Extract proband VCF..."
    #* extract the proband variants from the joint vcf file
    bcftools view -c1 -Oz -o./${params.run_id}.vcf.gz -s \${SAMPLE_P} ${params.run_id}.joint.vcf.gz

    echo "Fix VCF formatting..."
    cat ${projectDir}/resources/trio_pipeline/acceptable_header.txt | grep -v '#CHROM' > header.txt
    zcat ./${params.run_id}.vcf.gz | grep '#CHROM' | head -1 >> header.txt
    zcat ./${params.run_id}.vcf.gz | grep -v '#' | awk '
        BEGIN {FS = "\\t"} {
            \$7 = "PASS";
            \$8 = ".";
            gsub(":.*", "", \$9);
            gsub(":.*", "", \$10);
            gsub("\\|", "/", \$10);
            print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8 "\\t" \$9 "\\t" \$10
        }' > body.txt
    cat header.txt body.txt | bgzip > ./${params.run_id}_fixed.vcf.gz

    # create ${params.run_id}.inheritance.txt and ${params.run_id}.summary.txt
    trio_check_inheritance.GATK.py \
        ${params.run_id} \
        ./${params.run_id}.DeNovo.g.vcf \
        ${params.ref_ver} \
        \$SAMPLE_P \
        \$SAMPLE_F \
        \$SAMPLE_M


    """
}

process GENERATE_TRIO_FEATURES {
    publishDir "${params.outdir}/${params.run_id}/features_trio", mode: 'copy'

    input:
    path merged_compressed_scores
    path merged_matrix
    path inheritance

    output:
    path "${params.run_id}.triomatrix.tsv", emit: triomatrix

    script:
    """
    trio_merge_inheritance.py \
        $merged_compressed_scores \
        $merged_matrix \
        $inheritance \
        ./${params.run_id}.triomatrix.tsv
    """
}

process PREDICTION_TRIO {
    container "zhandongliulab/aim-lite-oldpython"
    publishDir "${params.outdir}/${params.run_id}/prediction_trio", mode: 'copy'

    input:
    path merged_compressed_scores
    path triomatrix

    output:
    path "*.trio.expanded.csv.gz"
    path "*.trio.prediction.csv"

    script:
    """
    run_final.py \
        ${projectDir}/resources/trio_pipeline/rf_trio_base.job \
        ${projectDir}/resources/trio_pipeline/features.csv \
        $triomatrix \
        ./${params.run_id}.trio.csv
    run_final.py \
        ${projectDir}/resources/trio_pipeline/rf_trio_NDG.job \
        ${projectDir}/resources/trio_pipeline/features_NDG.csv \
        $triomatrix \
        ./${params.run_id}.trio.NDG.csv


    # Generate ${params.run_id}.trio.prediction.csv
    trio_merge_rm.py ${params.run_id}

    merge_rm.py \
        ${params.run_id}.trio.prediction.csv \
        $merged_compressed_scores \
        ${params.run_id}.trio.expanded.csv.gz
    """
}
