process PREPARE_TRIO_VCF {
    container "broadinstitute/gatk"

    input:
    path vcf
    path chrmap_file
    path sample_ids

    output:
    path "${params.run_id}.joint.vcf.gz", emit: vcf
    path "${params.run_id}.joint.vcf.gz.tbi", emit: tbi

    script:
    """
    set -euo pipefail

    source "${sample_ids}"

    echo "Normalize joint VCF file, split multiallelic rows..."
    bcftools norm --multiallelics -both \
        -Oz -o "${params.run_id}.tmp.vcf.gz" \
        "$vcf"
    tabix -p vcf "${params.run_id}.tmp.vcf.gz"

    echo "Remove duplicate records..."
    bcftools norm --rm-dup none \
        -Oz -o "${params.run_id}.joint.vcf.gz" \
        "${params.run_id}.tmp.vcf.gz"

    tabix -p vcf "${params.run_id}.joint.vcf.gz"
    zcat "${params.run_id}.joint.vcf.gz" | grep '^#CHROM' > vcfheaders.txt

    echo "Rename trio samples in header..."
    trio_rename_VCFheader.py \
        "\${SAMPLE_P}" \
        "\${SAMPLE_F}" \
        "\${SAMPLE_M}" \
        vcfheaders.txt

    bcftools reheader --samples sample_ids.txt \
        -o "${params.run_id}.joint2.vcf.gz" \
        "${params.run_id}.joint.vcf.gz"

    mv "${params.run_id}.joint2.vcf.gz" "${params.run_id}.joint.vcf.gz"
    tabix -f -p vcf "${params.run_id}.joint.vcf.gz"

    echo "Rename chromosomes and filter contigs..."
    bcftools annotate --rename-chrs "$chrmap_file" -x ID "${params.run_id}.joint.vcf.gz" -Oz -o step1.vcf.gz
    tabix step1.vcf.gz
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y \
        step1.vcf.gz -Oz -o step1_1.vcf.gz

    mv step1_1.vcf.gz "${params.run_id}.joint.vcf.gz"
    tabix -f -p vcf "${params.run_id}.joint.vcf.gz"

    rm -f "${params.run_id}.tmp.vcf.gz" "${params.run_id}.tmp.vcf.gz.tbi" step1.vcf.gz* vcfheaders.txt sample_ids.txt
    """
}
