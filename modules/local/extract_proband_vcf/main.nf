process EXTRACT_PROBAND_VCF {
    container "broadinstitute/gatk"

    input:
    path joint_vcf
    path sample_ids

    output:
    path "${params.run_id}_fixed.vcf.gz", emit: vcf

    script:
    """
    set -euo pipefail

    source "${sample_ids}"

    echo "Extract proband VCF..."
    bcftools view -c1 -Oz -o "${params.run_id}.vcf.gz" -s "\${SAMPLE_P}" "${joint_vcf}"

    echo "Fix VCF formatting..."
    grep -v '#CHROM' "${projectDir}/resources/trio_pipeline/acceptable_header.txt" > header.txt
    zcat "${params.run_id}.vcf.gz" | grep '#CHROM' | head -1 >> header.txt
    zcat "${params.run_id}.vcf.gz" | grep -v '#' | awk '
        BEGIN {FS = "\\t"} {
            \$7 = "PASS";
            \$8 = ".";
            gsub(":.*", "", \$9);
            gsub(":.*", "", \$10);
            gsub("\\\\|", "/", \$10);
            print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8 "\\t" \$9 "\\t" \$10
        }' > body.txt
    cat header.txt body.txt | bgzip > "${params.run_id}_fixed.vcf.gz"

    rm -f header.txt body.txt "${params.run_id}.vcf.gz"
    """
}
