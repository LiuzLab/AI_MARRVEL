process ENSURE_VCF_NON_EMPTY {
    input:
    path vcf
    path tbi

    output:
    path vcf, emit: vcf

    script:
    """
    variant_count=\$(bcftools view -H $vcf | wc -l)
    if [ "\${variant_count}" -gt 0 ]; then
        echo "\${variant_count}" > ${params.run_id}.vcf_non_empty.txt
        echo "Final pre-processed VCF has \${variant_count} variants."
    else
        echo "All variants contained only GT=0/0 or were filtered. Nothing left to process." >&2
        echo "Provide a VCF with non-reference calls or relax upstream filtering thresholds." >&2
        exit 1
    fi
    """
}
