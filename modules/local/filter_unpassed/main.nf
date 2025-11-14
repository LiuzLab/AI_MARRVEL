process FILTER_UNPASSED {
    input:
    path vcf
    path tbi

    output:
    path "${params.run_id}.filt.vcf.gz", emit: vcf
    path "${params.run_id}.filt.vcf.gz.tbi", emit: tbi

    script:
    """
    #run quality filters
    bcftools filter  $vcf -i'FILTER == "PASS"' -Oz -o ${params.run_id}.filt.vcf.gz

    #check number of variants left after quality filtering
    variant_count=\$(bcftools view -H ${params.run_id}.filt.vcf.gz | wc -l)
    if [ "\${variant_count}" -gt 0 ]; then
        echo "Quality filtering completed successfully. Variants passing the quality filters: \${variant_count}"
    else
        echo "The VCF file doesn't have FILTER annotation, or all variants filtered."
        echo "Pipeline will proceed with unfiltered VCF file."
        cp $vcf ${params.run_id}.filt.vcf.gz
    fi
    tabix -p vcf ${params.run_id}.filt.vcf.gz
    """
}
