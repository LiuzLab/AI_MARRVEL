process FILTER_UNPASSED {
    input:
    path vcf
    path tbi
    path chrmap

    output:
    path "${params.run_id}.filt.vcf.gz", emit: vcf
    path "${params.run_id}.filt.vcf.gz.tbi", emit: tbi

    script:
    """
    # Annotate with new chromosomes and preserve original coordinates in ID
    bcftools annotate --rename-chrs $chrmap -x ID $vcf -Oz -o ${params.run_id}-annot

    # Annotate with new IDs based on CHROM, POS, REF, ALT
    bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' ${params.run_id}-annot -Oz -o ${params.run_id}-add-id.vcf.gz

    #run quality filters
    bcftools filter  ${params.run_id}-add-id.vcf.gz -i'FILTER == "PASS"' -Oz -o ${params.run_id}.filt.vcf.gz

    #check number of variants left after quality filtering
    variant_count=\$(bcftools view -H ${params.run_id}.filt.vcf.gz | wc -l)
    if [ "\${variant_count}" -gt 0 ]; then
        echo "Quality filtering completed successfully. Variants passing the quality filters: \${variant_count}"
    else
        echo "The VCF file doesn't have FILTER annotation, or all variants filtered."
        echo "Pipeline will proceed with unfiltered VCF file."
        cp ${params.run_id}-add-id.vcf.gz ${params.run_id}.filt.vcf.gz
    fi

    # drop sites that are homozygous reference (GT=0/0) across all samples
    bcftools view -i 'COUNT(FMT/GT!="0/0")>0' ${params.run_id}.filt.vcf.gz -Oz -o ${params.run_id}.filt.nonref.vcf.gz
    mv ${params.run_id}.filt.nonref.vcf.gz ${params.run_id}.filt.vcf.gz

    variant_count=\$(bcftools view -H ${params.run_id}.filt.vcf.gz | wc -l)
    if [ "\${variant_count}" -gt 0 ]; then
        echo "Removed homozygous reference records. Variants remaining: \${variant_count}"
    else
        echo "All variants contained only GT=0/0 across chromosomes. Nothing left to process." >&2
        echo "Provide a VCF with non-reference calls or relax upstream filtering thresholds." >&2
        exit 1
    fi
    tabix -p vcf ${params.run_id}.filt.vcf.gz
    """
}
