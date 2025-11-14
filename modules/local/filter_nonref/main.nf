process FILTER_NONREF_VARIANTS {
    input:
    path vcf
    path tbi

    output:
    path "${params.run_id}.filt.nonref.vcf.gz", emit: vcf
    path "${params.run_id}.filt.nonref.vcf.gz.tbi", emit: tbi

    script:
    """
    # drop sites that are homozygous reference (GT=0/0) across all samples
    bcftools view -i 'COUNT(FMT/GT!="0/0")>0' $vcf -Oz -o ${params.run_id}.filt.nonref.vcf.gz

    variant_count=\$(bcftools view -H ${params.run_id}.filt.nonref.vcf.gz | wc -l)
    echo "Removed homozygous reference records. Variants remaining: \${variant_count}"
    tabix -p vcf ${params.run_id}.filt.nonref.vcf.gz
    """
}
