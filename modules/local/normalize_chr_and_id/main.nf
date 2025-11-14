process NORMALIZE_CHR_AND_ID {
    input:
    path vcf
    path tbi
    path chrmap

    output:
    path "${params.run_id}.chrid.vcf.gz", emit: vcf
    path "${params.run_id}.chrid.vcf.gz.tbi", emit: tbi

    script:
    """
    # Annotate with new chromosomes and drop the existing ID
    bcftools annotate --rename-chrs $chrmap -x ID $vcf -Oz -o ${params.run_id}.chrid.tmp.vcf.gz

    # Annotate with new IDs based on CHROM, POS, REF, ALT
    bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' ${params.run_id}.chrid.tmp.vcf.gz -Oz -o ${params.run_id}.chrid.vcf.gz

    # Index normalized VCF and clean up temp file
    tabix -p vcf ${params.run_id}.chrid.vcf.gz
    rm -f ${params.run_id}.chrid.tmp.vcf.gz
    """
}
