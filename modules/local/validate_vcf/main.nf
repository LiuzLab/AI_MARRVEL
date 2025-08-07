process VALIDATE_VCF {
    container 'quay.io/biocontainers/vcftools:0.1.16--pl5321hdcf5f25_11'

    input:
    path vcf

    output:
    path "$vcf", emit: vcf

    """
    vcf-validator $vcf
    """
}
