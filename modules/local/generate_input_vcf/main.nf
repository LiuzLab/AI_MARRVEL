process GENERATE_INPUT_VCF {
    input:
    val variants

    output:
    path "input.vcf.gz", emit: vcf

    """
    generate_input_vcf.py $variants
    bgzip input.vcf
    """
}
