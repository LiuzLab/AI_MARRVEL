process SPLIT_VCF_BY_CHROMOSOME {
    input:
    path vcf

    output:
    path "chr*.vcf.gz", emit: chr_vcfs

    script:
    """
    # Get the list of chromosomes from the VCF file

    bcftools index ${vcf}

    bcftools query -f '%CHROM\n' ${vcf} | sort | uniq > chrom_list.txt
    # Split the VCF file by chromosome
    while read chrom; do
        bcftools view -r \${chrom} ${vcf} -Oz -o chr\${chrom}.vcf.gz
    done < chrom_list.txt
    """
}
