process FILTER_PROBAND {
    publishDir "${params.outdir}/${params.run_id}/vcf/", mode: 'copy'

    input:
    path vcf
    path tbi
    tuple(
        path(ref_gnomad_genome),
        path(ref_gnomad_genome_idx),
        path(ref_gnomad_exome),
        path(ref_gnomad_exome_idx),
    )

    output:
    path "${params.run_id}.filt.rmBL.vcf.gz", emit: vcf

    script:
    """
    mkdir -m777 isec_tmp1
    bcftools isec -p isec_tmp1 -w 1 -Oz \\
    $vcf $ref_gnomad_genome
    # tabix isec_tmp1/0000.vcf.gz

    mkdir -m777 isec_tmp2
    bcftools isec -p isec_tmp2 -w 1 -Oz \\
    $vcf $ref_gnomad_exome
    # tabix isec_tmp2/0000.vcf.gz

    mkdir -m777 isec_tmp3
    bcftools isec -p isec_tmp3 -Ov \\
    isec_tmp1/0000.vcf.gz isec_tmp2/0000.vcf.gz

    mv isec_tmp3/0002.vcf ${params.run_id}.filt.rmBL.vcf
    bgzip ${params.run_id}.filt.rmBL.vcf
    """

}
