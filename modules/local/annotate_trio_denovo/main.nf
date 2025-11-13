process ANNOTATE_TRIO_DENOVO {
    container "broadinstitute/gatk"

    input:
    path joint_vcf
    path joint_tbi
    path ped
    path fasta
    path fasta_index
    path fasta_dict

    output:
    path "${params.run_id}.DeNovo.g.vcf", emit: denovo_vcf

    script:
    """
    gatk VariantAnnotator \
        -A PossibleDeNovo \
        -R "${fasta}" \
        --pedigree "$ped" \
        -V "${joint_vcf}" \
        -O "${params.run_id}.DeNovo.g.vcf"
    """
}
