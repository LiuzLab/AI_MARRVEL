process FILTER_BED {
    input:
    path vcf
    path tbi
    path ref_filter_bed_file

    output:
    path "${params.run_id}.recode.vcf.gz", emit: vcf
    path "${params.run_id}.recode.vcf.gz.tbi", emit: tbi

    script:
    """
    if [[ -s ${ref_filter_bed_file} ]]; then  # Check if the bed filter file is not empty.
        awk '{gsub(/^chr/, ""); print}' ${ref_filter_bed_file} > bed
        bcftools filter --regions-file bed ${vcf} -Oz -o "${params.run_id}.recode.vcf.gz"
        tabix -p vcf "${params.run_id}.recode.vcf.gz"
    else
        cp ${vcf} "${params.run_id}.recode.vcf.gz"
        cp ${tbi} "${params.run_id}.recode.vcf.gz.tbi"
    fi
    """
}
