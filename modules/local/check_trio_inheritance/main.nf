process CHECK_TRIO_INHERITANCE {
    container "broadinstitute/gatk"

    input:
    path denovo_vcf
    path sample_ids

    output:
    path "${params.run_id}.inheritance.txt", emit: inheritance

    script:
    """
    set -euo pipefail

    source "${sample_ids}"

    trio_check_inheritance.GATK.py \
        "${params.run_id}" \
        "${denovo_vcf}" \
        "${params.ref_ver}" \
        "\${SAMPLE_P}" \
        "\${SAMPLE_F}" \
        "\${SAMPLE_M}"
    """
}
