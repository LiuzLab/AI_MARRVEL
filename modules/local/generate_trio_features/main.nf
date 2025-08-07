process GENERATE_TRIO_FEATURES {
    publishDir "${params.outdir}/${params.run_id}/features_trio", mode: 'copy'

    input:
    path merged_compressed_scores
    path merged_matrix
    path inheritance

    output:
    path "${params.run_id}.triomatrix.tsv", emit: triomatrix

    script:
    """
    trio_merge_inheritance.py \
        $merged_compressed_scores \
        $merged_matrix \
        $inheritance \
        ./${params.run_id}.triomatrix.tsv
    """
}
