process MERGE_PREDICTION_TO_TRANSCRIPT {
    publishDir "${params.outdir}/${params.run_id}/prediction/", mode: "copy"

    input:
    path default_prediction
    path merged_compressed_scores

    output:
    path "final_matrix_expanded/*.csv.gz", emit: "final_matrix_expanded"

    script:
    """
    mkdir final_matrix_expanded

    # Generate final matrix expanded
    merge_rm.py \
        $default_prediction \
        $merged_compressed_scores \
        final_matrix_expanded/${params.run_id}.expanded.csv.gz
    """
}
