process EXPAND_PREDICTION_AND_COLLAPSE {
    publishDir "${params.outdir}/${params.run_id}/prediction/", mode: "copy"

    input:
    path default_prediction
    path merged_compressed_scores

    output:
    path "final_matrix_expanded/*.transcript.csv.gz", emit: "transcriptdf"
    path "final_matrix_expanded/*.variant.csv", emit: "vardf"
    path "final_matrix_expanded/*.gene.csv", emit: "genedf"

    script:
    """
    mkdir final_matrix_expanded

    merge_rankmatrix.py \
        $default_prediction \
        $merged_compressed_scores \
        final_matrix_expanded/${params.run_id}.expanded.transcript.csv.gz

    $moduleDir/bin/generate_varrank.py \
        --input final_matrix_expanded/${params.run_id}.expanded.transcript.csv.gz \
        --output final_matrix_expanded/${params.run_id}.variant.csv

    $moduleDir/bin/generate_generank.py \
        --input final_matrix_expanded/${params.run_id}.variant.csv \
        --output final_matrix_expanded/${params.run_id}.gene.csv

    """
}
