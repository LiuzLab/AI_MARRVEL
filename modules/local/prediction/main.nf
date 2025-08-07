process PREDICTION {
    publishDir "${params.outdir}/${params.run_id}/prediction/", mode: "copy"

    input:
    path merged_matrix
    path merged_compressed_scores
    path ref_model_inputs_dir

    output:
    path "conf_4Model/${params.run_id}_default_predictions.csv", emit: "default_predictions"
    path "conf_4Model/${params.run_id}_recessive_predictions.csv", optional: true
    path "conf_4Model/${params.run_id}_nd_predictions.csv", optional: true
    path "conf_4Model/${params.run_id}_nd_recessive_predictions.csv", optional: true
    path "final_matrix_expanded/*.csv.gz"
    path "shap_outputs"

    script:
    """
    mkdir final_matrix_expanded
    mkdir conf_4Model

    run_final.py \
        $ref_model_inputs_dir/default/final_model.job \
        $ref_model_inputs_dir/default/features.csv \
        $merged_matrix \
        ${params.run_id}.default_prediction.csv

    # Generate final matrix expanded
    merge_rm.py \
        ${params.run_id}.default_prediction.csv \
        $merged_compressed_scores \
        final_matrix_expanded/${params.run_id}.expanded.csv.gz

    # Generate conf_4Model
    extraModel_main.py -id ${params.run_id}
    """
}
