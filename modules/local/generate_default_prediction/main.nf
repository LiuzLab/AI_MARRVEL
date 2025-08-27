process GENERATE_DEFAULT_PREDICTION {
    input:
    path merged_matrix
    path ref_model_inputs_dir

    output:
    path "${params.run_id}.default_prediction.csv", emit: "default_prediction"

    script:
    """
    # Generate default prediction output
    run_final.py \
        $ref_model_inputs_dir/default/final_model.job \
        $ref_model_inputs_dir/default/features.csv \
        $merged_matrix \
        ${params.run_id}.default_prediction.csv
    """
}