process GENERATE_EXTRAMODEL_PREDICTIONS_AND_SHAP_INTERPRETATION {
    label "big_mem"

    publishDir "${params.outdir}/${params.run_id}/prediction/", mode: "copy"

    input:
    path default_prediction
    path final_matrix_expanded
    path ref_model_inputs_dir

    output:
    path "conf_4Model/${params.run_id}_default_predictions.csv"
    path "conf_4Model/${params.run_id}_recessive_predictions.csv", optional: true
    path "conf_4Model/${params.run_id}_nd_predictions.csv", optional: true
    path "conf_4Model/${params.run_id}_nd_recessive_predictions.csv", optional: true
    path "shap_outputs"

    script:
    """
    # Generate 4 outputs under conf_4Model directory, and perform SHAP interpretation
    extraModel_main.py -id ${params.run_id}
    """
}