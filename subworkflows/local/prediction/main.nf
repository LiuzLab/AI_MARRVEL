include { GENERATE_DEFAULT_PREDICTION } from '../../../modules/local/generate_default_prediction'
include { MERGE_PREDICTION_TO_TRANSCRIPT } from '../../../modules/local/merge_prediction_to_transcript'
include { GENERATE_GENERANK } from '../../../modules/local/generate_generank'
include { GENERATE_EXTRAMODEL_PREDICTIONS_AND_SHAP_INTERPRETATION } from '../../../modules/local/generate_extramodel_predictions_and_shap_interpretation'


workflow PREDICTION {
    take:
    merged_matrix
    merged_compressed_scores
    ref_model_inputs_dir

    main:
    GENERATE_DEFAULT_PREDICTION(merged_matrix, ref_model_inputs_dir)
    MERGE_PREDICTION_TO_TRANSCRIPT(
        GENERATE_DEFAULT_PREDICTION.out.default_prediction,
        merged_compressed_scores,
    )

    GENERATE_GENERANK(
        MERGE_PREDICTION_TO_TRANSCRIPT.out.final_matrix_expanded,
    )

    GENERATE_EXTRAMODEL_PREDICTIONS_AND_SHAP_INTERPRETATION(
        GENERATE_DEFAULT_PREDICTION.out.default_prediction,
        MERGE_PREDICTION_TO_TRANSCRIPT.out.final_matrix_expanded,
        ref_model_inputs_dir,
    )

    emit:
    default_predictions = GENERATE_DEFAULT_PREDICTION.out.default_prediction
}