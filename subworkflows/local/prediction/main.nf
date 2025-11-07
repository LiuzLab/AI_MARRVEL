include { GENERATE_DEFAULT_PREDICTION } from '../../../modules/local/generate_default_prediction'
include { EXPAND_PREDICTION_AND_COLLAPSE } from '../../../modules/local/expand_prediction_and_collapse'
include { GENERATE_GENERANK } from '../../../modules/local/generate_generank'
include { GENERATE_EXTRAMODEL_PREDICTIONS_AND_SHAP_INTERPRETATION } from '../../../modules/local/generate_extramodel_predictions_and_shap_interpretation'


workflow PREDICTION {
    take:
    merged_matrix
    merged_compressed_scores
    ref_model_inputs_dir

    main:
    GENERATE_DEFAULT_PREDICTION(merged_matrix, ref_model_inputs_dir)
    EXPAND_PREDICTION_AND_COLLAPSE(
        GENERATE_DEFAULT_PREDICTION.out.default_prediction,
        merged_compressed_scores,
    )

    GENERATE_EXTRAMODEL_PREDICTIONS_AND_SHAP_INTERPRETATION(
        GENERATE_DEFAULT_PREDICTION.out.default_prediction,
        EXPAND_PREDICTION_AND_COLLAPSE.out.transcriptdf,
        ref_model_inputs_dir,
    )

    emit:
    default_predictions = GENERATE_DEFAULT_PREDICTION.out.default_prediction
}