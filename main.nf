nextflow.enable.dsl = 2

include { validateParameters }         from 'plugin/nf-schema'

include { showVersion}                 from "./modules/local/utils"
include { VCF_PRE_PROCESS_TRIO }       from "./modules/local/vcf_pre_process_trio";
include { GENERATE_TRIO_FEATURES }     from "./modules/local/generate_trio_features";
include { PREDICTION_TRIO }            from "./modules/local/prediction_trio";
include { PREDICTION }                 from "./modules/local/prediction";

include { PREPARE_DATA}                from "./subworkflows/local/prepare_data"
include { HANDLE_INPUT}                from "./subworkflows/local/handle_input"
include { VCF_PRE_PROCESS}             from "./subworkflows/local/vcf_pre_process"
include { GENERATE_SINGLETON_FEATURES} from "./subworkflows/local/generate_singleton_feature"

showVersion()
validateParameters()

workflow {
    data = PREPARE_DATA()
    chrmap_file = data.map { it.chrmap_file }
    ref_model_inputs_dir = data.map { it.ref_model_inputs_dir }
    fasta_tuple = data.map { it.fasta_tuple }

    (vcf, hpo) = HANDLE_INPUT()

    if (params.input_ped) {
        (vcf, inheritance) = VCF_PRE_PROCESS_TRIO(
            vcf,
            file(params.input_ped),
            fasta_tuple.map { it[0] },
            fasta_tuple.map { it[1] },
            fasta_tuple.map { it[2] },
            chrmap_file,
        )
    }

    if (!params.input_variant && !params.input_phenopacket) {
        vcf = VCF_PRE_PROCESS(
            vcf,
            data,
        )
    }

    GENERATE_SINGLETON_FEATURES(vcf, hpo, data)
    PREDICTION(
        GENERATE_SINGLETON_FEATURES.out.merged_matrix,
        GENERATE_SINGLETON_FEATURES.out.merged_compressed_scores,
        ref_model_inputs_dir,
    )

    if (params.input_ped) {
        GENERATE_TRIO_FEATURES(
            GENERATE_SINGLETON_FEATURES.out.merged_compressed_scores,
            PREDICTION.out.default_predictions,
            inheritance,
        )
        PREDICTION_TRIO(
            GENERATE_SINGLETON_FEATURES.out.merged_compressed_scores,
            GENERATE_TRIO_FEATURES.out.triomatrix,
        )
    }
}
