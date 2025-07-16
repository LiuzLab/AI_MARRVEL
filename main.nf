nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'

include {
    showVersion
} from "./modules/local/utils"

include {
    VCF_PRE_PROCESS_TRIO; GENERATE_TRIO_FEATURES; PREDICTION_TRIO
} from "./modules/local/trio"

include {
    PREPARE_DATA
} from "./subworkflows/local/prepare_data"

include {
    HANDLE_INPUT
} from "./subworkflows/local/handle_input"

include {
    VCF_PRE_PROCESS; GENERATE_SINGLETON_FEATURES; PREDICTION
} from "./subworkflows/local/singleton"

showVersion()
validateParameters()

workflow {
    data = PREPARE_DATA()
    chrmap_file = data.map { it.chrmap_file }
    ref_model_inputs_dir = data.map { it.ref_model_inputs_dir }
    fasta_tuple = data.map { it.fasta_tuple }

    (vcf, hpo) = HANDLE_INPUT()

    if (params.input_ped) {
        VCF_PRE_PROCESS_TRIO(
            vcf,
            file(params.input_ped),
            fasta_tuple.map { it[0] },
            fasta_tuple.map { it[1] },
            fasta_tuple.map { it[2] },
            chrmap_file,
        )
        vcf = VCF_PRE_PROCESS_TRIO.out.vcf
    }

    if (!params.input_variant && !params.input_phenopacket) {
        VCF_PRE_PROCESS(
            vcf,
            data,
        )
        vcf = VCF_PRE_PROCESS.out.vcf
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
            VCF_PRE_PROCESS_TRIO.out.inheritance,
        )
        PREDICTION_TRIO(
            GENERATE_SINGLETON_FEATURES.out.merged_compressed_scores,
            GENERATE_TRIO_FEATURES.out.triomatrix,
        )
    }
}
