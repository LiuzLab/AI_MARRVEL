nextflow.enable.dsl = 2

include {
    showUsage; showVersion; validateInputParams; addDependentParams
} from "./modules/local/utils"

include {
    BUILD_REFERENCE_INDEX
} from "./modules/local/singleton"

include {
    VCF_PRE_PROCESS_TRIO; GENERATE_TRIO_FEATURES; PREDICTION_TRIO
} from "./modules/local/trio"

include {
    VCF_PRE_PROCESS; GENERATE_SINGLETON_FEATURES; PREDICTION
} from "./subworkflows/local/singleton"

showUsage()
showVersion()
validateInputParams()
addDependentParams(params)

workflow {
    hpo = params.input_hpo
    if (params.input_ped) {
        BUILD_REFERENCE_INDEX()
        VCF_PRE_PROCESS_TRIO(
            file(params.input_vcf),
            file(params.input_ped),
            BUILD_REFERENCE_INDEX.out.fasta,
            BUILD_REFERENCE_INDEX.out.fasta_index,
            BUILD_REFERENCE_INDEX.out.fasta_dict,
            params.chrmap,
        )
        VCF_PRE_PROCESS(
            VCF_PRE_PROCESS_TRIO.out.vcf,
            BUILD_REFERENCE_INDEX.out.fasta,
            BUILD_REFERENCE_INDEX.out.fasta_index,
            BUILD_REFERENCE_INDEX.out.fasta_dict,
        )
        vcf = VCF_PRE_PROCESS.out.vcf
    } else if (params.input_vcf) {
        BUILD_REFERENCE_INDEX()
        VCF_PRE_PROCESS(
            file(params.input_vcf),
            BUILD_REFERENCE_INDEX.out.fasta,
            BUILD_REFERENCE_INDEX.out.fasta_index,
            BUILD_REFERENCE_INDEX.out.fasta_dict,
        )
        vcf = VCF_PRE_PROCESS.out.vcf
    } else if (params.input_variant) {
        GENERATE_INPUT_VCF(params.input_variant)
        vcf = GENERATE_INPUT_VCF.out.vcf
    }

    GENERATE_SINGLETON_FEATURES(vcf, hpo)
    PREDICTION(
        GENERATE_SINGLETON_FEATURES.out.merged_matrix,
        GENERATE_SINGLETON_FEATURES.out.merged_compressed_scores,
        file(params.ref_model_inputs_dir)
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
