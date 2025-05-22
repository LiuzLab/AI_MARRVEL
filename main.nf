nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'

include {
    showVersion; addDependentParams
} from "./modules/local/utils"

include {
    BUILD_REFERENCE_INDEX; BUILD_REFERENCE_INDEX as BUILD_REFERENCE_INDEX_TRIO; PHENOPACKET_TO_VARIANTS_AND_HPOS; GENERATE_INPUT_VCF; GENERATE_INPUT_VARIANTS
} from "./modules/local/singleton"

include {
    VCF_PRE_PROCESS_TRIO; GENERATE_TRIO_FEATURES; PREDICTION_TRIO
} from "./modules/local/trio"

include {
    VCF_PRE_PROCESS; GENERATE_SINGLETON_FEATURES; PREDICTION
} from "./subworkflows/local/singleton"

showVersion()
validateParameters()
addDependentParams(params)

workflow {
    skip_preprocess_vcf_flag = false

    if (params.input_phenopacket) {
        skip_preprocess_vcf_flag = true

        PHENOPACKET_TO_VARIANTS_AND_HPOS(params.input_phenopacket)
        GENERATE_INPUT_VCF(PHENOPACKET_TO_VARIANTS_AND_HPOS.out.variants)
        vcf = GENERATE_INPUT_VCF.out.vcf
        hpo = PHENOPACKET_TO_VARIANTS_AND_HPOS.out.hpo
    } else if (params.input_variant) {
        skip_preprocess_vcf_flag = true

        GENERATE_INPUT_VARIANTS(params.input_variant)
        GENERATE_INPUT_VCF(GENERATE_INPUT_VARIANTS.out.variants)
        vcf = GENERATE_INPUT_VCF.out.vcf
        hpo = params.input_hpo
    } else if (params.input_ped && params.input_vcf) {
        BUILD_REFERENCE_INDEX_TRIO()
        VCF_PRE_PROCESS_TRIO(
            file(params.input_vcf),
            file(params.input_ped),
            BUILD_REFERENCE_INDEX_TRIO.out.fasta,
            BUILD_REFERENCE_INDEX_TRIO.out.fasta_index,
            BUILD_REFERENCE_INDEX_TRIO.out.fasta_dict,
            params.chrmap,
        )
        vcf = VCF_PRE_PROCESS_TRIO.out.vcf
        hpo = params.input_hpo
    } else if (params.input_vcf) {
        vcf = file(params.input_vcf)
        hpo = params.input_hpo
    } else {
        error "No input VCF or phenopacket or variant provided."
    }

    if (!hpo) {
        hpo = file(moduleDir.resolve("./assets/NO_FILE"))
    }

    if (!skip_preprocess_vcf_flag) {
        BUILD_REFERENCE_INDEX()
        VCF_PRE_PROCESS(
            vcf,
            BUILD_REFERENCE_INDEX.out.fasta,
            BUILD_REFERENCE_INDEX.out.fasta_index,
            BUILD_REFERENCE_INDEX.out.fasta_dict,
        )
        vcf = VCF_PRE_PROCESS.out.vcf
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
