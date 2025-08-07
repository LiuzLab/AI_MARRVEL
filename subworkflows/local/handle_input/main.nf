include { PHENOPACKET_TO_VARIANTS_AND_HPOS } from "../../../modules/local/phenopacket_to_variants_and_hpo";
include { GENERATE_INPUT_VCF }               from "../../../modules/local/generate_input_vcf";
include { GENERATE_INPUT_VARIANTS }          from "../../../modules/local/generate_input_variants";

workflow HANDLE_INPUT {
    if (params.input_phenopacket) {
        PHENOPACKET_TO_VARIANTS_AND_HPOS(params.input_phenopacket)
        GENERATE_INPUT_VCF(PHENOPACKET_TO_VARIANTS_AND_HPOS.out.variants)
        vcf = GENERATE_INPUT_VCF.out.vcf
        hpo = PHENOPACKET_TO_VARIANTS_AND_HPOS.out.hpo
    } else if (params.input_variant) {
        GENERATE_INPUT_VARIANTS(params.input_variant)
        GENERATE_INPUT_VCF(GENERATE_INPUT_VARIANTS.out.variants)
        vcf = GENERATE_INPUT_VCF.out.vcf
        hpo = params.input_hpo
    } else if (params.input_vcf) {
        vcf = file(params.input_vcf)
        hpo = params.input_hpo
    } else {
        error "No input VCF or phenopacket or variant provided."
    }

    if (!hpo) {
        hpo = file(moduleDir.resolve("../../../assets/NO_FILE"))
    }

    emit:
    vcf = vcf
    hpo = hpo
}
