include { VCF_TO_VARIANTS }            from "../../../modules/local/vcf_to_variants";
include { VARIANTS_TO_ENSEMBL }        from "../../../modules/local/variants_to_ensembl";
include { ENSEMBL_TO_GENESYM }         from "../../../modules/local/ensembl_to_genesym";
include { GENESYM_TO_PHRANK }          from "../../../modules/local/genesym_to_phrank";

workflow PHRANK_SCORING {
    take:
    vcf
    hpo
    data

    main:
    ensembl_to_location_file = data.map { it.ensembl_to_location_file }
    ensembl_to_symbol_file = data.map { it.ensembl_to_symbol_file }
    phrank_tuple = data.map { it.phrank_tuple }
    VCF_TO_VARIANTS(vcf)
    VARIANTS_TO_ENSEMBL(
        VCF_TO_VARIANTS.out,
        ensembl_to_location_file,
    )
    ENSEMBL_TO_GENESYM(
        VARIANTS_TO_ENSEMBL.out,
        ensembl_to_symbol_file,
    )
    GENESYM_TO_PHRANK(
        ENSEMBL_TO_GENESYM.out,
        hpo,
        phrank_tuple,
    )

    emit:
    phrank = GENESYM_TO_PHRANK.out
}
