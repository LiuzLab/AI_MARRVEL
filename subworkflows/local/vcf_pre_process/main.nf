include { NORMALIZE_VCF }              from "../../../modules/local/normalize_vcf";
include { CONVERT_GVCF }               from "../../../modules/local/convert_gvcf";
include { FILTER_UNPASSED }            from "../../../modules/local/filter_unpassed";
include { FILTER_MITO_AND_UNKOWN_CHR } from "../../../modules/local/filter_mito_and_unkown_chr";
include { FILTER_BED }                 from "../../../modules/local/filter_bed";
include { FILTER_PROBAND }             from "../../../modules/local/filter_proband";

workflow VCF_PRE_PROCESS {
    take:
    input_vcf
    data

    main:
    chrmap_file = data.map { it.chrmap_file }
    ref_filter_bed_file = data.map { it.ref_filter_bed_file }
    fasta_tuple = data.map { it.fasta_tuple }
    gnomad_tuple = data.map { it.gnomad_tuple }

    NORMALIZE_VCF(input_vcf)

    CONVERT_GVCF(
        NORMALIZE_VCF.out.vcf,
        NORMALIZE_VCF.out.tbi,
        fasta_tuple,
        chrmap_file
    )
    FILTER_UNPASSED(
        CONVERT_GVCF.out.vcf,
        CONVERT_GVCF.out.tbi,
        chrmap_file
    )
    FILTER_MITO_AND_UNKOWN_CHR(
        FILTER_UNPASSED.out.vcf,
        FILTER_UNPASSED.out.tbi,
    )
    FILTER_BED(
        FILTER_MITO_AND_UNKOWN_CHR.out.vcf,
        FILTER_MITO_AND_UNKOWN_CHR.out.tbi,
        ref_filter_bed_file,
    )
    FILTER_PROBAND(
        FILTER_BED.out.vcf,
        FILTER_BED.out.tbi,
        gnomad_tuple,
    )

    emit:
    vcf = FILTER_PROBAND.out.vcf
}
