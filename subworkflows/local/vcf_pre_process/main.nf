include { ENSURE_BGZIP_VCF }           from "../../../modules/local/normalize_vcf";
include { CONVERT_GVCF }               from "../../../modules/local/convert_gvcf";
include { NORMALIZE_CHR_AND_ID }       from "../../../modules/local/normalize_chr_and_id";
include { FILTER_UNPASSED }            from "../../../modules/local/filter_unpassed";
include { FILTER_NONREF_VARIANTS }     from "../../../modules/local/filter_nonref";
include { FILTER_MITO_AND_UNKNOWN_CHR } from "../../../modules/local/filter_mito_and_unknown_chr";
include { FILTER_BED }                 from "../../../modules/local/filter_bed";
include { FILTER_BLACKLIST_VARIANTS }  from "../../../modules/local/filter_blacklist_variants";
include { ENSURE_VCF_NON_EMPTY }       from "../../../modules/local/ensure_vcf_non_empty";

workflow VCF_PRE_PROCESS {
    take:
    input_vcf
    data

    main:
    chrmap_file = data.map { it.chrmap_file }
    ref_filter_bed_file = data.map { it.ref_filter_bed_file }
    fasta_tuple = data.map { it.fasta_tuple }
    gnomad_tuple = data.map { it.gnomad_tuple }

    ENSURE_BGZIP_VCF(input_vcf)

    CONVERT_GVCF(
        ENSURE_BGZIP_VCF.out.vcf,
        ENSURE_BGZIP_VCF.out.tbi,
        fasta_tuple,
        chrmap_file
    )
    NORMALIZE_CHR_AND_ID(
        CONVERT_GVCF.out.vcf,
        CONVERT_GVCF.out.tbi,
        chrmap_file
    )
    FILTER_UNPASSED(
        NORMALIZE_CHR_AND_ID.out.vcf,
        NORMALIZE_CHR_AND_ID.out.tbi,
    )
    FILTER_NONREF_VARIANTS(
        FILTER_UNPASSED.out.vcf,
        FILTER_UNPASSED.out.tbi,
    )
    FILTER_MITO_AND_UNKNOWN_CHR(
        FILTER_NONREF_VARIANTS.out.vcf,
        FILTER_NONREF_VARIANTS.out.tbi,
    )
    FILTER_BED(
        FILTER_MITO_AND_UNKNOWN_CHR.out.vcf,
        FILTER_MITO_AND_UNKNOWN_CHR.out.tbi,
        ref_filter_bed_file,
    )
    FILTER_BLACKLIST_VARIANTS(
        FILTER_BED.out.vcf,
        FILTER_BED.out.tbi,
        gnomad_tuple,
    )
    ENSURE_VCF_NON_EMPTY(
        FILTER_BLACKLIST_VARIANTS.out.vcf,
        FILTER_BLACKLIST_VARIANTS.out.tbi,
    )

    emit:
    vcf = ENSURE_VCF_NON_EMPTY.out.vcf
}
