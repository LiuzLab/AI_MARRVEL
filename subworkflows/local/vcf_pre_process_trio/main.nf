include { PARSE_TRIO_PED }         from "../../../modules/local/parse_trio_ped";
include { PREPARE_TRIO_VCF }       from "../../../modules/local/prepare_trio_vcf";
include { ANNOTATE_TRIO_DENOVO }   from "../../../modules/local/annotate_trio_denovo";
include { EXTRACT_PROBAND_VCF }    from "../../../modules/local/extract_proband_vcf";
include { CHECK_TRIO_INHERITANCE } from "../../../modules/local/check_trio_inheritance";

workflow VCF_PRE_PROCESS_TRIO {
    take:
    vcf
    ped
    fasta
    fasta_index
    fasta_dict
    chrmap_file

    main:
    PARSE_TRIO_PED(ped)
    sample_ids = PARSE_TRIO_PED.out.sample_ids

    PREPARE_TRIO_VCF(
        vcf,
        chrmap_file,
        sample_ids,
    )

    ANNOTATE_TRIO_DENOVO(
        PREPARE_TRIO_VCF.out.vcf,
        PREPARE_TRIO_VCF.out.tbi,
        ped,
        fasta,
        fasta_index,
        fasta_dict,
    )

    EXTRACT_PROBAND_VCF(
        PREPARE_TRIO_VCF.out.vcf,
        sample_ids,
    )

    CHECK_TRIO_INHERITANCE(
        ANNOTATE_TRIO_DENOVO.out.denovo_vcf,
        sample_ids,
    )

    emit:
    vcf = EXTRACT_PROBAND_VCF.out.vcf
    inheritance = CHECK_TRIO_INHERITANCE.out.inheritance
}
