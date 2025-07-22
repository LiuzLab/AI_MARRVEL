include {
    VALIDATE_VCF; NORMALIZE_VCF; GENERATE_INPUT_VCF;
    VCF_TO_VARIANTS; VARIANTS_TO_ENSEMBL; ENSEMBL_TO_GENESYM; GENESYM_TO_PHRANK;
    CONVERT_GVCF; FILTER_UNPASSED; FILTER_MITO_AND_UNKOWN_CHR; FILTER_BED; FILTER_PROBAND;
    SPLIT_VCF_BY_CHROMOSOME; ANNOTATE_BY_VEP; HPO_SIM; ANNOTATE_BY_MODULES;
    ANNOTATE_TIER; JOIN_PHRANK; MERGE_SCORES_BY_CHROMOSOME; PREDICTION
} from "../../../modules/local/singleton"

workflow VCF_PRE_PROCESS {
    take:
    input_vcf
    data

    main:
    chrmap_file = data.map { it.chrmap_file }
    ref_filter_bed_file = data.map { it.ref_filter_bed_file }
    fasta_tuple = data.map { it.fasta_tuple }
    gnomad_tuple = data.map { it.gnomad_tuple }

    VALIDATE_VCF(input_vcf)
    NORMALIZE_VCF(VALIDATE_VCF.out.vcf)

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

workflow GENERATE_SINGLETON_FEATURES {
    take:
    vcf
    hpo
    data

    main:
    ref_annot_dir = data.map { it.ref_annot_dir }
    ref_var_tier_dir = data.map { it.ref_var_tier_dir }
    ref_merge_expand_dir = data.map { it.ref_merge_expand_dir }
    ref_mod5_diffusion_dir = data.map { it.ref_mod5_diffusion_dir }
    omim_tuple = data.map { it.omim_tuple }
    vep_tuple = data.map { it.vep_tuple }

    SPLIT_VCF_BY_CHROMOSOME(vcf)
    ANNOTATE_BY_VEP(
        SPLIT_VCF_BY_CHROMOSOME.out.chr_vcfs.flatten(),
        vep_tuple,
    )

    HPO_SIM(
        hpo,
        omim_tuple,
    )
    ANNOTATE_BY_MODULES (
        ANNOTATE_BY_VEP.out.vep_output,
        HPO_SIM.out.hgmd_sim,
        HPO_SIM.out.omim_sim,
        ref_annot_dir,
    )

    PHRANK_SCORING(vcf, hpo, data)
    ANNOTATE_TIER (
        ANNOTATE_BY_MODULES.out.scores,
        PHRANK_SCORING.out,
        ref_annot_dir,
        ref_var_tier_dir,
    )
    JOIN_PHRANK (
        ANNOTATE_BY_MODULES.out.scores,
        PHRANK_SCORING.out,
        ref_merge_expand_dir,
    )

    MERGE_SCORES_BY_CHROMOSOME(
        PHRANK_SCORING.out,
        ANNOTATE_TIER.out.tier.collect(),
        JOIN_PHRANK.out.compressed_scores.collect(),
        ref_annot_dir,
        ref_mod5_diffusion_dir,
        ref_merge_expand_dir,
    )

    emit:
    merged_matrix = MERGE_SCORES_BY_CHROMOSOME.out.merged_matrix
    merged_compressed_scores = MERGE_SCORES_BY_CHROMOSOME.out.merged_compressed_scores
}
