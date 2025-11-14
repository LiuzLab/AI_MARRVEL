include { SPLIT_VCF_BY_CHROMOSOME }    from "../../../modules/local/split_vcf_by_chromosome";
include { ANNOTATE_BY_VEP }            from "../../../modules/local/annotate_by_vep";
include { HPO_SIM }                    from "../../../modules/local/hpo_sim";
include { ANNOTATE_BY_MODULES }        from "../../../modules/local/annotate_by_modules";
include { ANNOTATE_TIER }              from "../../../modules/local/annotate_tier";
include { JOIN_PHRANK }                from "../../../modules/local/join_phrank";
include { MERGE_SCORES_BY_CHROMOSOME } from "../../../modules/local/merge_scores_by_chromosome";

include { PHRANK_SCORING }                from "../../../subworkflows/local/phrank_scoring";

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
