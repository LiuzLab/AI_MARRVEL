include {
    addDependentParams
} from "../../../modules/local/utils"

include {
    VALIDATE_VCF; NORMALIZE_VCF; GENERATE_INPUT_VCF;
    VCF_TO_VARIANTS; VARIANTS_TO_ENSEMBL; ENSEMBL_TO_GENESYM; GENESYM_TO_PHRANK;
    CONVERT_GVCF; FILTER_UNPASSED; FILTER_MITO_AND_UNKOWN_CHR; FILTER_BED; FILTER_PROBAND;
    SPLIT_VCF_BY_CHROMOSOME; ANNOTATE_BY_VEP; HPO_SIM; ANNOTATE_BY_MODULES;
    ANNOTATE_TIER; JOIN_PHRANK; MERGE_SCORES_BY_CHROMOSOME; PREDICTION
} from "../../../modules/local/singleton"


addDependentParams(params)

workflow VCF_PRE_PROCESS {
    take:
    input_vcf
    fasta
    fasta_index
    fasta_dict

    main:
    VALIDATE_VCF(input_vcf)
    NORMALIZE_VCF(VALIDATE_VCF.out.vcf)

    CONVERT_GVCF(
        NORMALIZE_VCF.out.vcf,
        NORMALIZE_VCF.out.tbi,
        fasta,
        fasta_index,
        fasta_dict,
        params.chrmap
    )
    FILTER_UNPASSED(
        CONVERT_GVCF.out.vcf,
        CONVERT_GVCF.out.tbi,
        params.chrmap
    )
    FILTER_MITO_AND_UNKOWN_CHR(
        FILTER_UNPASSED.out.vcf,
        FILTER_UNPASSED.out.tbi,
    )
    FILTER_BED(
        FILTER_MITO_AND_UNKOWN_CHR.out.vcf,
        FILTER_MITO_AND_UNKOWN_CHR.out.tbi,
        moduleDir.resolve(params.ref_filter_bed),
    )
    FILTER_PROBAND(
        FILTER_BED.out.vcf,
        FILTER_BED.out.tbi,
        params.ref_gnomad_genome,
        params.ref_gnomad_genome_idx,
        params.ref_gnomad_exome,
        params.ref_gnomad_exome_idx
    )

    emit:
    vcf = FILTER_PROBAND.out.vcf
}

workflow PHRANK_SCORING {
    take:
    vcf
    hpo

    main:
    VCF_TO_VARIANTS(vcf)
    VARIANTS_TO_ENSEMBL(VCF_TO_VARIANTS.out, params.ref_loc)
    ENSEMBL_TO_GENESYM(VARIANTS_TO_ENSEMBL.out, params.ref_to_sym)
    GENESYM_TO_PHRANK(ENSEMBL_TO_GENESYM.out,
                    hpo,
                    params.phrank_dagfile,
                    params.phrank_disease_annotation,
                    params.phrank_gene_annotation,
                    params.phrank_disease_gene)

    emit:
    phrank = GENESYM_TO_PHRANK.out
}

workflow GENERATE_SINGLETON_FEATURES {
    take:
    vcf
    hpo

    main:
    SPLIT_VCF_BY_CHROMOSOME(vcf)
    ANNOTATE_BY_VEP(
        SPLIT_VCF_BY_CHROMOSOME.out.chr_vcfs.flatten(),
        params.vep_dir_cache,
        params.vep_dir_plugins,
        params.vep_custom_gnomad,
        params.vep_custom_clinvar,
        params.vep_custom_hgmd,
        params.vep_plugin_revel,
        params.vep_plugin_spliceai_snv,
        params.vep_plugin_spliceai_indel,
        params.vep_plugin_cadd,
        params.vep_plugin_dbnsfp,
        file(params.vep_idx)
    )

    HPO_SIM(hpo,
            params.omim_hgmd_phen,
            params.omim_obo,
            params.omim_genemap2,
            params.omim_pheno)
    ANNOTATE_BY_MODULES (
        ANNOTATE_BY_VEP.out.vep_output,
        HPO_SIM.out.hgmd_sim,
        HPO_SIM.out.omim_sim,
        file(params.ref_annot_dir)
    )

    PHRANK_SCORING(vcf, hpo)
    ANNOTATE_TIER (
        ANNOTATE_BY_MODULES.out.scores,
        PHRANK_SCORING.out,
        file(params.ref_annot_dir),
        file(params.ref_var_tier_dir),
    )
    JOIN_PHRANK (
        ANNOTATE_BY_MODULES.out.scores,
        PHRANK_SCORING.out,
        file(params.ref_merge_expand_dir)
    )

    MERGE_SCORES_BY_CHROMOSOME(
        PHRANK_SCORING.out,
        ANNOTATE_TIER.out.tier.collect(),
        JOIN_PHRANK.out.compressed_scores.collect(),
        file(params.ref_annot_dir),
        file(params.ref_mod5_diffusion_dir),
        file(params.ref_merge_expand_dir)
    )

    emit:
    merged_matrix = MERGE_SCORES_BY_CHROMOSOME.out.merged_matrix
    merged_compressed_scores = MERGE_SCORES_BY_CHROMOSOME.out.merged_compressed_scores
}
