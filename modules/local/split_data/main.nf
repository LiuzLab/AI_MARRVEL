process SPLIT_DATA {
    input:
    path "data_except_vep/*"
    path "data_only_vep/*"
    path bed_filter_file // NOTE(JL): Not tested yet
    val exome_filter // NOTE(JL): Not tested yet

    output:
    path "data_except_vep/bcf_annotate/chrmap.txt", emit: chrmap_file
    path "ensembl_to_location.txt", emit: ensembl_to_location_file
    path "data_except_vep/phrank/${params.ref_ver}/ensembl_to_symbol.txt", emit: ensembl_to_symbol_file
    path "ref_filter.bed", emit: ref_filter_bed_file

    path "data_except_vep/annotate", emit: ref_annot_dir
    path "data_except_vep/var_tier", emit: ref_var_tier_dir
    path "data_except_vep/merge_expand", emit: ref_merge_expand_dir
    path "data_except_vep/mod5_diffusion", emit: ref_mod5_diffusion_dir
    path "data_except_vep/predict_new", emit: ref_predict_new_dir
    path "data_except_vep/model_inputs", emit: ref_model_inputs_dir

    // for phrank
    tuple(
        path("data_except_vep/phrank/${params.ref_ver}/child_to_parent.txt"), // dagfile
        path("data_except_vep/phrank/${params.ref_ver}/disease_to_pheno.txt"), // disease_annotation
        path("data_except_vep/phrank/${params.ref_ver}/gene_to_phenotype.txt"), // gene_annotation
        path("data_except_vep/phrank/${params.ref_ver}/disease_to_gene.txt"), // disease_gene
        emit: phrank_tuple
    )

    // OMIM and HPO
    tuple(
        path("data_except_vep/omim_annotate/${params.ref_ver}/HGMD_phen.tsv"), // omim_hgmd_phen
        path("data_except_vep/omim_annotate/hp.obo"), // omim_obo
        path("data_except_vep/omim_annotate/${params.ref_ver}/genemap2_v2022.rds"), // omim_genemap2
        path("data_except_vep/omim_annotate/${params.ref_ver}/HPO_OMIM.tsv"), // omim_pheno
        emit: omim_tuple
    )

    // GNOMAD VCF
    tuple(
        path("data_except_vep/filter_vep/${params.ref_ver}/gnomad.${params.ref_ver}.blacklist.genomes.vcf.gz"), // ref_gnomad_genome
        path("data_except_vep/filter_vep/${params.ref_ver}/gnomad.${params.ref_ver}.blacklist.genomes.vcf.gz.tbi"), // ref_gnomad_genome_idx
        path("data_except_vep/filter_vep/${params.ref_ver}/gnomad.${params.ref_ver}.blacklist.exomes.vcf.gz"), // ref_gnomad_exome
        path("data_except_vep/filter_vep/${params.ref_ver}/gnomad.${params.ref_ver}.blacklist.exomes.vcf.gz.tbi"), // ref_gnomad_exome_idx
        emit: gnomad_tuple
    )

    // VEP
    tuple(
        path("data_only_vep/vep/${params.ref_ver}/"), // vep_dir_cache
        path("data_only_vep/vep/${params.ref_ver}/Plugins/"), // vep_dir_plugins
        path("gnomad.genomes*.vcf.gz", arity: 1), // vep_custom_gnomad
        path("data_only_vep/vep/${params.ref_ver}/clinvar_20220730.vcf.gz"), // vep_custom_clinvar
        path("data_only_vep/vep/${params.ref_ver}/HGMD_Pro_2022.2_${params.ref_ver}.vcf.gz"), // vep_custom_hgmd
        path("new_tabbed_revel*.tsv.gz", arity: 1), // vep_plugin_revel
        path("data_only_vep/vep/${params.ref_ver}/spliceai_scores.masked.snv.${params.ref_ver}.vcf.gz"), // vep_plugin_spliceai_snv
        path("data_only_vep/vep/${params.ref_ver}/spliceai_scores.masked.indel.${params.ref_ver}.vcf.gz"), // vep_plugin_spliceai_indel
        path("*_whole_genome_*.tsv.gz", arity: 1), // vep_plugin_cadd
        path("dbNSFP*.gz", arity: 1), // vep_plugin_dbnsfp
        path("data_only_vep/vep/${params.ref_ver}/*.tbi"), // vep_idx
        emit: vep_tuple
    )

    script:
    """
    # tree .

    ref_assembly=\$( [ "${params.ref_ver}" = "hg19" ] && echo "grch37" || echo "grch38" )
    ln -s data_except_vep/phrank/${params.ref_ver}/\${ref_assembly}_symbol_to_location.txt ./ensembl_to_location.txt

    if [[ -s "${bed_filter_file}" ]]; then # Check if the bed filter file is not empty. 
        ln -s "${bed_filter_file}" ref_filter.bed
    elif [[ "${exome_filter}" == "true" ]]; then # Check if the exome filter is true.
        ln -s "./data_except_vep/filter_exonic/${params.ref_ver}.bed" ref_filter.bed
    else
        touch ref_filter.bed
    fi

    # VEP
    if [[ "${params.ref_ver}" == "hg19" ]]; then
        vep_gnomad_name="gnomad.genomes.r2.1.sites.grch37_noVEP.vcf.gz"
        vep_cadd_name="hg19_whole_genome_SNVs.tsv.gz"
        vep_dbnsfp_name="dbNSFP4.3a_grch37.gz"
        vep_revel_name="new_tabbed_revel.tsv.gz"
    else
        vep_gnomad_name="gnomad.genomes.GRCh38.v3.1.2.sites.vcf.gz"
        vep_cadd_name="hg38_whole_genome_SNV.tsv.gz"
        vep_dbnsfp_name="dbNSFP4.1a_grch38.gz"
        vep_revel_name="new_tabbed_revel_grch38.tsv.gz"
    fi

    ln -s "./data_only_vep/vep/${params.ref_ver}/\${vep_gnomad_name}" .
    ln -s "./data_only_vep/vep/${params.ref_ver}/\${vep_cadd_name}" .
    ln -s "./data_only_vep/vep/${params.ref_ver}/\${vep_dbnsfp_name}" .
    ln -s "./data_only_vep/vep/${params.ref_ver}/\${vep_revel_name}" .

    """
}
