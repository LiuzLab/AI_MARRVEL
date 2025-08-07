process ANNOTATE_BY_VEP {
    container 'ensemblorg/ensembl-vep:release_104.3'
    tag "${vcf.simpleName}"
    publishDir "${params.outdir}/${params.run_id}/vep/", mode: "copy"

    input:
    path vcf
    tuple(
        path(vep_dir_cache),
        path(vep_dir_plugins),
        path(vep_custom_gnomad),
        path(vep_custom_clinvar),
        path(vep_custom_hgmd),
        path(vep_plugin_revel),
        path(vep_plugin_spliceai_snv),
        path(vep_plugin_spliceai_indel),
        path(vep_plugin_cadd),
        path(vep_plugin_dbnsfp),
        path(vep_idx),
    )

    output:
    path "${vcf.baseName}-vep.txt", emit: vep_output

    script:
    def ref_assembly = (params.ref_ver == 'hg38') ? 'GRCh38' : 'GRCh37'
    """
    /opt/vep/src/ensembl-vep/vep \\
        --dir_cache ${vep_dir_cache} \\
        --dir_plugins ${vep_dir_plugins} \\
        --fork ${task.cpus} --everything --format vcf \\
        --cache --offline --tab --force_overwrite \\
        --species homo_sapiens --assembly ${ref_assembly} \\
        --custom ${vep_custom_gnomad},gnomADg,vcf,exact,0,AF,AF_popmax,controls_nhomalt \\
        --custom ${vep_custom_clinvar},clinvar,vcf,exact,0,CLNREVSTAT,CLNSIG,CLNSIGCONF \\
        --custom ${vep_custom_hgmd},hgmd,vcf,exact,0,CLASS,GENE,PHEN,RANKSCORE \\
        --af_gnomad --plugin REVEL,${vep_plugin_revel},ALL \\
        --plugin SpliceAI,snv=${vep_plugin_spliceai_snv},indel=${vep_plugin_spliceai_indel},cutoff=0.5 \\
        --plugin CADD,${vep_plugin_cadd},ALL \\
        --plugin dbNSFP,${vep_plugin_dbnsfp},ALL \\
        --individual all --output_file ${vcf.baseName}-vep.txt --input_file $vcf \\
        --buffer_size 50
    """
}
