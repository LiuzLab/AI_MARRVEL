
def showUsage() {
    if (params.help) {
        def helpFile = file(params.usage_file)  // Specify your Markdown file path here
        if (helpFile.exists()) {
            println helpFile.text
        } else {
            println """
            Sorry something went wrong, usage file not found!
            Please vist our website for more info : https://ai-marrvel.readthedocs.io/en/latest/
            """
        }
        exit 0
    }
}

def showVersion() {
    if (!params.version) {
        return
    }

    println "1.1.0"
    exit 0
}

def validateInputParams() {
    def checkPathParamMap = [
        "input_hpo": params.input_hpo,
        "ref_dir"  : params.ref_dir,
        "ref_ver"  : params.ref_ver,
        "outdir"   : params.outdir,
    ]

    checkPathParamMap.each { paramName, paramValue ->
        if (paramValue)
            return

        println("Input parameter '${paramName}' not specified or is null!")
        println("To see usage and available parameters run `nextflow run main.nf --help`")
        exit 1
    }

    if (params.input_variant && params.input_vcf) {
        println("Error: Cannot use '--input_variant' with --input_vcf'")
        println("To see usage and available parameters run `nextflow run main.nf --help`")
        exit 1
    }

    if (params.input_vcf && !params.input_vcf.endsWith(".vcf") && !params.input_vcf.endsWith(".vcf.gz")) {
        println("Error: '--input_vcf' value '${params.input_vcf}' should be a VCF file (.vcf) or (.vcf.gz)")
        println("To see usage and available parameters run `nextflow run main.nf --help`")
        exit 1
    }

    if (params.input_variant && params.input_variant.count('-') != 3) {
        println("Error: '--input_variant' should be formated like 'X-47038564-T-C'")
        println("To see usage and available parameters run `nextflow run main.nf --help`")
        exit 1
    }

    if (!file(params.ref_dir).isDirectory()) {
        println("Error: '--ref_dir' should be an directory.")
        println("To see usage and available parameters run `nextflow run main.nf --help`")
        exit 1
    }

    if (!params.ref_ver.equals("hg19") && !params.ref_ver.equals("hg38") ) {
        println("Error: '--ref_ver' value ${params.ref_ver} should be either set to 'hg19' or 'hg38'.")
        println("To see usage and available parameters run `nextflow run main.nf --help`")
        exit 1
    }
}

def addDependentParams(params) {
    params.ref_assembly = params.ref_ver == "hg19" ? "grch37" : "grch38"

    // for data dependency
    params.chrmap = "${params.ref_dir}/bcf_annotate/chrmap.txt"

    params.ref_loc = "${params.ref_dir}/phrank/${params.ref_ver}/${params.ref_assembly}_symbol_to_location.txt"
    params.ref_to_sym = "${params.ref_dir}/phrank/${params.ref_ver}/ensembl_to_symbol.txt"
    params.ref_sorted_sym = "${params.ref_dir}/phrank/${params.ref_ver}/gene_to_symbol_sorted.txt"

    // FILTER BED
    // EXONIC FILTER BED
    params.ref_exonic_filter_bed = "${params.ref_dir}/filter_exonic/${params.ref_ver}.bed"
    if (params.bed_filter) {
        params.ref_filter_bed = params.bed_filter
    } else if (params.exome_filter) {
        params.ref_filter_bed = params.ref_exonic_filter_bed
    } else {
        params.ref_filter_bed = "/dev/null"
    }

    // for phrank
    params.phrank_dagfile = "${params.ref_dir}/phrank/${params.ref_ver}/child_to_parent.txt"
    params.phrank_disease_annotation = "${params.ref_dir}/phrank/${params.ref_ver}/disease_to_pheno.txt"
    params.phrank_gene_annotation = "${params.ref_dir}/phrank/${params.ref_ver}/gene_to_phenotype.txt"
    params.phrank_disease_gene = "${params.ref_dir}/phrank/${params.ref_ver}/disease_to_gene.txt"

    // OMIM and HPO
    params.omim_hgmd_phen = "${params.ref_dir}/omim_annotate/${params.ref_ver}/HGMD_phen.tsv" 
    params.omim_obo = "${params.ref_dir}/omim_annotate/hp.obo"
    params.omim_genemap2 = "${params.ref_dir}/omim_annotate/${params.ref_ver}/genemap2_v2022.rds"
    params.omim_pheno = "${params.ref_dir}/omim_annotate/${params.ref_ver}/HPO_OMIM.tsv"

    // GNOMAD VCF
    params.ref_gnomad_genome = "${params.ref_dir}/filter_vep/${params.ref_ver}/gnomad.${params.ref_ver}.blacklist.genomes.vcf.gz"
    params.ref_gnomad_genome_idx = "${params.ref_gnomad_genome}.tbi"
    params.ref_gnomad_exome = "${params.ref_dir}/filter_vep/${params.ref_ver}/gnomad.${params.ref_ver}.blacklist.exomes.vcf.gz"
    params.ref_gnomad_exome_idx = "${params.ref_gnomad_exome}.tbi"

    // VEP 
    params.vep_dbnsfp_name = params.ref_ver == "hg19" ? "dbNSFP4.3a_grch37.gz" : "dbNSFP4.1a_grch38.gz"
    params.vep_gnomad_name = params.ref_ver == "hg19" ? "gnomad.genomes.r2.1.sites.grch37_noVEP.vcf.gz" : "gnomad.genomes.GRCh38.v3.1.2.sites.vcf.gz"
    params.vep_cadd_name = params.ref_ver == "hg19" ? "hg19_whole_genome_SNVs.tsv.gz" : "hg38_whole_genome_SNV.tsv.gz"

    params.vep_dir_cache = "${params.ref_dir}/vep/${params.ref_ver}/"
    params.vep_dir_plugins = "${params.ref_dir}/vep/${params.ref_ver}/Plugins/"
    params.vep_custom_gnomad = "${params.ref_dir}/vep/${params.ref_ver}/${params.vep_gnomad_name}" 
    params.vep_custom_clinvar = "${params.ref_dir}/vep/${params.ref_ver}/clinvar_20220730.vcf.gz"
    params.vep_custom_hgmd = "${params.ref_dir}/vep/${params.ref_ver}/HGMD_Pro_2022.2_${params.ref_ver}.vcf.gz"
    params.vep_plugin_revel = "${params.ref_dir}/vep/${params.ref_ver}/new_tabbed_revel_${params.ref_assembly}.tsv.gz" // changed for hg19
    params.vep_plugin_spliceai_snv = "${params.ref_dir}/vep/${params.ref_ver}/spliceai_scores.masked.snv.${params.ref_ver}.vcf.gz"
    params.vep_plugin_spliceai_indel = "${params.ref_dir}/vep/${params.ref_ver}/spliceai_scores.masked.indel.${params.ref_ver}.vcf.gz"
    params.vep_plugin_cadd = "${params.ref_dir}/vep/${params.ref_ver}/${params.vep_cadd_name}" // changed for hg19
    params.vep_plugin_dbnsfp = "${params.ref_dir}/vep/${params.ref_ver}/${params.vep_dbnsfp_name}"
    params.vep_idx = "${params.ref_dir}/vep/${params.ref_ver}/*.tbi"

    params.ref_annot_dir = "${params.ref_dir}/annotate"
    params.ref_var_tier_dir = "${params.ref_dir}/var_tier"
    params.ref_merge_expand_dir = "${params.ref_dir}/merge_expand"
    params.ref_mod5_diffusion_dir = "${params.ref_dir}/mod5_diffusion"
    params.ref_predict_new_dir = "${params.ref_dir}/predict_new"
    params.ref_model_inputs_dir = "${params.ref_dir}/model_inputs"

    // Documentation
    params.usage_file = "${projectDir}/docs/source/nf_usage.txt"
    params.script_chunking = "${projectDir}/scripts/split_chunks.py"
    params.script_annot = "${projectDir}/scripts/annotation/*.py"
}
