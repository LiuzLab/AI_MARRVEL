nextflow.enable.dsl = 2

// Process to handle the VCF file
process INDEX_VCF {
    input:
    path vcf
    
    output:
    path "input.vcf.gz"
    path "input.vcf.gz.tbi"

    script:
    """
    # Determine the file type
    INPUT_VCF_TYPE="\$(file -b $vcf)"

    if echo "\${INPUT_VCF_TYPE}" | grep -q 'symbolic link to'; then
        SYM_LINK="\$(readlink -f $vcf)"
        INPUT_VCF_TYPE="\$(file -b \${SYM_LINK})"
    fi

 
    if echo "\${INPUT_VCF_TYPE}" | grep -q 'BGZF'; then
        echo "The file is in BGZF format, ready for tabix."
        cp $vcf input.vcf.gz
    elif echo "\${INPUT_VCF_TYPE}" | grep -q 'gzip compressed data'; then
        echo "GZIP format detected, converting to BGZF."
        gunzip -c $vcf | bgzip > input.vcf.gz
    elif echo "\${INPUT_VCF_TYPE}" | grep -q 'ASCII text'; then
        echo "Plain VCF file detected, compressing and indexing."
        bgzip -c $vcf > input.vcf.gz
    else
        echo "The file $vcf does not exist or is not a recognized format."
        exit 1
    fi

    tabix -p vcf input.vcf.gz
    """
}

process VCF_PRE_PROCESS {
    input:
    path vcf
    path tbi
    path chrmap

    output:
    path "${params.run_id}.filt.vcf.gz"

    script:
    """
    # Annotate with new chromosomes and preserve original coordinates in ID
    bcftools annotate --rename-chrs $chrmap -x ID $vcf -Oz -o ${params.run_id}-annot

    # Annotate with new IDs based on CHROM, POS, REF, ALT
    bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' ${params.run_id}-annot -Oz -o ${params.run_id}-add-id.vcf.gz

    #run quality filters 
    bcftools filter  ${params.run_id}-add-id.vcf.gz -i'FILTER == "PASS"' -Oz -o ${params.run_id}.filt.vcf.gz

    #check number of variants left
    variant_count=\$(bcftools view -H ${params.run_id}.filt.vcf.gz | wc -l)
    if [ "\${variant_count}" -gt 0 ]; then
        echo "Quality filtering completed successfully. Variants passing the filters: \${variant_count}"
    else
        echo "The VCF file doesn't have FILTER annotation, or all variants filtered."
        echo "Pipeline will proceed with unfiltered VCF file."
        cp ${params.run_id}-add-id.vcf.gz ${params.run_id}.filt.vcf.gz
    fi
    """
}


process REMOVE_MITO_AND_UNKOWN_CHR {
    publishDir "${params.outdir}/vcf/", mode: 'copy'
    input:
    path vcf

    output:
    path "${params.run_id}.filt.rmMT.vcf.gz"
    path "${params.run_id}.filt.rmMT.vcf.gz.tbi"


    script:
    """
    tabix -p vcf $vcf
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $vcf -o ${params.run_id}.filt.rmMT.vcf

    bgzip ${params.run_id}.filt.rmMT.vcf
    tabix ${params.run_id}.filt.rmMT.vcf.gz
    """
}


process ANNOT_PHRANK {
    input:
    path vcf

    output:
    path "${params.run_id}-var-filt.txt"

    script:
    """
    zcat $vcf | awk 'substr(\$0, 1, 1) != "#"' | cut -f1,2,4,5 | sed 's/\t/:/g' > ${params.run_id}-var.txt
    cat ${params.run_id}-var.txt | sed 's/chr//g' | sort -u > ${params.run_id}-var-filt.txt
    """
}


process ANNOT_ENSMBLE {
    input:
    path vcf
    path ref 

    output:
    path "${params.run_id}-ensmbl.txt"

    script:
    """
    location_to_gene.py $vcf $ref | \\
     sed 's/:/\\t/g' | sed 's/X\\t/23\\t/g' | sed 's/Y\\t/24\\t/g' | \\
     sed 's/MT\\t/25\\t/g' > ${params.run_id}-ensmbl.txt
    """
}


process TO_GENE_SYM {
    input:
    path ensmbl
    path ref_to_sym
    path ref_sorted_sym

    output:
    path "${params.run_id}-gene.txt"

    script:
    """
    cat $ensmbl | sort -k5,5 | join -1 5 -2 1 - $ref_to_sym  | sed 's/ /\\t/g' | cut -f2- > genesym.txt
    cat genesym.txt | cut -f5 | sort -u | join -t\$'\\t' -1 1 -2 2 - $ref_sorted_sym | cut -f2 | sort -u > ${params.run_id}-gene.txt
    """
}


process PHRANK_SCORING {
    publishDir "${params.outdir}/phrank/", mode: 'copy'

    input:
    path gene
    path hpo
    path dagfile
    path disease_annotation
    path gene_annotation
    path disease_gene

    output:
    path "${params.run_id}.phrank.txt" 

    script:
    """
    run_phrank.py \\
        $gene $hpo $dagfile $disease_annotation $gene_annotation $disease_gene > ${params.run_id}.phrank.txt
    """
}


process HPO_SIM {
    input:
    path hpo
    path omim_hgmd_phen
    path omim_obo
    path omim_genemap2
    path omim_pheno

    output:
    path "${params.run_id}-cz"
    path "${params.run_id}-dx"

    script:
    """
    phenoSim.R $hpo $omim_hgmd_phen $omim_obo $omim_genemap2 $omim_pheno \\
        ${params.run_id}-cz ${params.run_id}-dx
    """

}

process FILTER_PROBAND {
    publishDir "${params.outdir}/vcf/", mode: 'copy'
    
    input:
    path vcf
    path tbi
    path ref_gnomad_genome
    path ref_gnomad_genome_idx
    path ref_gnomad_exome
    path ref_gnomad_exome_idx

    output:
    path "${params.run_id}.filt.rmBL.vcf"

    script:
    """
    mkdir -m777 isec_tmp1
    bcftools isec -p isec_tmp1 -w 1 -Oz \\
    $vcf $ref_gnomad_genome
    # tabix isec_tmp1/0000.vcf.gz

    mkdir -m777 isec_tmp2
    bcftools isec -p isec_tmp2 -w 1 -Oz \\
    $vcf $ref_gnomad_exome
    # tabix isec_tmp2/0000.vcf.gz

    mkdir -m777 isec_tmp3
    bcftools isec -p isec_tmp3 -Ov \\
    isec_tmp1/0000.vcf.gz isec_tmp2/0000.vcf.gz

    mv isec_tmp3/0002.vcf ${params.run_id}.filt.rmBL.vcf
    """

}

process VEP_ANNOTATE {
    publishDir "${params.outdir}/vep/", mode: "copy"

    input:
    path vcf
    path vep_dir_cache
    path vep_dir_plugins
    path vep_custom_gnomad
    path vep_custom_clinvar
    path vep_custom_hgmd
    path vep_plugin_revel
    path vep_plugin_spliceai_snv
    path vep_plugin_spliceai_indel
    path vep_plugin_cadd
    path vep_plugin_dbnsfp
    path vep_idx

    output:
    path "${params.run_id}-vep.txt"

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
        --individual all --output_file ${params.run_id}-vep.txt --input_file $vcf
    """
}

process FEATURE_ENGINEERING_PART1 {
    input:
    path vep
    path hgmd_sim
    path omim_sim
    path ref_annot_dir
    // not sure why projectDir is not working

    output:
    path "${params.run_id}_scores.csv"

    script:
    """
    AIM_FREE_RAM=\$(free -g | awk 'NR==2{printf \$7}')

    split_vep_chunks.py $vep \$AIM_FREE_RAM

    while read -r INDEX LINEH LINEA LINEB
    do
        sed -n -e "\${LINEH}p" -e "\${LINEA},\${LINEB}p" $vep > vep-\${INDEX}.txt

        feature.py \\
            -patientHPOsimiOMIM $omim_sim \\
            -patientHPOsimiHGMD $hgmd_sim \\
            -varFile vep-\${INDEX}.txt \\
            -inFileType vepAnnotTab \\
            -patientFileType one \\
            -genomeRef ${params.ref_ver} \\
            -diseaseInh AD \\
            -modules curate,conserve
        
        if [ \${INDEX} -gt 1 ]; then
            sed -n "2,\$p" scores.csv > scores_\${INDEX}.csv
        else
            mv scores.csv scores_\${INDEX}.csv
        fi
    done < vep_split.txt


    for INDEX in \$(cut -d\$'\\t' -f1 vep_split.txt)
    do
        cat scores_\${INDEX}.csv
    done > ${params.run_id}_scores.csv
    """
}

process FEATURE_ENGINEERING_PART2 {
    input:
    path scores, stageAs: "scores.csv"
    path phrank
    path ref_annot_dir
    path ref_var_tier_dir
    path ref_merge_expand_dir
    path ref_mod5_diffusion_dir

    output:
    path "${params.run_id}.matrix.txt"
    path "scores.txt.gz"

    script:
    """
    VarTierDiseaseDBFalse.R ${params.ref_ver}
    generate_new_matrix_2.py ${params.run_id} ${params.ref_ver}
    """
}

process PREDICTION {
    publishDir "${params.outdir}/prediction/", mode: "copy"

    input:
    path matrix
    path scores
    path ref_predict_new_dir
    path ref_model_inputs_dir

    output:
    path "conf_4Model/*.csv"
    path "conf_4Model/integrated/*.csv"

    script:
    """
    mkdir final_matrix_expanded
    mkdir conf_4Model
    run_final.py ${params.run_id}
    merge_rm.py ${params.run_id}
    extraModel_main.py -id ${params.run_id}
    """
}


workflow { 
    INDEX_VCF(params.input_vcf)

    VCF_PRE_PROCESS(INDEX_VCF.out, params.chrmap)

    REMOVE_MITO_AND_UNKOWN_CHR(VCF_PRE_PROCESS.out)

    ANNOT_PHRANK(VCF_PRE_PROCESS.out)

    ANNOT_ENSMBLE(ANNOT_PHRANK.out, params.ref_loc)

    TO_GENE_SYM(ANNOT_ENSMBLE.out, params.ref_to_sym, params.ref_sorted_sym)

    PHRANK_SCORING( TO_GENE_SYM.out, 
                    params.input_hpo, 
                    params.phrank_dagfile,
                    params.phrank_disease_annotation,
                    params.phrank_gene_annotation,
                    params.phrank_disease_gene)

    HPO_SIM(params.input_hpo,
            params.omim_hgmd_phen,
            params.omim_obo,
            params.omim_genemap2,
            params.omim_pheno)

    FILTER_PROBAND(
        REMOVE_MITO_AND_UNKOWN_CHR.out[0],
        REMOVE_MITO_AND_UNKOWN_CHR.out[1],
        params.ref_gnomad_genome,
        params.ref_gnomad_genome_idx,
        params.ref_gnomad_exome,
        params.ref_gnomad_exome_idx
    )

    VEP_ANNOTATE(
        FILTER_PROBAND.out,
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

    FEATURE_ENGINEERING_PART1 ( // will rename it once we have analyzed/review the part
        VEP_ANNOTATE.out,
        HPO_SIM.out[0],
        HPO_SIM.out[1],
        file(params.ref_annot_dir)
    )

    FEATURE_ENGINEERING_PART2 (
        FEATURE_ENGINEERING_PART1.out[0],
        PHRANK_SCORING.out,
        file(params.ref_annot_dir),
        file(params.ref_var_tier_dir),
        file(params.ref_merge_expand_dir),
        file(params.ref_mod5_diffusion_dir)
    )

    PREDICTION( 
        FEATURE_ENGINEERING_PART2.out[0],
        FEATURE_ENGINEERING_PART2.out[1],
        file(params.ref_predict_new_dir),
        file(params.ref_model_inputs_dir)
    )
}
