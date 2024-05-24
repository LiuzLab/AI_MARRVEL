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
    bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' ${params.run_id}-annot -Oz -o add-id.vcf.gz

    #run quality filters 
    bcftools filter add-id.vcf.gz -i'FILTER == "PASS"' -Oz -o ${params.run_id}.filt.vcf.gz

    #check number of variants left
    variant_count=\$(bcftools view -H ${params.run_id}.filt.vcf.gz | wc -l)
    if [ "\${variant_count}" -gt 0 ]; then
        echo "Quality filtering completed successfully. Variants passing the filters: \${variant_count}"
    else
        echo "The VCF file doesn't have FILTER annotation, or all variants filtered."
        echo "Pipeline will proceed with unfiltered VCF file."
        cp  add-id.vcf.gz ${params.run_id}.filt.vcf.gz
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
    path "*-var-filt.txt"

    script:
    """
    zcat $vcf | awk 'substr(\$0, 1, 1) != "#"' | cut -f1,2,4,5 | sed 's/\t/:/g' > var.txt
    cat var.txt | sed 's/chr//g' | sort -u > ${params.run_id}-var-filt.txt
    """
}


process ANNOT_ENSMBLE {
    input:
    path ref 
    path vcf

    output:
    path "*-ensmbl.txt"

    script:
    """
    python2.7 ${workflow.projectDir}/scripts/phrank/src/location_to_gene.py $vcf $ref | \\
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
    path "*-gene.txt"

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
    python3.8 ${workflow.projectDir}/scripts/phrank/src/run_phrank.py \\
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
    path "*-cz"
    path "*-dx"

    script:
    """
    Rscript ${workflow.projectDir}/scripts/phenoSim.R $hpo $omim_hgmd_phen $omim_obo $omim_genemap2 $omim_pheno \\
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
    path "*.filt.rmBL.vcf"

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
    cpus 10
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
    path "*-vep.txt"

    script:
    def ref_assembly = (params.ref_ver == 'hg38') ? 'GRCh38' : 'GRCh37'
    """
    /opt/vep/src/ensembl-vep/vep \\
        --dir_cache ${vep_dir_cache} \\
        --dir_plugins ${vep_dir_plugins} \\
        --fork ${task.cpus} --everything --format vcf \\
        --cache --offline --tab --force_overwrite \\
        --species homo_sapiens --assembly ${ref_assembly} \\
        --custom ${vep_custom_gnomad},gnomADg,vcf,exact,0,AF,AF_popmax,controls_nhomal \\
        --custom ${vep_custom_clinvar},clinvar,vcf,exact,0,CLNREVSTAT,CLNSIG,CLNSIGCONF \\
        --custom ${vep_custom_hgmd},hgmd,vcf,exact,0,CLASS,GENE,PHEN,RANKSCORE \\
        --af_gnomad --plugin REVEL,${vep_plugin_revel},ALL \\
        --plugin SpliceAI,snv=${vep_plugin_spliceai_snv},indel=${vep_plugin_spliceai_indel},cutoff=0.5 \\
        --plugin CADD,${vep_plugin_cadd},ALL \\
        --plugin dbNSFP,${vep_plugin_dbnsfp},ALL \\
        --individual all --output_file ${params.run_id}-vep.txt --input_file $vcf
    """
}

process FEATURE_ENGINEERING {
    input:
    path vep
    path omim_sim
    path hgmd_sim
    path ref_annot_dir
    // not sure why projectDir is not working
    path script_chunking
    path script_annot

    output:
    path "*_scores_*.txt"

    script:
    """
    #AIM_FREE_RAM=\$(free -g | awk 'NR==2{printf \$7}')

    python3.8 $script_chunking $vep 256

    while read -r INDEX LINEH LINEA LINEB
    do
        sed -n -e "\${LINEH}p" -e "\${LINEA},\${LINEB}p" $vep > vep-\${INDEX}.txt

        python3.8 main.py \\
            -outPrefix r1 \\
            -patientID $vep \\
            -patientHPOsimiOMIM $omim_sim \\
            -patientHPOsimiHGMD $hgmd_sim \\
            -varFile vep-\${INDEX}.txt \\
            -inFileType vepAnnotTab \\
            -patientFileType one \\
            -genomeRef ${params.ref_ver} \\
            -diseaseInh AD \\
            -modules curate,conserve
        
        if [ \$INDEX -gt 1 ]; then
            sed -n "2,\$p" scores.csv > scores_\$INDEX.csv
            sed -n "2,\$p" r1_scores.txt > r1_scores_\$INDEX.txt
        else
            mv scores.csv scores_\$INDEX.csv
            mv r1_scores.csv r1_scores\$INDEX.csv
        fi
    done < vep_split.txt
    """
}

workflow { 
    INDEX_VCF(params.input_vcf)

    VCF_PRE_PROCESS(INDEX_VCF.out, params.chrmap)

    REMOVE_MITO_AND_UNKOWN_CHR(VCF_PRE_PROCESS.out)

    ANNOT_PHRANK(VCF_PRE_PROCESS.out)

    ANNOT_ENSMBLE(params.ref_dir, ANNOT_PHRANK.out)

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

    FEATURE_ENGINEERING(
        VEP_ANNOTATE.out,
        HPO_SIM.out[0],
        HPO_SIM.out[1],
        params.ref_annot_dir,
        params.script_chunking,
        file(params.script_annot)
    )
}
