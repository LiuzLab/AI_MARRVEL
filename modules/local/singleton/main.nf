process PHENOPACKET_TO_VARIANTS_AND_HPOS {
    input:
    path phenopacket_json

    output:
    path "${params.run_id}.variants.txt", emit: variants
    path "${params.run_id}.hpos.txt", emit: hpo

    """
    jq -r '.phenotypicFeatures[] | select(.excluded != true) | .type.id' $phenopacket_json > ${params.run_id}.hpos.txt
    jq -r '.interpretations[].diagnosis.genomicInterpretations[].variantInterpretation.variationDescriptor.vcfRecord | "\\(.chrom | sub("^chr";""))_\\(.pos)_\\(.ref)_\\(.alt)"' $phenopacket_json > ${params.run_id}.variants.unsorted.txt
    sort -t'_' -k1,1V -k2,2n -k3,3 -k4,4 ${params.run_id}.variants.unsorted.txt > ${params.run_id}.variants.txt
    """
}

process GENERATE_INPUT_VARIANTS {
    input:
    val input_variant

    output:
    path "${params.run_id}.variants.txt", emit: variants

    """
    echo $input_variant > ${params.run_id}.variants.txt
    """
}

process GENERATE_INPUT_VCF {
    input:
    val variants

    output:
    path "input.vcf.gz", emit: vcf

    """
    generate_input_vcf.py $variants
    bgzip input.vcf
    """
}

process VALIDATE_VCF {
    container 'quay.io/biocontainers/vcftools:0.1.16--pl5321hdcf5f25_11'

    input:
    path vcf

    output:
    path "$vcf", emit: vcf

    """
    vcf-validator $vcf
    """
}

// Process to handle the VCF file
process NORMALIZE_VCF {
    input:
    path vcf

    output:
    path "input.vcf.gz", emit: vcf
    path "input.vcf.gz.tbi", emit: tbi

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
    elif echo "\${INPUT_VCF_TYPE}" | grep -q -e 'ASCII text' -e 'Unicode text' -e 'Variant Call Format'; then
        echo "Plain VCF file detected, compressing and indexing."
        bgzip -c $vcf > input.vcf.gz
    else
        echo "The file $vcf does not exist or is not a recognized format."
        exit 1
    fi

    tabix -p vcf input.vcf.gz
    """
}

process FILTER_BED {
    input:
    path vcf
    path tbi
    path ref_filter_bed

    output:
    path "${params.run_id}.recode.vcf.gz", emit: vcf
    path "${params.run_id}.recode.vcf.gz.tbi", emit: tbi

    script:
    """
    if [ -f ${ref_filter_bed} ]; then
        awk '{gsub(/^chr/, ""); print}' ${ref_filter_bed} > bed
        bcftools filter --regions-file bed ${vcf} -Oz -o "${params.run_id}.recode.vcf.gz"
        tabix -p vcf "${params.run_id}.recode.vcf.gz"
    else
        cp ${vcf} "${params.run_id}.recode.vcf.gz"
        cp ${tbi} "${params.run_id}.recode.vcf.gz.tbi"
    fi
    """
}

process BUILD_REFERENCE_INDEX {
    container "broadinstitute/gatk"
    storeDir "${params.storedir}/general/reference_index/"

    output:
    path "final_${params.ref_ver}.fa", emit: fasta
    path "final_${params.ref_ver}.fa.fai", emit: fasta_index
    path "final_${params.ref_ver}.dict", emit: fasta_dict

    script:
    """
    wget --quiet http://hgdownload.soe.ucsc.edu/goldenPath/${params.ref_ver}/bigZips/${params.ref_ver}.fa.gz
    gunzip ${params.ref_ver}.fa.gz
    sed 's/>chr/>/g' ${params.ref_ver}.fa > num_prefix_${params.ref_ver}.fa
    samtools faidx num_prefix_${params.ref_ver}.fa 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M > final_${params.ref_ver}.fa
    samtools faidx final_${params.ref_ver}.fa
    gatk CreateSequenceDictionary -R final_${params.ref_ver}.fa
    """
}

process CONVERT_GVCF {
    container "broadinstitute/gatk"

    input:
    path vcf
    path tbi
    path fasta
    path fasta_index
    path fasta_dict
    path chrmap_file

    output:
    path "${params.run_id}.nog.vcf.gz", emit: vcf
    path "${params.run_id}.nog.vcf.gz.tbi", emit: tbi


    script:
    """
    # Define input/output paths and reference genome
    reference_genome="$fasta"
    input_file="$vcf"
    output_file="${params.run_id}.nog.vcf.gz"

    # Step 0: Check for <NON_REF>
    zcat "\$input_file" | head -n 10000 | grep -q "<NON_REF>" || { echo "It's not gVCF"; cp $vcf ${params.run_id}.nog.vcf.gz; cp $tbi ${params.run_id}.nog.vcf.gz.tbi; exit 0; }

    # Step 1: Annotate and remove ID field
    echo "Step 1: Annotating and removing ID field..."
    bcftools annotate --rename-chrs "$chrmap_file" -x ID "\$input_file" -Oz -o step1.vcf.gz
    tabix step1.vcf.gz
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y step1.vcf.gz -Oz -o step1_1.vcf.gz

    # Step 2: Sort the VCF file
    echo "Step 2: Sorting the VCF file..."
    bcftools sort step1_1.vcf.gz -Oz -o step2.vcf.gz

    # Step 2.1: Index step2.vcf.gz with tabix
    echo "Indexing step2.vcf.gz with tabix..."
    tabix -p vcf step2.vcf.gz

    # Step 3: Genotype GVCFs with GATK
    echo "Step 3: Running GenotypeGVCFs with GATK..."
    gatk GenotypeGVCFs -R "\$reference_genome" -V step2.vcf.gz -O step3.vcf.gz --allow-old-rms-mapping-quality-annotation-data

    # Step 4: Filter based on quality with GATK
    echo "Step 4: Running VariantFiltration with GATK..."
    gatk VariantFiltration -V step3.vcf.gz -O step4.vcf.gz --filter-expression "QUAL < 30.0" --filter-name "LowQual"

    # Rename the final output file
    echo "Renaming the final output file..."
    mv step4.vcf.gz "\$output_file"

    # Index the final output using tabix
    echo "Indexing the final output with tabix..."
    tabix -p vcf "\$output_file"

    # Display the number of non-comment lines (ignore lines starting with #)
    lines=\$(zcat "\$output_file" | grep -v '^#' | wc -l)
    echo "File: \$output_file has \$lines lines (excluding comment lines)."

    # Clean up intermediate files
    echo "Cleaning up intermediate files..."
    rm -f step*.vcf.gz*
    """
}


process FILTER_UNPASSED {
    input:
    path vcf
    path tbi
    path chrmap

    output:
    path "${params.run_id}.filt.vcf.gz", emit: vcf
    path "${params.run_id}.filt.vcf.gz.tbi", emit: tbi

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
    tabix -p vcf ${params.run_id}.filt.vcf.gz
    """
}


process FILTER_MITO_AND_UNKOWN_CHR {
    publishDir "${params.outdir}/${params.run_id}/vcf/", mode: 'copy'
    input:
    path vcf
    path tbi

    output:
    path "${params.run_id}.filt.rmMT.vcf.gz", emit: vcf
    path "${params.run_id}.filt.rmMT.vcf.gz.tbi", emit: tbi


    script:
    """
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $vcf -o ${params.run_id}.filt.rmMT.vcf

    bgzip ${params.run_id}.filt.rmMT.vcf
    tabix ${params.run_id}.filt.rmMT.vcf.gz
    """
}


process VCF_TO_VARIANTS {
    input:
    path vcf

    output:
    path "${params.run_id}-var-filt.txt", emit: var

    script:
    """
    zcat $vcf | awk 'substr(\$0, 1, 1) != "#"' | cut -f1,2,4,5 | sed 's/\t/:/g' > ${params.run_id}-var.txt
    cat ${params.run_id}-var.txt | sed 's/chr//g' | sort -u > ${params.run_id}-var-filt.txt
    """
}


process VARIANTS_TO_ENSEMBL {
    input:
    path var
    path ref

    output:
    path "${params.run_id}-ensmbl.txt"

    script:
    """
    location_to_gene.py $var $ref | \\
     sed 's/:/\\t/g' | sed 's/X\\t/23\\t/g' | sed 's/Y\\t/24\\t/g' | \\
     sed 's/MT\\t/25\\t/g' > ${params.run_id}-ensmbl.txt
    """
}


process ENSEMBL_TO_GENESYM {
    input:
    path ensmbl
    path ref_to_sym

    output:
    path "${params.run_id}-gene.txt", emit: gene

    script:
    """
    # Generate sorted gene to symbol file
    sort -t \$'\t' -k2,2 $ref_to_sym > sorted_ref_to_sym.txt

    cat $ensmbl | sort -k5,5 | join -1 5 -2 1 - $ref_to_sym  | sed 's/ /\\t/g' | cut -f2- > genesym.txt
    cat genesym.txt | cut -f5 | sort -u | join -t\$'\\t' -1 1 -2 2 - sorted_ref_to_sym.txt | cut -f2 | sort -u > ${params.run_id}-gene.txt
    """
}


process GENESYM_TO_PHRANK {
    publishDir "${params.outdir}/${params.run_id}/phrank/", mode: 'copy'

    input:
    path gene
    path hpo
    path dagfile
    path disease_annotation
    path gene_annotation
    path disease_gene

    output:
    path "${params.run_id}.phrank.txt", emit: phrank

    script:
    """
    run_phrank.py \\
        $gene $hpo $dagfile $disease_annotation $gene_annotation $disease_gene > ${params.run_id}.phrank.txt
    """
}


process HPO_SIM {
    container 'zhandongliulab/aim-lite-r'

    input:
    path hpo, name: "input.hpos.txt"
    path omim_hgmd_phen
    path omim_obo
    path omim_genemap2
    path omim_pheno

    output:
    path "${params.run_id}-cz", emit: hgmd_sim
    path "${params.run_id}-dx", emit: omim_sim

    script:
    """
    cp $hpo input.copied.hpos.txt
    if [[ -z \$(egrep 'HP:[0-9]{7}' $hpo) ]] ; then
        echo "HP:0000001" > input.copied.hpos.txt
    fi

    phenoSim.R input.copied.hpos.txt $omim_hgmd_phen $omim_obo $omim_genemap2 $omim_pheno \\
        ${params.run_id}-cz ${params.run_id}-dx
    """

}

process FILTER_PROBAND {
    publishDir "${params.outdir}/${params.run_id}/vcf/", mode: 'copy'

    input:
    path vcf
    path tbi
    path ref_gnomad_genome
    path ref_gnomad_genome_idx
    path ref_gnomad_exome
    path ref_gnomad_exome_idx

    output:
    path "${params.run_id}.filt.rmBL.vcf.gz", emit: vcf

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
    bgzip ${params.run_id}.filt.rmBL.vcf
    """

}

process SPLIT_VCF_BY_CHROMOSOME {
    input:
    path vcf

    output:
    path "chr*.vcf.gz", emit: chr_vcfs

    script:
    """
    # Get the list of chromosomes from the VCF file

    bcftools index ${vcf}

    bcftools query -f '%CHROM\n' ${vcf} | sort | uniq > chrom_list.txt
    # Split the VCF file by chromosome
    while read chrom; do
        bcftools view -r \${chrom} ${vcf} -Oz -o chr\${chrom}.vcf.gz
    done < chrom_list.txt
    """
}

process ANNOTATE_BY_VEP {
    container 'ensemblorg/ensembl-vep:release_104.3'
    tag "${vcf.simpleName}"
    publishDir "${params.outdir}/${params.run_id}/vep/", mode: "copy"

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

process ANNOTATE_BY_MODULES {
    tag "${vep.simpleName}"

    input:
    path vep
    path hgmd_sim, stageAs: "hgmd_sim.tsv"
    path omim_sim, stageAs: "omim_sim.tsv"
    path ref_annot_dir

    output:
    path "${vep.baseName}_scores.csv", emit: scores

    script:
    """
    feature.py \\
        -patientHPOsimiOMIM $omim_sim \\
        -patientHPOsimiHGMD $hgmd_sim \\
        -varFile ${vep} \\
        -inFileType vepAnnotTab \\
        -patientFileType one \\
        -genomeRef ${params.ref_ver} \\
        -diseaseInh AD \\
        -modules curate,conserve

        mv scores.csv ${vep.baseName}_scores.csv
    """
}

process ANNOTATE_TIER {
    container 'zhandongliulab/aim-lite-r'
    tag "${scores.simpleName}"

    input:
    path scores
    path phrank
    path ref_annot_dir
    path ref_var_tier_dir

    output:
    path "${scores.simpleName}_Tier.v2.tsv", emit: tier

    script:
    """
    mv $scores scores.csv
    VarTierDiseaseDBFalse.R ${params.ref_ver}
    mv Tier.v2.tsv ${scores.simpleName}_Tier.v2.tsv
    """
}

process JOIN_PHRANK {
    tag "${scores.simpleName}"

    input:
    path scores
    path phrank
    path ref_merge_expand_dir

    output:
    path "${scores.simpleName}_scores.txt.gz", emit: compressed_scores

    script:
    """
    mv $scores scores.csv
    generate_new_matrix_2.py ${params.run_id} ${params.ref_ver}
    mv scores.txt.gz  ${scores.simpleName}_scores.txt.gz
    """
}

process MERGE_SCORES_BY_CHROMOSOME {
    publishDir "${params.outdir}/${params.run_id}/merged", mode: "copy"

    input:
    path phrank
    path tier, stageAs: "*_Tier.v2.tsv"
    path compressed_scores, stageAs: "*_scores.txt.gz"

    path ref_annot_dir
    path ref_mod5_diffusion_dir
    path ref_merge_expand_dir

    output:
    path "${params.run_id}.matrix.txt", emit: merged_matrix
    path "scores.txt.gz", emit: merged_compressed_scores

    script:
    """
    # Merge matrices
    awk 'FNR==1 && NR!=1{next;}{print}' ${tier} > Tier.v2.tsv

    # Merge compressed scores
    header_written=false

    for file in ${compressed_scores}; do
        if [ "\$header_written" = false ]; then
            # Write the first file's content, including the header
            zcat \${file} | gzip > scores.txt.gz
            header_written=true
        else
            # Skip the header and append the rest of the data
            zcat \${file} | tail -n +2 | gzip >> scores.txt.gz
        fi
    done

    post_processing.py ${params.run_id} ${params.ref_ver}
    """
}

process PREDICTION {
    publishDir "${params.outdir}/${params.run_id}/prediction/", mode: "copy"

    input:
    path merged_matrix
    path merged_compressed_scores
    path ref_model_inputs_dir

    output:
    path "conf_4Model/${params.run_id}_default_predictions.csv", emit: "default_predictions"
    path "conf_4Model/${params.run_id}_recessive_predictions.csv", optional: true
    path "conf_4Model/${params.run_id}_nd_predictions.csv", optional: true
    path "conf_4Model/${params.run_id}_nd_recessive_predictions.csv", optional: true
    path "conf_4Model/integrated/*.csv"
    path "final_matrix_expanded/*.csv.gz"
    path "shap_outputs"

    script:
    """
    mkdir final_matrix_expanded
    mkdir conf_4Model

    run_final.py \
        $ref_model_inputs_dir/default/final_model.job \
        $ref_model_inputs_dir/default/features.csv \
        $merged_matrix \
        ${params.run_id}.default_prediction.csv

    # Generate final matrix expanded
    merge_rm.py \
        ${params.run_id}.default_prediction.csv \
        $merged_compressed_scores \
        final_matrix_expanded/${params.run_id}.expanded.csv.gz

    # Generate conf_4Model
    extraModel_main.py -id ${params.run_id}
    """
}
