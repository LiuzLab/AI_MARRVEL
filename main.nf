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
    input:
    path vcf

    output:
    path "${params.run_id}.filt.rmMT.vcf"


    script:
    """
    tabix -p vcf $vcf
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $vcf -o ${params.run_id}.filt.rmMT.vcf
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
}