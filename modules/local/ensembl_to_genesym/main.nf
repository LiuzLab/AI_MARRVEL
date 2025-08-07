process ENSEMBL_TO_GENESYM {
    input:
    path ensmbl
    path ensembl_to_symbol_file

    output:
    path "${params.run_id}-gene.txt", emit: gene

    script:
    """
    # Generate sorted gene to symbol file
    sort -t \$'\t' -k2,2 $ensembl_to_symbol_file > sorted_ensembl_to_symbol.txt

    cat $ensmbl | sort -k5,5 | join -1 5 -2 1 - $ensembl_to_symbol_file  | sed 's/ /\\t/g' | cut -f2- > genesym.txt
    cat genesym.txt | cut -f5 | sort -u | join -t\$'\\t' -1 1 -2 2 - sorted_ensembl_to_symbol.txt | cut -f2 | sort -u > ${params.run_id}-gene.txt
    """
}
