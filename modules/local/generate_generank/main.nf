
process GENERATE_GENERANK {
    publishDir "${params.outdir}/${params.run_id}/prediction/generank", mode: "copy"

    input:
    path final_matrix_expanded

    output:
    path "${params.run_id}.gene.csv"

    """
    $moduleDir/bin/generate_generank.py \
        --input ${final_matrix_expanded} \
        --output ${params.run_id}.gene.csv
    """
}
