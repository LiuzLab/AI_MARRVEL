process GENERATE_INPUT_VARIANTS {
    input:
    val input_variant

    output:
    path "${params.run_id}.variants.txt", emit: variants

    """
    echo $input_variant > ${params.run_id}.variants.txt
    """
}
