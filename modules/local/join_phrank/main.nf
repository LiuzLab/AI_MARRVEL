process JOIN_PHRANK {
    tag "${scores.simpleName}"
    label "big_mem"

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