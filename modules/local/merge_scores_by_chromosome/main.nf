process MERGE_SCORES_BY_CHROMOSOME {
    label "big_mem"
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
