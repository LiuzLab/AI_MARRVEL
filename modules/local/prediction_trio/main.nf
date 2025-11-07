process PREDICTION_TRIO {
    container "zhandongliulab/aim-lite-oldpython"
    publishDir "${params.outdir}/${params.run_id}/prediction_trio", mode: 'copy'

    input:
    path merged_compressed_scores
    path triomatrix

    output:
    path "*.trio.expanded.csv.gz"
    path "*.trio.prediction.csv"

    script:
    """
    run_final.py \
        ${projectDir}/resources/trio_pipeline/rf_trio_base.job \
        ${projectDir}/resources/trio_pipeline/features.csv \
        $triomatrix \
        ./${params.run_id}.trio.csv
    run_final.py \
        ${projectDir}/resources/trio_pipeline/rf_trio_NDG.job \
        ${projectDir}/resources/trio_pipeline/features_NDG.csv \
        $triomatrix \
        ./${params.run_id}.trio.NDG.csv


    # Generate ${params.run_id}.trio.prediction.csv
    trio_merge_rankmatrix.py ${params.run_id}

    merge_rankmatrix.py \
        ${params.run_id}.trio.prediction.csv \
        $merged_compressed_scores \
        ${params.run_id}.trio.expanded.csv.gz
    """
}
