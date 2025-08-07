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
