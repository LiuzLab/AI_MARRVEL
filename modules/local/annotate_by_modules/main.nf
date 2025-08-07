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
        ${params.impact_filter ? "-enableLIT" : ""} \\
        -diseaseInh AD \\
        -modules curate,conserve \\
        -inFileType vepAnnotTab \\
        -patientFileType one \\
        -patientHPOsimiOMIM ${omim_sim} \\
        -patientHPOsimiHGMD ${hgmd_sim} \\
        -varFile ${vep} \\
        -genomeRef ${params.ref_ver}

    mv scores.csv ${vep.baseName}_scores.csv
    """
}
