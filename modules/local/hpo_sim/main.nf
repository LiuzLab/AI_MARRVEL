process HPO_SIM {
    container 'zhandongliulab/aim-lite-r'

    input:
    path hpo, name: "input.hpos.txt"
    tuple(
        path(omim_hgmd_phen),
        path(omim_obo),
        path(omim_genemap2),
        path(omim_pheno),
    )

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
