import groovy.json.JsonOutput

process GENERATE_MANIFEST_JSON {
    publishDir "${params.outdir}/${params.run_id}/manifest/", mode: 'copy'

    input:
    path ref_dir

    output:
    path "manifest.json"
    path "git.hash"
    path "data_dependencies_tree.hash"
    path "data_dependencies_tree.txt"

    script:
    
    """
cat <<EOF >> manifest.json
${JsonOutput.prettyPrint(JsonOutput.toJson(workflow.manifest))}
EOF
   
    (cd $ref_dir && find . -type f -printf '%p\\t%s\\n') \\
    | sort \\
    | tee data_dependencies_tree.txt \\
    | md5sum >> data_dependencies_tree.hash

    git -C ${projectDir} rev-parse HEAD >> git.hash
    """
}
