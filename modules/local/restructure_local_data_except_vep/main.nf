process RESTRUCTURE_LOCAL_DATA_EXCEPT_VEP {
    input:
    path local_data_path

    output:
    path "*"

    """
    find $local_data_path -follow -mindepth 1 -maxdepth 1 -not -name "vep" -exec ln -s {} . \\;
    """
}
