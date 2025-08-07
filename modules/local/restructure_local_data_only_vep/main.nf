process RESTRUCTURE_LOCAL_DATA_ONLY_VEP {
    input:
    path local_data_path

    output:
    path "vep"

    """
    ln -s $local_data_path/vep .
    """
}
