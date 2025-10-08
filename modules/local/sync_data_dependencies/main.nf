process SYNC_DATA_DEPENDENCIES {
    label 'no_container'

    input:
    val  s3_bucket_data_name
    path local_data_dependencies_path
    val  string_local_data_dependencies_path


    output:
    path local_data_dependencies_path

    script:
    """
    if ! command -v aws >/dev/null 2>&1; then
        echo "[ERROR] AWS CLI not found. Install: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html"
        exit 1
    fi

    if [ ! -e "$local_data_dependencies_path" ]; then
        echo "[ERROR] Directory not found: $local_data_dependencies_path $string_local_data_dependencies_path"
        echo "To create/sync it, run:\n"
        echo "  aws s3 sync s3://$s3_bucket_data_name $local_data_dependencies_path --no-sign-request\n"
        exit 1
    fi

    echo "[SYNC] syncing s3://$s3_bucket_data_name to $local_data_dependencies_path $string_local_data_dependencies_path"
    aws s3 sync s3://$s3_bucket_data_name $local_data_dependencies_path --only-show-errors
    """
}
