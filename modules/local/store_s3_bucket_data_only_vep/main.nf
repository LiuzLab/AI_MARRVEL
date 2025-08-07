process STORE_S3_BUCKET_DATA_ONLY_VEP {
    storeDir "./${params.storedir}/s3_bucket_data_only_vep/v0"

    input:
    val s3_bucket_data_name

    output:
    path "vep"

    // NOTE(JL): The path may be changed to s3://aim-data-dependencies-public-vep-v0 in future
    """
    aws s3 sync s3://$s3_bucket_data_name/vep ./vep
    """
}
