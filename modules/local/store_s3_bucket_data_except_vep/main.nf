process STORE_S3_BUCKET_DATA_EXCEPT_VEP {
    tag "Takes 5 minutes to run"
    storeDir "./${params.storedir}/s3_bucket_data_except_vep/${params.s3_bucket_data_name}"

    output:
    path "*"

    """
    aws s3 sync s3://${params.s3_bucket_data_name} . --exclude "vep/*"
    """
}
