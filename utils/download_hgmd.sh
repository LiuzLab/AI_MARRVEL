#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <data_directory_prefix> <s3_bucket_id>"
    exit 1
fi

# Assign arguments to variables
DATA_DIR_PREFIX=$1
S3_BUCKET_ID=$2

# Validate the S3 bucket ID
if [[ "$S3_BUCKET_ID" != "aim-data-dependencies-2.0" && "$S3_BUCKET_ID" != "aim-data-dependencies-2.0-public" ]]; then
    echo "Error: S3 bucket ID must be either 'aim-data-dependencies-2.0' or 'aim-data-dependencies-2.0-public'."
    exit 1
fi


# Define the list of files to download
FILES=(
    "merge_expand/hg19/hgmd_c.tsv.gz"
    "merge_expand/hg19/hgmd_nc.tsv.gz"
    "merge_expand/hg38/hgmd_c.tsv.gz"
    "merge_expand/hg38/hgmd_nc.tsv.gz"
    "omim_annotate/hg19/HGMD_phen.tsv"
    "omim_annotate/hg38/HGMD_phen.tsv"
    "vep/hg19/HGMD_Pro_2022.2_hg19.vcf.gz"
    "vep/hg19/HGMD_Pro_2022.2_hg19.vcf.gz.tbi"
    "vep/hg38/HGMD_Pro_2022.2_hg38.vcf.gz"
    "vep/hg38/HGMD_Pro_2022.2_hg38.vcf.gz.tbi"
)

# Loop through each file and download it from S3
for FILE in "${FILES[@]}"; do
    # Construct the full S3 path and the local file path
    S3_PATH="s3://${S3_BUCKET_ID}/${FILE}"
    LOCAL_PATH="${DATA_DIR_PREFIX}/${FILE}"

    # Create the directory if it does not exist
    mkdir -p "$(dirname "${LOCAL_PATH}")"

    # Download the file
    echo "Downloading ${S3_PATH} to ${LOCAL_PATH}..."
    aws s3 cp "${S3_PATH}" "${LOCAL_PATH}"

    if [ $? -ne 0 ]; then
        echo "Failed to download ${S3_PATH}"
        exit 1
    fi
done

echo "All files downloaded successfully."