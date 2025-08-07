process CONVERT_GVCF {
    container "broadinstitute/gatk"

    input:
    path vcf
    path tbi
    tuple(
        path(fasta),
        path(fasta_index),
        path(fasta_dict),
    )
    path chrmap_file

    output:
    path "${params.run_id}.nog.vcf.gz", emit: vcf
    path "${params.run_id}.nog.vcf.gz.tbi", emit: tbi


    script:
    """
    # Define input/output paths and reference genome
    reference_genome="$fasta"
    input_file="$vcf"
    output_file="${params.run_id}.nog.vcf.gz"

    # Step 0: Check for <NON_REF>
    zcat "\$input_file" | head -n 10000 | grep -q "<NON_REF>" || { echo "It's not gVCF"; cp $vcf ${params.run_id}.nog.vcf.gz; cp $tbi ${params.run_id}.nog.vcf.gz.tbi; exit 0; }

    # Step 1: Annotate and remove ID field
    echo "Step 1: Annotating and removing ID field..."
    bcftools annotate --rename-chrs "$chrmap_file" -x ID "\$input_file" -Oz -o step1.vcf.gz
    tabix step1.vcf.gz
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y step1.vcf.gz -Oz -o step1_1.vcf.gz

    # Step 2: Sort the VCF file
    echo "Step 2: Sorting the VCF file..."
    bcftools sort step1_1.vcf.gz -Oz -o step2.vcf.gz

    # Step 2.1: Index step2.vcf.gz with tabix
    echo "Indexing step2.vcf.gz with tabix..."
    tabix -p vcf step2.vcf.gz

    # Step 3: Genotype GVCFs with GATK
    echo "Step 3: Running GenotypeGVCFs with GATK..."
    gatk GenotypeGVCFs -R "\$reference_genome" -V step2.vcf.gz -O step3.vcf.gz --allow-old-rms-mapping-quality-annotation-data

    # Step 4: Filter based on quality with GATK
    echo "Step 4: Running VariantFiltration with GATK..."
    gatk VariantFiltration -V step3.vcf.gz -O step4.vcf.gz --filter-expression "QUAL < 30.0" --filter-name "LowQual"

    # Rename the final output file
    echo "Renaming the final output file..."
    mv step4.vcf.gz "\$output_file"

    # Index the final output using tabix
    echo "Indexing the final output with tabix..."
    tabix -p vcf "\$output_file"

    # Display the number of non-comment lines (ignore lines starting with #)
    lines=\$(zcat "\$output_file" | grep -v '^#' | wc -l)
    echo "File: \$output_file has \$lines lines (excluding comment lines)."

    # Clean up intermediate files
    echo "Cleaning up intermediate files..."
    rm -f step*.vcf.gz*
    """
}
