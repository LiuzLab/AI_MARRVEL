process PHENOPACKET_TO_VARIANTS_AND_HPOS {
    input:
    path phenopacket_json

    output:
    path "${params.run_id}.variants.txt", emit: variants
    path "${params.run_id}.hpos.txt", emit: hpo

    """
    jq -r '.phenotypicFeatures[] | select(.excluded != true) | .type.id' $phenopacket_json > ${params.run_id}.hpos.txt
    jq -r '.interpretations[].diagnosis.genomicInterpretations[].variantInterpretation.variationDescriptor.vcfRecord | "\\(.chrom | sub("^chr";""))_\\(.pos)_\\(.ref)_\\(.alt)"' $phenopacket_json > ${params.run_id}.variants.unsorted.txt
    sort -t'_' -k1,1V -k2,2n -k3,3 -k4,4 ${params.run_id}.variants.unsorted.txt > ${params.run_id}.variants.txt
    """
}
