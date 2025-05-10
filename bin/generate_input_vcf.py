#!/usr/bin/env python3.8

import string
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("variants_file", type=str)

args = parser.parse_args()

vcf_header = """
##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
""".strip()

vcf_row_template = string.Template("""
$chrom	$pos	${chrom}_${pos}_${ref}_${alt}	$ref	$alt	.	.	.	GT:AD:DP:GQ:PL	0/1:7,5:12:99:142,0,214
""".strip())

def main(variants_file):
    with open(variants_file, "r") as f:
        variants = [line.strip() for line in f if line.strip()]

    vcf_rows = []
    for variant in variants:
        chrom, pos, ref, alt = variant.split("_")
        vcf_row = vcf_row_template.substitute(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
        )
        vcf_rows.append(vcf_row)

    with open("input.vcf", "w", encoding="ascii") as text_file:
        text_file.write(vcf_header + "\n")
        text_file.write("\n".join(vcf_rows))

if __name__ == "__main__":
    main(**vars(args))
