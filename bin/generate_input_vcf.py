#!/usr/bin/env python3.8

import string
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("variant", type=str)

args = parser.parse_args()

vcf_template = string.Template("""
##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
##FILTER=<ID=PASS,Description="All filters passed">
##reference=file:///staging/human/reference/b37/b37.fa.default/reference.bin
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
$chrom	$pos	.	$ref	$alt	.	.	.	GT:AD:DP:GQ:PL	0/1:7,5:12:99:142,0,214
""".strip())

def main(variant):
    chrom, pos, ref, alt = variant.split("-")
    with open("input.vcf", "w", encoding="ascii") as text_file:
        text_file.write(vcf_template.substitute(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
        ))

if __name__ == "__main__":
    main(**vars(args))
