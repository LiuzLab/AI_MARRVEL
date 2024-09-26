#!/usr/bin/env python3.8

import pandas as pd
import sys
import io

def main():
    inputvcf = sys.stdin.read()
    df = pd.read_csv(io.StringIO(inputvcf), sep="\t", comment="#", header=None)
    df = df.rename(columns={
        0: '#CHROM',
        1: 'POS',
        2: 'ID',
        3: 'REF',
        4: 'ALT',
        5: 'QUAL',
        6: 'FILTER',
        7: 'INFO',
        8: 'FORMAT',
    })
    df[["INFO"]] = "."
    df = df[df.FILTER == "PASS"]
    vcfbody = df.to_csv(sep="\t", index=False)

    vcfformatheader = "\n".join(line for line in inputvcf.splitlines() if line.startswith("##FORMAT"))
    vcfcontigheader = "\n".join(f"##contig=<ID={chr},length={gdf.POS.max()}>" for (chr, gdf) in df.groupby("#CHROM"))

    vcfheader = f"""
##fileformat=VCFv4.2
{vcfformatheader}
{vcfcontigheader}
    """.strip()

    vcfcontent = f"""
{vcfheader}
{vcfbody}
    """.strip()

    print(vcfcontent)


if __name__ == "__main__":
    main()
