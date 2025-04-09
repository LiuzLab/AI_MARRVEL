#!/usr/bin/env python3.8
# coding: utf-8

import os
import pandas as pd
import numpy as np
import sys
from add_c_nc import add_c_nc


def main():
    ### read score file ###
    score = pd.read_csv("scores.csv")

    ### add more clin and hgmd features ###
    merged = add_c_nc(score, sys.argv[2])

    path_phrank = sys.argv[1] + ".phrank.txt"

    ### get original coordinates ###
    merged["varId"] = merged["varId"].apply(lambda x: x.split("_E")[0])

    phr = pd.read_csv(path_phrank, sep="\t", names=["ENSG", "phrank"])
    merged = merged.merge(phr, left_on="geneEnsId", right_on="ENSG", how="left")

    ### save intermediate file for merging ###
    merged.to_csv("scores.txt.gz", compression="gzip", sep="\t")


if __name__ == "__main__":
    main()