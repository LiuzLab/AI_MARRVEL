#!/usr/bin/env python3.8

import os
import sys
import pandas as pd
# from marrvel_score_recalc import *
from fillna_tier import feature_engineering as fillna
from mod5_diffusion import diffusion, diffuseSample

# from labelling import add_label
from simple_repeat_anno import simple_repeat_anno


def main():
    ### run diffusion module 5 ###
    merged = pd.read_csv("scores.txt.gz", sep="\t")
    tier_f = pd.read_csv("Tier.v2.tsv", sep="\t").sort_index()

    phrank_empty = True
    path_phrank = sys.argv[1] + ".phrank.txt"
    if os.path.getsize(path_phrank) != 0:
        mod5 = diffuseSample(sys.argv[1], merged, ".")
        phrank_empty = False

    ### run Chaozhong's feature engineering code ###
    iff = fillna(merged, tier_f)

    iff.insert(
        loc=0, column="diffuse_Phrank_STRING", value=0 if phrank_empty else mod5["diffuse_Phrank_STRING"]
    )

    ### remove chr 26 ###
    iff = iff.loc[~iff.index.str.startswith("26")]
    # iff.to_csv(sys.argv[1]+".out.txt", sep="\t")

    ### run chaozhong's simple repeat annotation ###
    iff = simple_repeat_anno(
        sys.argv[1], iff, f"merge_expand/{sys.argv[2]}/simpleRepeats.{sys.argv[2]}.bed"
    )
    iff = iff.drop(columns=["varId_dash"])

    ### delete intermediate files ###
    # os.remove(sys.argv[1]+".bed")
    # os.remove(sys.argv[1]+".sp.bed")

    ### save final file for prediction ###
    iff.to_csv(sys.argv[1] + ".matrix.txt", sep="\t")


if __name__ == "__main__":
    main()