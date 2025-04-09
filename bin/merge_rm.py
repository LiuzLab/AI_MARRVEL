#!/usr/bin/env python3.8
import pandas as pd
import sys

feature_path = sys.argv[1]
annot_path = sys.argv[2]
output_path = sys.argv[3]

# read feature matrix and score file
df = pd.read_csv(feature_path)
annot = pd.read_csv(annot_path, sep="\t", compression="gzip")


# get columns from score file
annot = annot[
    [
        "varId",
        "varId_dash",
        "geneSymbol",
        "geneEnsId",
        "rsId",
        "HGVSc",
        "HGVSp",
        "IMPACT",
        "Consequence",
        "phenoList",
        "phenoInhList",
        "clin_code",
        "clinvarCondition",
    ]
]

# get original coordinates
annot["varId"] = annot["varId"].apply(lambda x: x.split("_E")[0])

# rename for final output
annot = annot.rename(
    columns={"varId": "origId", "IMPACT": "IMPACT_text", "clin_code": "clinvarSignDesc"}
)

# merge
test = df.merge(annot, right_on="varId_dash", left_on="Unnamed: 0", how="left")

test.to_csv(
    output_path,
    index=False,
    compression="gzip",
)
