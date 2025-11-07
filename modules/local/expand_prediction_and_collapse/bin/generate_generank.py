#!/usr/bin/env python

import argparse
import pandas as pd
import ast

IMPACT_DTYPE = pd.CategoricalDtype(
    categories=["MODIFIER", "LOW", "MODERATE", "HIGH"],
    ordered=True
)

ZYG_MAP = {
    1: "Heterozygous",
    2: "Homozygous",
}


def get_vardf(filename):
    if filename.endswith(".csv"):
        vardf = pd.read_csv(filename)
    else:
        vardf = pd.read_parquet(filename)

    vardf["Consequence"] = \
        vardf.Consequence.map(ast.literal_eval)

    vardf["phenoDict"] = \
        vardf.phenoDict.map(ast.literal_eval)

    return vardf


def get_genedf(vardf):
    vardf = vardf[
        (vardf.clinvarSignDesc.str.contains("Pathogenic")
            | vardf.clinvarSignDesc.str.contains("Likely_pathogenic")
            | vardf.clinvarSignDesc.str.contains("Uncertain_significance"))
        | (vardf.zyg == "Homozygous")
        | (vardf.IMPACT_text == "HIGH")
        | ((vardf.IMPACT_text == "MODERATE") | (vardf.CADD_PHRED >= 15))
    ]

    genedf_base = vardf.groupby("geneSymbol")[[
        "geneEnsId", "varId", "clinvarSignDesc",
        "clinvarCondition", "CADD_PHRED",
        "zyg", "gnomadAF", "predict"
    ]].agg(lambda x: list(set(x) - {'-'}))

    genedf_maximpact = \
        vardf.groupby("geneSymbol")[["IMPACT_text"]].max()
    genedf_flatuniqueconsequence = \
        vardf.groupby("geneSymbol")[["Consequence"]].agg(
            lambda x: list(
                set(xxx for xx in x for xxx in xx if xxx is not None)
            )
        )
    genedf_mergedphenodict = \
        vardf.groupby("geneSymbol")[["phenoDict"]].agg(
            lambda x: {
                k: (v or "Unknown")
                for phenoDict in x for k, v in phenoDict.items()
            }
        )

    genedf = genedf_base \
        .join(genedf_maximpact) \
        .join(genedf_flatuniqueconsequence) \
        .join(genedf_mergedphenodict)

    genedf["geneEnsId"] = \
        genedf.geneEnsId.apply(max)
    genedf["predict"] = \
        genedf.predict.apply(max)

    genedf = genedf.sort_values(["predict", "IMPACT_text"], ascending=False)
    genedf = genedf.reset_index()

    return genedf


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        required=True,
        help="Input variant-level CSV file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output gene-level CSV file"
    )
    args = parser.parse_args()

    vardf = get_vardf(args.input)
    genedf = get_genedf(vardf)

    genedf.to_csv(args.output, index=False)
