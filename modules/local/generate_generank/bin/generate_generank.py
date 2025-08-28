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


def get_transcriptdf(filename):
    transcriptdf = pd.read_csv(filename, low_memory=False)
    transcriptdf = transcriptdf.rename(columns={"Unnamed: 0": "varId"})
    transcriptdf = transcriptdf[transcriptdf.geneSymbol != "-"]

    transcriptdf["clinvarCondition"] = \
        transcriptdf["clinvarCondition"].fillna("-")
    transcriptdf["Consequence"] = \
        transcriptdf["Consequence"].str.split(",")

    transcriptdf["IMPACT_text"] = \
        transcriptdf["IMPACT_text"].astype(IMPACT_DTYPE)
    transcriptdf["zyg"] = \
        transcriptdf.zyg.map(ZYG_MAP).fillna("Unknown")

    transcriptdf["phenoList"] = \
        transcriptdf.phenoList.map(ast.literal_eval)
    transcriptdf["phenoInhList"] = \
        transcriptdf.phenoInhList.map(ast.literal_eval)
    transcriptdf["phenoDict"] = \
        transcriptdf[["phenoList", "phenoInhList"]].apply(
            lambda x: {k: v for k, v in zip(x.iloc[0], x.iloc[1])},
            axis=1
        )
    transcriptdf["phenoCount"] = \
        transcriptdf.phenoList.apply(len)

    return transcriptdf


def get_genedf(transcriptdf):
    transcriptdf = transcriptdf.copy()

    # Extract one variant per gene with the highest impact and phenotype count
    vardf = transcriptdf.sort_values(
            ['varId', 'IMPACT_text', 'phenoCount'],
            ascending=[True, False, False],
        ) \
        .drop_duplicates('varId', keep='first') \
        .sort_values('varId')

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
    genedf["isTransHeterozygote"] = \
        genedf.varId.apply(lambda x: "Yes" if len(x) > 1 else "No")

    genedf = genedf.sort_values(["predict", "IMPACT_text"], ascending=False)
    genedf = genedf.reset_index()

    return genedf


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        required=True,
        help="Input transcript-level CSV file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output gene-level CSV file"
    )
    args = parser.parse_args()

    transcriptdf = get_transcriptdf(args.input)
    genedf = get_genedf(transcriptdf)

    genedf.to_csv(args.output, index=False)
