
process GENERATE_GENERANK {
    publishDir "${params.outdir}/${params.run_id}/prediction/generank", mode: "copy"

    input:
    path final_matrix_expanded

    output:
    path "${params.run_id}.genedf.csv"

    """
    #!/usr/bin/env python

    import pandas as pd
    import ast

    def get_transcriptdf():
        df = pd.read_csv(f"$final_matrix_expanded", low_memory=False)
        df = df.rename(columns={"Unnamed: 0": "varId"})
        df = df[df.geneSymbol != "-"]

        df["IMPACT_text"] = df["IMPACT_text"].astype(pd.CategoricalDtype(categories=["MODIFIER", "LOW", "MODERATE", "HIGH"], ordered=True))

        df["clinvarCondition"] = df["clinvarCondition"].fillna("-")
        df["Consequence"] = df["Consequence"].str.split(",")
        df[["phenoList", "phenoInhList"]] = df[["phenoList", "phenoInhList"]].applymap(lambda x: ast.literal_eval(x))
        df["zyg"] = df.zyg.apply(lambda x: "Heterozygous" if x == 1 else "Homozygous" if x == 2 else "Unknown")

        df["phenoDict"] = df[["phenoList", "phenoInhList"]].apply(lambda x: {k: v for k, v in zip(x.iloc[0], x.iloc[1])}, axis=1)
        df["phenoCount"] = df.phenoList.apply(len)

        return df

    def get_genedf(df):
        df = df.copy()
        vardf = df.sort_values(
                ['varId', 'IMPACT_text', 'phenoCount'],
                ascending=[True, False, False],
            ) \
            .drop_duplicates('varId', keep='first') \
            .sort_values('varId')

        genedf1 = vardf.groupby("geneSymbol")[["geneEnsId", "varId", "clinvarSignDesc", "clinvarCondition", "CADD_PHRED", "zyg", "gnomadAF", "predict"]].agg(lambda x: list(set(x) - {'-'}))
        genedf1["isTransHeterozygote"] = genedf1.varId.apply(lambda x: "Yes" if len(x) > 1 else "No")
        genedf2 = df.groupby("geneSymbol")[["IMPACT_text"]].max()
        genedf3 = df.groupby("geneSymbol")[["Consequence"]].agg(lambda x: list(set(xxx for xx in x for xxx in xx if xxx is not None)))
        genedf4 = df.groupby("geneSymbol")[["phenoDict"]].agg(lambda x: {k: (v or "Unknown") for phenoDict in x for k, v in phenoDict.items()})
        genedf = genedf1.join(genedf2).join(genedf3).join(genedf4)
        genedf["predict"] = genedf.predict.apply(max)
        genedf["geneEnsId"] = genedf.geneEnsId.apply(max)
        genedf = genedf.sort_values(["predict", "IMPACT_text"], ascending=False)
        genedf = genedf.reset_index()

        return genedf

    if __name__ == "__main__":
        transcriptdf = get_transcriptdf()
        genedf = get_genedf(transcriptdf)

        genedf.to_csv("./${params.run_id}.genedf.csv", index=False)
    """
}
