import json
from typing import List

import pandas as pd


def load_gnomad_metrics(path: str = "annotate/anno_hg19/gnomad.v2.1.1.lof_metrics.by_gene.txt") -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_hpo_similarity(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_clinvar_gene(path: str = "annotate/anno_hg19/gene_clinvar.csv") -> pd.DataFrame:
    clinvar_gene_df = pd.read_csv(path, sep=",")
    clinvar_gene_df.sort_values("symbol", inplace=True)
    clinvar_gene_df.set_index(["symbol"], inplace=True, drop=False)
    return clinvar_gene_df


def load_omim_gene(path: str = "annotate/anno_hg19/gene_omim.json") -> pd.DataFrame:
    with open(path) as handle:
        omim_gene_list = json.load(handle)
    return pd.DataFrame(omim_gene_list)


def load_omim_alleles(path: str = "annotate/anno_hg19/omim_alleric_variants.json") -> List[dict]:
    with open(path) as handle:
        return json.load(handle)


def _normalize_chr_column(df: pd.DataFrame, column: str = "Chr") -> pd.DataFrame:
    df[column] = df[column].replace("X", 23)
    df[column] = df[column].replace("Y", 24)
    df[column] = df[column].replace("MT", 25)
    df[column] = df[column].replace("GL.*", 26, regex=True)
    df[column] = df[column].astype(int)
    return df


def load_dgv(genome_ref: str) -> pd.DataFrame:
    path = "annotate/anno_hg38/dgv.csv" if genome_ref == "hg38" else "annotate/anno_hg19/dgv.csv"
    dgv_df = pd.read_csv(path, sep=",")
    dgv_df = dgv_df.fillna(0)
    dgv_df.columns = ["Chr", "Start", "Stop"] + dgv_df.columns.tolist()[3:]
    dgv_df["Start"] = dgv_df["Start"].astype(int)
    dgv_df["Stop"] = dgv_df["Stop"].astype(int)
    dgv_df = _normalize_chr_column(dgv_df)
    return dgv_df


def load_decipher(genome_ref: str) -> pd.DataFrame:
    path = "annotate/anno_hg38/decipher.csv" if genome_ref == "hg38" else "annotate/anno_hg19/decipher.csv"
    decipher_df = pd.read_csv(path, sep=",")
    decipher_df = decipher_df.fillna(0)
    decipher_df.columns = ["Chr", "Start", "Stop"] + decipher_df.columns.tolist()[3:]
    decipher_df["Start"] = decipher_df["Start"].astype(int)
    decipher_df["Stop"] = decipher_df["Stop"].astype(int)
    decipher_df = _normalize_chr_column(decipher_df)
    return decipher_df

