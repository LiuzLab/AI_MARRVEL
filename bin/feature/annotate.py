from typing import Optional

import numpy as np
import pandas as pd

from annotation.utils_for_marrvel_flatfile_module_3 import getAnnotateInfoRows_3


def _add_dgv_annotations(annotate_info_df: pd.DataFrame, dgv_df: pd.DataFrame) -> pd.DataFrame:
    dgv_coords = dgv_df[["Chr", "Start", "Stop"]].values.astype("int32")

    def _annotate_row(var_obj: pd.Series):
        dgv_var_found = 0
        dgv_type = "-"
        dgv_subtype = "-"

        chrom_val = int(var_obj.chrom)
        start_val = int(var_obj.start)
        stop_val = int(var_obj.stop)

        mask = (dgv_coords[:, 0] == chrom_val) & (dgv_coords[:, 1] <= start_val) & (dgv_coords[:, 2] >= stop_val)
        if np.any(mask):
            match = dgv_df[mask].iloc[0]
            dgv_var_found = 1
            dgv_type = match["type"]
            dgv_subtype = match["subType"]

        return pd.Series(
            {
                "dgvDictList": [],
                "dgvTypeList": [dgv_type],
                "dgvSubtypeList": [dgv_subtype],
                "dgvVarFound": dgv_var_found,
            }
        )

    res_df = annotate_info_df.apply(_annotate_row, axis=1)
    annotate_info_df[res_df.columns] = res_df
    return annotate_info_df


def annotate_variants(
    transcript_df: pd.DataFrame,
    genome_ref: str,
    module_list: list,
    clinvar_gene_df: pd.DataFrame,
    clinvar_allele_df,
    omim_gene_sorted_df: pd.DataFrame,
    omim_allele_list,
    hgmd_hpo_gene_sorted_df: pd.DataFrame,
    decipher_sorted_df: pd.DataFrame,
    gnomad_metrics_sorted_df: pd.DataFrame,
    dgv_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    annotate_info_df = getAnnotateInfoRows_3(
        transcript_df,
        genome_ref,
        clinvar_gene_df,
        clinvar_allele_df,
        omim_gene_sorted_df,
        omim_allele_list,
        hgmd_hpo_gene_sorted_df,
        module_list,
        decipher_sorted_df,
        gnomad_metrics_sorted_df,
    )

    if "conserve" in module_list and dgv_df is not None:
        annotate_info_df = _add_dgv_annotations(annotate_info_df, dgv_df)

    return annotate_info_df

