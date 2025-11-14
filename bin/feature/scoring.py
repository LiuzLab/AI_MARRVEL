from typing import List

import pandas as pd

from annotation.utils_for_symMatch import clinVarSymMatch, hgmdSymMatch, omimSymMatch
from annotation.utils_1 import getCurationScore
from annotation.utils_for_marrvel_flatfile_module_2 import getConservationScore
from annotation.marrvel_score_recalc import (
    load_raw_matrix,
    omimCurate,
    hgmdCurate,
    clinvarCurate,
    conservationCurate,
)


def apply_curate_scores(
    annotate_info_df: pd.DataFrame,
    omim_hpo_score_df: pd.DataFrame,
    hgmd_hpo_score_acc_df: pd.DataFrame,
    hgmd_hpo_score_gene_df: pd.DataFrame,
    in_file_type: str,
) -> pd.DataFrame:
    rows: List[dict] = []
    for _, var_obj in annotate_info_df.iterrows():
        omimSymMatch(var_obj, omim_hpo_score_df, in_file_type)
        hgmdSymMatch(var_obj, hgmd_hpo_score_acc_df, hgmd_hpo_score_gene_df)
        clinVarSymMatch(var_obj, in_file_type)
        curation_scores = getCurationScore(var_obj)
        rows.append(
            {
                "curationScoreOMIM": curation_scores[0],
                "curationScoreHGMD": curation_scores[1],
                "curationScoreClinVar": curation_scores[2],
                "curationScoreTotal": curation_scores[3],
                "symptomScore": var_obj.symptomScore,
                "symptomName": var_obj.symptomName,
                "omimSymMatchFlag": var_obj.omimSymMatchFlag,
                "omimSymptomSimScore": var_obj.omimSymptomSimScore,
                "hgmdSymptomScore": var_obj.hgmdSymptomScore,
                "hgmdSymMatchFlag": var_obj.hgmdSymMatchFlag,
                "hgmdSymptomSimScore": var_obj.hgmdSymptomSimScore,
            }
        )

    result_df = pd.DataFrame(rows)
    annotate_info_df[result_df.columns] = result_df
    return annotate_info_df


def apply_conservation_scores(annotate_info_df: pd.DataFrame, disease_inh: str) -> pd.DataFrame:
    rows: List[dict] = []
    for _, var_obj in annotate_info_df.iterrows():
        ret_list = getConservationScore(var_obj, disease_inh)
        rows.append(
            {
                "conservationScoreGnomad": ret_list[0],
                "conservationScoreDGV": ret_list[1],
                "conservationScoreOELof": ret_list[2],
            }
        )

    result_df = pd.DataFrame(rows)
    annotate_info_df[result_df.columns] = result_df
    return annotate_info_df


def recalculate_scores(annotate_info_df: pd.DataFrame) -> pd.DataFrame:
    score = load_raw_matrix(annotate_info_df)
    score = omimCurate(score)
    score = hgmdCurate(score)
    score = clinvarCurate(score)
    score = conservationCurate(score)
    return score
