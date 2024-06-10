import pandas as pd
from pandas.api.types import is_numeric_dtype
import numpy as np

def comp(file_a, file_b, sort = False, skipped_col = []):
    print(file_a)
    print(file_b)
    df_a = pd.read_csv(file_a) if "csv" in file_a else pd.read_csv(file_a, sep="\t")
    df_b = pd.read_csv(file_b) if "csv" in file_b else pd.read_csv(file_b, sep="\t")

    if sort:
        def sort_df(df):
            return df.sort_values(by=df.columns.tolist()).reset_index(drop=True)

        df_a = sort_df(df_a)
        df_b = sort_df(df_b)

    assert df_a.shape == df_b.shape, f"The dimensions has to be the same. {df_a.shape} {df_b.shape}"
    assert all(df_a.columns == df_b.columns), "The columns have to be identical each other."

    for i, (c, t) in enumerate(zip(df_a.columns, df_a.dtypes)):
        if c in skipped_col: 
            continue
        
        va = df_a[c].values
        vb = df_b[c].values

        if is_numeric_dtype(t):
            diff = np.isclose(va, vb)
        else:
            diff = va == vb

        if not all(diff):
            if all(pd.isnull(va[~diff])) and all(pd.isnull(vb[~diff])):
                diff[~diff] = True

        idx = np.where(~diff)

        assert all(diff), f"The {c} column ({'numeric' if is_numeric_dtype(t) else 'non-numeric'}) has {sum(~diff)} mismatches.\n{idx}\n{va[~diff]}\n{vb[~diff]}"

    print("Two files are identical!") 
    return None




file_a = "work/ee/14e106eee477c2b81ca0179455e845/Tier.v2.tsv"
file_b = "/home/hwan/AI_MARRVEL/test_data/output_all_latest_docker/tier-test-false/12345_Tier.v2.tsv"
comp(file_a, file_b, False)

file_a = "work/ee/14e106eee477c2b81ca0179455e845/1.matrix.txt"
file_b = "/home/hwan/AI_MARRVEL/test_data/output_all_latest_docker_check_simple_repeat/12345.matrix.txt"

comp(file_a, file_b, False, ["diffuse_Phrank_STRING", "identifier"])

file_a = "work/b7/98771d7a080c1794613a028c5e8a3d/final_matrix_expanded/1.expanded.csv.gz"
file_b = "/home/hwan/AI_MARRVEL/test_data/output_all_latest_docker_check_simple_repeat/final_matrix_expanded/12345.expanded.csv.gz"

comp(file_a, file_b, False, ["diffuse_Phrank_STRING", "identifier"])



file_a = "work/b7/98771d7a080c1794613a028c5e8a3d/1.csv"
file_b = "/home/hwan/AI_MARRVEL/test_data/output_all_latest_docker_check_simple_repeat/final_matrix/12345.csv"

comp(file_a, file_b, False, ["diffuse_Phrank_STRING", "identifier"])



file_a = "work/2d/f671cde44671f7989140c3772a841f/conf_4Model/1_nd_predictions.csv"
file_b = "/home/hwan/AI_MARRVEL/test_data/output_all_latest_docker_check_simple_repeat/12345_nd_predictions.csv"

comp(file_a, file_b, True, ["diffuse_Phrank_STRING", "diffuse_Phrank_STRING_1", "diffuse_Phrank_STRING_2", "identifier"])

file_a = "work/2d/f671cde44671f7989140c3772a841f/conf_4Model/1_recessive_predictions.csv"
file_b = "/home/hwan/AI_MARRVEL/test_data/output_all_latest_docker_check_simple_repeat/12345_recessive_predictions.csv"

comp(file_a, file_b, True, ["diffuse_Phrank_STRING", "diffuse_Phrank_STRING_1", "diffuse_Phrank_STRING_2", "identifier"])

file_a = "work/2d/f671cde44671f7989140c3772a841f/conf_4Model/1_nd_recessive_predictions.csv"
file_b = "/home/hwan/AI_MARRVEL/test_data/output_all_latest_docker_check_simple_repeat/12345_nd_recessive_predictions.csv"

comp(file_a, file_b, True, ["diffuse_Phrank_STRING", "diffuse_Phrank_STRING_1", "diffuse_Phrank_STRING_2", "identifier"])

file_a = "work/2d/f671cde44671f7989140c3772a841f/conf_4Model/integrated/1_integrated.csv"
file_b = "/home/hwan/AI_MARRVEL/test_data/output_all_latest_docker_check_simple_repeat/12345_integrated.csv"

comp(file_a, file_b, True, ["diffuse_Phrank_STRING", "diffuse_Phrank_STRING_1", "diffuse_Phrank_STRING_2", "identifier"])
