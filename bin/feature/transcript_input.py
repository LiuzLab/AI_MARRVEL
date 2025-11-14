import time
from typing import Tuple

import pandas as pd


RENAME_MAP = {
    "GERP++_RS": "GERPpp_RS",
    "GERP++_NR": "GERPpp_NR",
    "fathmm-MKL_coding_score": "fathmm_MKL_coding_score",
    "M-CAP_score": "M_CAP_score",
}


def _count_header_rows(var_file: str) -> int:
    num_header_skip = 0
    with open(var_file, "r") as handle:
        for line in handle:
            if line.startswith("##"):
                num_header_skip += 1
            else:
                break
    return num_header_skip


def filter_low_impact_transcripts(transcript_df: pd.DataFrame) -> pd.DataFrame:
    def _filter_group(group: pd.DataFrame) -> pd.DataFrame:
        mask = group.IMPACT.isin(["HIGH", "MODERATE"])
        if mask.any():
            return group.loc[mask]
        return group

    return (
        transcript_df.groupby("#Uploaded_variation", group_keys=False)
        .apply(_filter_group)
        .reset_index(drop=True)
    )


def _rename_columns(transcript_df: pd.DataFrame) -> None:
    rename_map = {src: dst for src, dst in RENAME_MAP.items() if src in transcript_df.columns}
    if rename_map:
        transcript_df.rename(columns=rename_map, inplace=True)


def load_transcripts(var_file: str, enable_lit: bool) -> Tuple[pd.DataFrame, float, int]:
    num_header_skip = _count_header_rows(var_file)
    print("input annoatated varFile:", var_file)
    t1 = time.time()
    transcript_df = pd.read_csv(var_file, sep="\t", skiprows=num_header_skip)

    if enable_lit:
        transcript_df = filter_low_impact_transcripts(transcript_df)

    print("shape:", transcript_df.shape)
    t2 = time.time()

    _rename_columns(transcript_df)
    input_read_time = t2 - t1
    input_num_rows = len(transcript_df.index)
    return transcript_df, input_read_time, input_num_rows

