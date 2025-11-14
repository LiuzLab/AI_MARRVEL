import time
from typing import Dict

import pandas as pd

from feature.annotate import annotate_variants
from feature.data_loaders import (
    load_clinvar_gene,
    load_decipher,
    load_dgv,
    load_gnomad_metrics,
    load_hpo_similarity,
    load_omim_alleles,
    load_omim_gene,
)
from feature.scoring import apply_conservation_scores, apply_curate_scores, recalculate_scores
from feature.transcript_input import load_transcripts


def _prepare_sorted_frames(
    gnomad_metrics_gene_df: pd.DataFrame,
    decipher_df: pd.DataFrame,
    omim_gene_df: pd.DataFrame,
    hgmd_hpo_score_df: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    sorted_frames: Dict[str, pd.DataFrame] = {
        "decipher": pd.DataFrame(),
        "gnomad": pd.DataFrame(),
        "omim_gene": pd.DataFrame(),
        "hgmd_gene": pd.DataFrame(),
        "hgmd_acc": pd.DataFrame(),
    }

    if not decipher_df.empty:
        sorted_frames["decipher"] = decipher_df.set_index(["Chr", "Start", "Stop"]).sort_index()
    if not gnomad_metrics_gene_df.empty:
        sorted_frames["gnomad"] = gnomad_metrics_gene_df.groupby("gene").first().sort_index()
    if not omim_gene_df.empty:
        sorted_frames["omim_gene"] = omim_gene_df.set_index("geneSymbol").sort_index()
    if not hgmd_hpo_score_df.empty:
        sorted_frames["hgmd_gene"] = hgmd_hpo_score_df.groupby("gene_sym").first().sort_index()
        sorted_frames["hgmd_acc"] = hgmd_hpo_score_df.groupby("acc_num").first().sort_index()

    return sorted_frames


def _load_curate_inputs(args):
    omim_hpo_score_df = pd.DataFrame()
    hgmd_hpo_score_df = pd.DataFrame()
    clinvar_gene_df = pd.DataFrame()
    omim_gene_df = pd.DataFrame()
    omim_allele_list = []

    omim_hpo_score_df = load_hpo_similarity(args.patientHPOsimiOMIM)
    print("patientHPOsimi-OMIM dimension:", omim_hpo_score_df.shape)

    hgmd_hpo_score_df = load_hpo_similarity(args.patientHPOsimiHGMD)
    print("patientHPOsimi-HGMD dimension:", hgmd_hpo_score_df.shape)

    # legacy code uses the hg19 file for both references
    clinvar_gene_df = load_clinvar_gene()
    omim_gene_df = load_omim_gene()
    omim_allele_list = load_omim_alleles()

    return (
        omim_hpo_score_df,
        hgmd_hpo_score_df,
        clinvar_gene_df,
        omim_gene_df,
        omim_allele_list,
    )


def run_pipeline(args) -> None:
    print("input file:", args.varFile)
    print("type of input file:", args.inFileType)
    print("low impact transcripts enabled: ", args.enableLIT)

    start_time = time.time()

    gnomad_metrics_gene_df = load_gnomad_metrics()

    (
        omim_hpo_score_df,
        hgmd_hpo_score_df,
        clinvar_gene_df,
        omim_gene_df,
        omim_allele_list,
    ) = _load_curate_inputs(args)

    clinvar_allele_df = []

    print("reading DGV flat file")
    dgv_df = load_dgv(args.genomeRef)
    print("finsihed reading DGV")

    print("reading Decipher flat file")
    decipher_df = load_decipher(args.genomeRef)
    print("finsihed reading DECIPHER")

    transcript_df, input_read_time, input_num_rows = load_transcripts(args.varFile, args.enableLIT)

    sorted_frames = _prepare_sorted_frames(
        gnomad_metrics_gene_df, decipher_df, omim_gene_df, hgmd_hpo_score_df
    )

    annotate_info_df = annotate_variants(
        transcript_df,
        args.genomeRef,
        clinvar_gene_df,
        clinvar_allele_df,
        sorted_frames["omim_gene"],
        omim_allele_list,
        sorted_frames["hgmd_gene"],
        sorted_frames["decipher"],
        sorted_frames["gnomad"],
        dgv_df,
    )

    annotate_info_df = apply_curate_scores(
        annotate_info_df,
        omim_hpo_score_df,
        sorted_frames["hgmd_acc"],
        sorted_frames["hgmd_gene"],
        args.inFileType,
    )

    annotate_info_df = apply_conservation_scores(annotate_info_df, args.diseaseInh)

    end_time = time.time()
    process_time = end_time - start_time

    print("pipeline time:", process_time)
    with open("log.txt", "w") as handle:
        handle.write(
            "Process time:"
            + str(process_time)
            + " seconds and in mins:"
            + str(process_time / 60)
            + "\n"
        )
    print("log file name:", "log.txt")
    print("input read time:", input_read_time)
    print("input num rows:", input_num_rows)
    print("Score re-calculation:")
    score = recalculate_scores(annotate_info_df)
    score.to_csv("scores.csv", index=False)
