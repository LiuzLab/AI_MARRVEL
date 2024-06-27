import os
import itertools
import pandas as pd
import numpy as np
from glob import glob
from tqdm import tqdm, trange
import logging
from os.path import exists
import shutil

def process_sample(data_folder, sample_id, default_pred, labeling=False):
    if labeling:
        raise 'The below code was not tested with real data with labeling=True'

    recessive_folder = f"{data_folder}/recessive_matrix/"
    if not os.path.exists(recessive_folder):
        os.mkdir(recessive_folder)

    tmp_folder = f"{data_folder}/recessive_matrix/tmp/"
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)

    # read feature matrix for single var
    feature_fn = f"{sample_id}.csv"
    # feature_df = []
    # for feature_fn in feature_fns:
    if "csv" in feature_fn:
        feature_df = pd.read_csv(feature_fn, index_col=0, sep=",")
    else:
        feature_df = pd.read_csv(feature_fn, index_col=0, sep="\t")
        # feature_df.append(df)
    # feature_df = pd.concat(feature_df)
    feature_df = feature_df.loc[~feature_df.index.duplicated(keep="first")]

    # Group variants by Ens IDs
    expanded_fn = f"final_matrix_expanded/{sample_id}.expanded.csv.gz"
    # for expanded_fn in expanded_fns:
    if "csv" in expanded_fn:
        df = pd.read_csv(expanded_fn, sep=",", index_col=0, compression="infer")
    else:
        df = pd.read_csv(expanded_fn, sep="\t", index_col=0, compression="infer")

    if df.shape[0] > 100000:
        df = df.loc[df["IMPACT.from.Tier"] > 1]

    varId_geneEnsId_df = df[['geneEnsId']]
    varId_geneEnsId_df.index.name = 'varId'

    geneEnsId_varId_df = varId_geneEnsId_df[varId_geneEnsId_df.geneEnsId.str.startswith("ENSG")] \
        .reset_index() \
        .sort_values(['geneEnsId', 'varId']) \
        .drop_duplicates()

    # Generate variant pairs for each gene
    gene_var_pairs = []
    for geneEnsId, grouped_df in geneEnsId_varId_df.groupby('geneEnsId'):
        # Choose upto 6 the most meaningful variants.
        if grouped_df.shape[0] > 6:
            grouped_df = grouped_df \
                .join(default_pred, on='varId') \
                .sort_values(["predict", "IMPACT.from.Tier"], ascending=[False, False], kind="stable") \
                .head(6)

        gene_var_pairs += [
            {"geneEnsId": geneEnsId, "varId1": pair[0][1].varId, "varId2": pair[1][1].varId}
            for pair in itertools.product(grouped_df.iterrows(), grouped_df.iterrows())
        ]

    # Remove all the duplicated variant pairs.
    gene_var_pairs_df = pd.DataFrame(gene_var_pairs)
    gene_var_pairs_df = gene_var_pairs_df.drop_duplicates(['varId1', 'varId2'])

    # Use only subset columns of features
    subset_feature_names = "diffuse_Phrank_STRING,hgmdSymptomScore,omimSymMatchFlag,hgmdSymMatchFlag,clinVarSymMatchFlag,omimGeneFound,omimVarFound,hgmdGeneFound,hgmdVarFound,clinVarVarFound,clinVarGeneFound,clinvarNumP,clinvarNumLP,clinvarNumLB,clinvarNumB,dgvVarFound,decipherVarFound,curationScoreHGMD,curationScoreOMIM,curationScoreClinVar,conservationScoreDGV,omimSymptomSimScore,hgmdSymptomSimScore,GERPpp_RS,gnomadAF,gnomadAFg,LRT_score,LRT_Omega,phyloP100way_vertebrate,gnomadGeneZscore,gnomadGenePLI,gnomadGeneOELof,gnomadGeneOELofUpper,IMPACT,CADD_phred,CADD_PHRED,DANN_score,REVEL_score,fathmm_MKL_coding_score,conservationScoreGnomad,conservationScoreOELof,Polyphen2_HDIV_score,Polyphen2_HVAR_score,SIFT_score,zyg,FATHMM_score,M_CAP_score,MutationAssessor_score,ESP6500_AA_AF,ESP6500_EA_AF,hom,hgmd_rs,spliceAImax,nc_ClinVar_Exp,nc_HGMD_Exp,nc_isPLP,nc_isBLB,c_isPLP,c_isBLB,nc_CLNREVSTAT,c_CLNREVSTAT,nc_RANKSCORE,c_RANKSCORE,CLASS,phrank,isB/LB,isP/LP,cons_transcript_ablation,cons_splice_acceptor_variant,cons_splice_donor_variant,cons_stop_gained,cons_frameshift_variant,cons_stop_lost,cons_start_lost,cons_transcript_amplification,cons_inframe_insertion,cons_inframe_deletion,cons_missense_variant,cons_protein_altering_variant,cons_splice_region_variant,cons_splice_donor_5th_base_variant,cons_splice_donor_region_variant,c_ClinVar_Exp_Del_to_Missense,c_ClinVar_Exp_Different_pChange,c_ClinVar_Exp_Same_pChange,c_HGMD_Exp_Del_to_Missense,c_HGMD_Exp_Different_pChange,c_HGMD_Exp_Same_pChange,c_HGMD_Exp_Stop_Loss,c_HGMD_Exp_Start_Loss,IMPACT.from.Tier,TierAD,TierAR,TierAR.adj,No.Var.HM,No.Var.H,No.Var.M,No.Var.L,AD.matched,AR.matched,recessive,dominant,simple_repeat".split(',')
    if labeling:
        subset_feature_names.append('is_strong')
    subset_feature_df = feature_df[subset_feature_names]

    # Join feature matrix for each variant pairs.
    recessive_feature_df = gene_var_pairs_df.join(
        subset_feature_df.add_suffix('_1'),
        on="varId1"
    ).join(
        subset_feature_df.add_suffix('_2'),
        on="varId2"
    )

    # Remove unnecessary var_i == var_j cases because of zygosity
    recessive_feature_df = recessive_feature_df[
        ~((recessive_feature_df.varId1 == recessive_feature_df.varId2) &
          (recessive_feature_df.zyg_1 != 2))
    ]

    # Calculate variant distance
    recessive_feature_df['var_dist'] = (
        recessive_feature_df.varId1.str.split('-').apply(lambda x: x[1]).astype(float)
        - recessive_feature_df.varId2.str.split('-').apply(lambda x: x[1]).astype(float)
    ).abs()

    # Calculate label
    if labeling:
        recessive_feature_df['is_causal'] = (
            (recessive_feature_df.varId1 == recessive_feature_df.varId2) &
            (recessive_feature_df.zyg_1 == 2) &
            (recessive_feature_df.is_strong_1 == 1)
        ) | (
            (recessive_feature_df.varId1 != recessive_feature_df.varId2) &
            (recessive_feature_df.zyg_1 == 1) &
            (recessive_feature_df.zyg_2 == 1) &
            (recessive_feature_df.is_strong_1 == 1) &
            (recessive_feature_df.is_strong_2 == 1)
        )

    # Create pair id as the legacy did
    recessive_feature_df = recessive_feature_df.set_index(recessive_feature_df.varId1 + '_' + recessive_feature_df.varId2)

    # Drop the intermediate columns
    recessive_feature_df = recessive_feature_df.drop(
        columns=['geneEnsId', 'varId1', 'varId2']
    )
    if labeling:
        recessive_feature_df = recessive_feature_df.drop(
            columns=['is_strong_1', 'is_strong_2']
        )

    # Cast all the values to float as the legacy did
    recessive_feature_df = recessive_feature_df.astype(float)

    # Sort before saving
    recessive_feature_df = recessive_feature_df.sort_index()

    recessive_feature_df.to_csv(f"{recessive_folder}/{sample_id}.csv")
