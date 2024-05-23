#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np



def add_c_nc(score, ref):

    #read data files
    if ref == 'hg38':
        clin_c=pd.read_csv("/run/data_dependencies/merge_expand/hg38/clin_c.tsv.gz", sep="\t")
        clin_nc=pd.read_csv("/run/data_dependencies/merge_expand/hg38/clin_nc.tsv.gz", sep="\t")
        hgmd_nc=pd.read_csv("/run/data_dependencies/merge_expand/hg38/hgmd_nc.tsv.gz", sep="\t")
        hgmd_c=pd.read_csv("/run/data_dependencies/merge_expand/hg38/hgmd_c.tsv.gz", sep="\t")
    else:
        clin_c=pd.read_csv("/run/data_dependencies/merge_expand/hg19/clin_c.tsv.gz", sep="\t")
        clin_nc=pd.read_csv("/run/data_dependencies/merge_expand/hg19/clin_nc.tsv.gz", sep="\t")
        hgmd_nc=pd.read_csv("/run/data_dependencies/merge_expand/hg19/hgmd_nc.tsv.gz", sep="\t")
        hgmd_c=pd.read_csv("/run/data_dependencies/merge_expand/hg19/hgmd_c.tsv.gz", sep="\t")
    

    a = score.pos.values
    ac = score.chrom.values
    temp=score[["varId"]]

    bh = clin_nc.new_stop.values
    bl = clin_nc.new_start.values
    bc= clin_nc.new_chr.values

    i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh) &(ac[:, None]==bc))

    cln=pd.concat([
        temp.loc[i, :].reset_index(drop=True),
        clin_nc.loc[j, :].reset_index(drop=True)
    ], axis=1)

    # Take into account When HGMD data base is empty
    if hgmd_nc.shape[0] == 0:
        pass
        
    else:
        #merge by region
        a = score.pos.values
        ac = score.chrom.values
        temp=score[["varId"]]
        bh = hgmd_nc.new_stop.values
        bl = hgmd_nc.new_start.values
        bc= hgmd_nc.new_chr.values
        i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh) &(ac[:, None]==bc))

        hgmd=pd.concat([
            temp.loc[i, :].reset_index(drop=True),
            hgmd_nc.loc[j, :].reset_index(drop=True)
        ], axis=1)


    merged=score.merge(cln, how="left", on="varId")
    if hgmd_nc.shape[0] == 0:
        merged['nc_HGMD_Exp'] = np.NaN
        merged['nc_RANKSCORE'] = np.NaN
    else:
        merged=merged.merge(hgmd, how="left", on="varId")
    
    if hgmd_c.shape[0] == 0:
        merged['c_HGMD_Exp'] = np.NaN
        merged['c_RANKSCORE'] = np.NaN
        merged['CLASS'] = np.NaN
    else:
        merged=merged.merge(hgmd_c, how="left", left_on=["chrom", "pos", "ref", "alt"], right_on=["new_chr", "new_pos", "ref", "alt"])
    merged=merged.merge(clin_c, how="left", left_on=["chrom", "pos", "ref", "alt"], right_on=["new_chr", "new_pos", "ref", "alt"])

    return merged

