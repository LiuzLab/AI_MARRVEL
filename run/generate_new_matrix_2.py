#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import sys

#from marrvel_score_recalc import *
from add_c_nc import add_c_nc
from fillna_tier import feature_engineering as fillna
from mod5_diffusion import diffusion, diffuseSample
#from labelling import add_label
from simple_repeat_anno import simple_repeat_anno

### read score file ###
score=pd.read_csv('/out/rami-test/%s_scores.csv'%(sys.argv[1]))

### add more clin and hgmd features ###
merged = add_c_nc(score, sys.argv[2])

### run diffusion module 5 ###
if os.path.getsize(f"/out/phrank/{sys.argv[1]}.txt") == 0:
    phrank_empty = True
else:
    phrank_empty = False
    mod5 = diffuseSample(sys.argv[1], merged, "/out/phrank/")

### get original coordinates ###
merged["varId"]=merged["varId"].apply(lambda x: x.split("_E")[0])

### read tier ###
tier_f=pd.read_csv("/out/tier-test-false/"+sys.argv[1]+"_Tier.v2.tsv", sep="\t")

### add phrank calculated elsewhere for final combination ###
if phrank_empty:
    merged['phrank'] = 0
else:
    phr=pd.read_csv("/out/phrank/"+sys.argv[1]+".txt", sep="\t", names=["ENSG", "phrank"])
    merged=merged.merge(phr,left_on="geneEnsId", right_on="ENSG", how="left")

### save intermediate file for merging ###
merged.to_csv("/out/scores/"+sys.argv[1]+".txt.gz", compression="gzip", sep="\t")

### run Chaozhong's feature engineering code ###
iff=fillna(merged,tier_f)

### add diffusion module data ###
if phrank_empty:
    iff.insert(loc=0, column='diffuse_Phrank_STRING', value=0)
else:
    iff.insert(loc=0, column='diffuse_Phrank_STRING', value=mod5['diffuse_Phrank_STRING'])

### remove chr 26 ###
iff= iff.loc[~iff.index.str.startswith('26')]
#iff.to_csv(sys.argv[1]+".out.txt", sep="\t")

### run chaozhong's simple repeat annotation ###
iff=simple_repeat_anno(sys.argv[1], iff, f"/run/data_dependencies/merge_expand/{sys.argv[2]}/simpleRepeats.{sys.argv[2]}.bed")

### delete intermediate files ###
os.remove("/out/"+sys.argv[1]+".bed")
os.remove("/out/"+sys.argv[1]+".sp.bed")

### save final file for prediction ###
iff.to_csv("/out/"+sys.argv[1]+".matrix.txt", sep="\t")
