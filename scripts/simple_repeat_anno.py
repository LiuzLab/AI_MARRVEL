import pandas as pd
import numpy as np
import os

#annotates a feature matrix with a simple repeat column

def simple_repeat_anno(sample_id, feature_file, sp_bed_file):
    sample = feature_file
    sample_bed = pd.DataFrame({'chr':sample.index.str.split('-').str[0],
                           'start':sample.index.str.split('-').str[1], 
                           'end':sample.index.str.split('-').str[1],
                           'ref':sample.index.str.split('-').str[2],
                           'alt':sample.index.str.split('-').str[3]})
    sample_bed.to_csv('/out/%s.bed'%(sample_id),index=False,header=False,sep='\t')

    # bedtools needed
    os.system('/run/bedtools intersect -a /out/%s.bed -b %s -wa > /out/%s.sp.bed'%(sample_id, sp_bed_file, sample_id))

    if os.stat('/out/%s.sp.bed'%(sample_id)).st_size == 0:
        sample['simple_repeat'] = 0
    else:
        sample_sp = pd.read_csv('/out/%s.sp.bed'%(sample_id),sep='\t',header=None)
        sample_sp['var_dash'] = sample_sp[0].astype(str) + '-' + sample_sp[2].astype(str) + '-' + sample_sp[3].astype(str) + '-' + sample_sp[4].astype(str)
        sample['simple_repeat'] = sample.index.isin(sample_sp['var_dash'].to_numpy()).astype(int)

    return sample



