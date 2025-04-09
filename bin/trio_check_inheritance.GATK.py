#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys


sampleID = sys.argv[1]
vcf_file = sys.argv[2]
ref_genome = sys.argv[3]
proband_id = sys.argv[4]
paternal_id = sys.argv[5]
maternal_id = sys.argv[6]


Y_autosomal = {'hg19':{'cut1':10001,
                       'cut2':2649520,
                       'cut3':59034050,
                       'cut4':59363566},
               'hg38':{'cut1':10001,
                       'cut2':2781479,
                       'cut3':56887903,
                       'cut4':57217415}}


n = 0
with open(vcf_file,'r') as F:
    for line in F:
        if line.startswith('##'):
            n += 1
        else:
            break

vcf_df = pd.read_csv(vcf_file, sep='\t', skiprows=n)
##I removed this flag because the vcfs we're dealing with don't have the PASS flag...
#vcf_df = vcf_df.loc[vcf_df['FILTER']=='PASS',:]

probandID = [coln for coln in vcf_df.columns if proband_id in coln][0]
paternalID = [coln for coln in vcf_df.columns if paternal_id in coln][0]
maternalID = [coln for coln in vcf_df.columns if maternal_id in coln][0]


# Preprocessing
vcf_df['varid'] = vcf_df['#CHROM'].astype(str) + '_' + vcf_df['POS'].astype(str) + '_' + vcf_df['REF'] + '_' + vcf_df['ALT']
vcf_df['Patient'] = vcf_df.loc[:,probandID].str.split(':').str[0].str.replace('|','/')
vcf_df['Mother'] = vcf_df.loc[:,maternalID].str.split(':').str[0].str.replace('|','/')
vcf_df['Father'] = vcf_df.loc[:,paternalID].str.split(':').str[0].str.replace('|','/')
vcf_df = vcf_df.loc[~vcf_df['Patient'].isin(['./.','0/0']),:]
vcf_df['Pattern'] = 'Left'
#print(vcf_df)
#print(vcf_df['Patient'])
#print(vcf_df['Patient'].str.split('/',expand=True))
#print(vcf_df[['Patient_0','Patient_1']])
#quit()
vcf_df[['Patient_0','Patient_1']] = vcf_df['Patient'].str.split('/',expand=True)
vcf_df[['Mother_0','Mother_1']] = vcf_df['Mother'].str.split('/',expand=True)
vcf_df[['Father_0','Father_1']] = vcf_df['Father'].str.split('/',expand=True)


# make GT like 1/2, 0/4 to several lines with only 0 or 1
'''
multivar_df = vcf_df.loc[vcf_df['ALT'].str.contains(','),:]
vcf_df = vcf_df.loc[~vcf_df['ALT'].str.contains(','),vcf_df.columns[0:16]]
multivar_df.loc[:,'ALT_list'] = multivar_df['ALT'].str.split(',')
multivar_df.loc[:,'ALT_int'] = multivar_df['ALT_list'].apply(lambda x: list(range(1,len(x)+1)))

multivar_df.loc[:,'PGT_0_list'] = multivar_df.apply(lambda x: np.isin(np.array(x.ALT_int).astype(str),x.Patient_0).astype(int), axis=1)
multivar_df.loc[:,'PGT_1_list'] = multivar_df.apply(lambda x: np.isin(np.array(x.ALT_int).astype(str),x.Patient_1).astype(int), axis=1)

multivar_df.loc[:,'FGT_0_list'] = multivar_df.apply(lambda x: np.isin(np.array(x.ALT_int).astype(str),x.Father_0).astype(int), axis=1)
multivar_df.loc[:,'FGT_1_list'] = multivar_df.apply(lambda x: np.isin(np.array(x.ALT_int).astype(str),x.Father_1).astype(int), axis=1)

multivar_df.loc[:,'MGT_0_list'] = multivar_df.apply(lambda x: np.isin(np.array(x.ALT_int).astype(str),x.Mother_0).astype(int), axis=1)
multivar_df.loc[:,'MGT_1_list'] = multivar_df.apply(lambda x: np.isin(np.array(x.ALT_int).astype(str),x.Mother_1).astype(int), axis=1)

multivar_df.loc[:,'PGT_list'] = multivar_df.apply(lambda x: ['%s/%s'%(x.PGT_0_list[i],x.PGT_1_list[i]) for i in range(len(x.ALT_int))], axis=1)
multivar_df.loc[:,'FGT_list'] = multivar_df.apply(lambda x: ['%s/%s'%(x.FGT_0_list[i],x.FGT_1_list[i]) for i in range(len(x.ALT_int))], axis=1)
multivar_df.loc[:,'MGT_list'] = multivar_df.apply(lambda x: ['%s/%s'%(x.MGT_0_list[i],x.MGT_1_list[i]) for i in range(len(x.ALT_int))], axis=1)

multivar_new = pd.DataFrame(np.ones((0,16)), columns=multivar_df.columns[0:16])

for i in range(multivar_df.shape[0]):
    n_lines = len(multivar_df.loc[multivar_df.index[i],'ALT_list'])
    multivar_new_tmp = multivar_df.iloc[np.repeat(i,n_lines),0:16]
    multivar_new_tmp['ALT'] = multivar_df.loc[multivar_df.index[i],'ALT_list']
    multivar_new_tmp['Patient'] = multivar_df.loc[multivar_df.index[i],'PGT_list']
    multivar_new_tmp['Mother'] = multivar_df.loc[multivar_df.index[i],'MGT_list']
    multivar_new_tmp['Father'] = multivar_df.loc[multivar_df.index[i],'FGT_list']
    multivar_new = pd.concat([multivar_new, multivar_new_tmp], axis=0)
    multivar_new.index = np.arange(multivar_new.shape[0])

multivar_new['varid'] = multivar_new['#CHROM'].astype(str) + '_' + multivar_new['POS'].astype(str) + '_' + multivar_new['REF'] + '_' + multivar_new['ALT']

vcf_df = pd.concat([vcf_df, multivar_new], axis=0)
'''

vcf_df.index = np.arange(vcf_df.shape[0])
vcf_df = vcf_df.loc[~vcf_df['Patient'].isin(['./.','0/0']),:]


print('Total number of variants:', vcf_df.shape[0])


# het condition
bool_M = vcf_df['Patient'].isin(['0/1','1/0']) & (vcf_df['Mother'].str.find('1') != -1) & (vcf_df['Father'].str.find('1') == -1)
bool_F = vcf_df['Patient'].isin(['0/1','1/0']) & (vcf_df['Father'].str.find('1') != -1) & (vcf_df['Mother'].str.find('1') == -1)

# more specific if one parent is 1/1
bool_UI = vcf_df['Patient'].isin(['0/1','1/0']) & (vcf_df['Father'].str.find('1') != -1) & (vcf_df['Mother'].str.find('1') != -1)
bool_M = bool_M | (bool_UI & (vcf_df['Father']!='1/1') & (vcf_df['Mother']=='1/1'))
bool_F = bool_F | (bool_UI & (vcf_df['Father']=='1/1') & (vcf_df['Mother']!='1/1'))

bool_N = vcf_df['Patient'].isin(['0/1','1/0']) & (vcf_df['Father'].str.find('1') == -1) & (vcf_df['Mother'].str.find('1') == -1)

vcf_df.loc[bool_UI,'Pattern'] = 'UI'
vcf_df.loc[bool_M,'Pattern'] = 'M'
vcf_df.loc[bool_F,'Pattern'] = 'F'
vcf_df.loc[bool_N,'Pattern'] = 'N'

# unkown cases
vcf_df.loc[(vcf_df['Father']=='./.') & (vcf_df['Mother']=='./.') & bool_N, 'Pattern'] = 'U'



# homo condition
bool_NM = vcf_df['Patient'].isin(['1/1']) & (vcf_df['Mother'].str.find('1') != -1) & (vcf_df['Father'].str.find('1') == -1)
bool_NF = vcf_df['Patient'].isin(['1/1']) & (vcf_df['Mother'].str.find('1') == -1) & (vcf_df['Father'].str.find('1') != -1)
bool_B = vcf_df['Patient'].isin(['1/1']) & (vcf_df['Mother'].str.find('1') != -1) & (vcf_df['Father'].str.find('1') != -1)
bool_N = vcf_df['Patient'].isin(['1/1']) & (vcf_df['Mother'].str.find('1') == -1) & (vcf_df['Father'].str.find('1') == -1)

vcf_df.loc[bool_NM,'Pattern'] = 'MN'
vcf_df.loc[bool_NF,'Pattern'] = 'NF'
vcf_df.loc[bool_B,'Pattern'] = 'MF'
vcf_df.loc[bool_N,'Pattern'] = 'NN'

# Unkown cases - fix: if the other parent is 0/1 1/1
vcf_df.loc[(vcf_df['Mother']=='./.') & bool_NF, 'Pattern'] = 'UF'
vcf_df.loc[(vcf_df['Father']=='./.') & bool_NM, 'Pattern'] = 'MU'
vcf_df.loc[(vcf_df['Father']=='./.') & (vcf_df['Mother']=='./.') & bool_N, 'Pattern'] = 'UU'




# Y chromosome non-psudoautosomal region
# https://en.wikipedia.org/wiki/Pseudoautosomal_region
bool_non_pseudo = (vcf_df['#CHROM'] == 'Y') & ((vcf_df['POS'] < Y_autosomal[ref_genome]['cut1']) | ((vcf_df['POS']>Y_autosomal[ref_genome]['cut2']) & (vcf_df['POS'] < Y_autosomal[ref_genome]['cut3'])) | (vcf_df['POS']>Y_autosomal[ref_genome]['cut4']))
bool_F = vcf_df['Patient'].isin(['0/1','1/0']) & (vcf_df['Father'].str.find('1') != -1) & bool_non_pseudo
bool_N = vcf_df['Patient'].isin(['0/1','1/0']) & (vcf_df['Father'].str.find('1') == -1) & bool_non_pseudo

vcf_df.loc[bool_F,'Pattern'] = 'F'
vcf_df.loc[bool_N,'Pattern'] = 'N'

bool_F = vcf_df['Patient'].isin(['1/1']) & (vcf_df['Father'].str.find('1') != -1) & bool_non_pseudo
bool_N = vcf_df['Patient'].isin(['1/1']) & (vcf_df['Father'].str.find('1') == -1) & bool_non_pseudo

vcf_df.loc[bool_F,'Pattern'] = 'F'
vcf_df.loc[bool_N,'Pattern'] = 'N'



# GATK results
bool_hidenovo = vcf_df['INFO'].str.find('hiConfDeNovo=%s'%(proband_id)) != -1
bool_lowdenovo = vcf_df['INFO'].str.find('loConfDeNovo=%s'%(proband_id)) != -1

vcf_df['GATK'] = 'Inherited'
vcf_df.loc[bool_hidenovo, 'GATK'] = 'hiConfDeNovo'
vcf_df.loc[bool_lowdenovo, 'GATK'] = 'loConfDeNovo'

inheritance = vcf_df.loc[:,['varid','Patient','Mother','Father','Pattern','GATK']].copy()
inheritance.to_csv('./%s.inheritance.txt'%(sampleID), index=False)

inheritance_summary = (inheritance['GATK'] + '_' + inheritance['Pattern']).value_counts()
inheritance_summary.to_csv('./%s.summary.txt'%(sampleID))
print((inheritance['GATK'] + '_' + inheritance['Pattern']).value_counts())
