import pandas as pd
import numpy as np
from glob import glob
import sys


def matchID_in_inheritance(score_file, singleton_rank, inheritance_file, mode='test'):

    score_df = pd.read_csv(score_file,
                           sep='\t', compression='gzip')
    score_df = score_df.loc[~(score_df['chrom'].isin([25,26]))]
    rank_df = pd.read_csv(singleton_rank, index_col=0)
    inh_df = pd.read_csv(inheritance_file)
    print('score_df:', score_df.shape)
    print('rank_df:', rank_df.shape)
    print('inh_df:', inh_df.shape)
    
    inh_df = inh_df.loc[~(inh_df['varid'].str.startswith('GL'))].copy()
    inh_df['varid'] = inh_df['varid'].str.replace('\.0','', regex=True)
    inh_df['varid'] = inh_df['varid'].str.replace('X','23')
    inh_df['varid'] = inh_df['varid'].str.replace('Y','24')
    
    ### To accommodate the situation where hg38 files have "chr" in the chromosome naming: Remove 'chr' prefix from varId in inheritance file ###
    inh_df['varid'] = inh_df['varid'].str.replace('^chr', '', regex=True)

    
    score_df_2 = score_df.loc[:,['varId','geneEnsId']].copy()
    score_df_2['varId'] = score_df_2['varId'].str.split('_').str[0] + '_'+\
                          score_df_2['varId'].str.split('_').str[1] + '_'+\
                          score_df_2['varId'].str.split('_').str[2] + '_'+\
                          score_df_2['varId'].str.split('_').str[3]
    score_df_2 = score_df_2.loc[~score_df_2.duplicated()]
    score_df_2 = score_df_2.groupby(by='varId', as_index=False).agg({'geneEnsId' : ','.join})  
    score_df_2['geneEnsId'] = score_df_2['geneEnsId'].str.replace('-','').str.rstrip(',')
    n_rows = score_df_2.shape[0]
    n_varid = len(np.unique(score_df_2['varId']))

    # check if varId is unique
    assert n_rows == n_varid
    score_df_2.index = score_df_2['varId']
    
    inh_df = inh_df.loc[inh_df['varid'].isin(score_df_2['varId'])].copy()
    inh_df['geneEnsId'] = score_df_2.loc[inh_df['varid'], 'geneEnsId'].tolist()

    assert np.all(inh_df['varid'].isin(rank_df.index))


    if mode == 'test':
        return inh_df

    else:
        if np.all(rank_df.loc[rank_df['is_causal']==1].index.isin(inh_df['varid'])):
            return inh_df
        else:
            sys.exit('%s: not all causal variants in inheritance file!'%(ID))


def matrix_inh(inh_df, rank_file):

    # Load in rank file
    rank_df = pd.read_csv(rank_file, index_col=0)

    if 'confidence' in rank_df.columns.tolist():
        rank_df = rank_df.drop(columns=['confidence','confidence level'])
    features = rank_df.columns[:-4].tolist() + ['ranking']
    features_name = features[:-2] + ['rf_score','rf_rank']
    rank_df = rank_df.loc[:,features].copy()
    rank_df.columns = features_name

    # Load in inheritance file
    inh_df.index = inh_df['varid']
    # feature engineering
    inh_df.loc[inh_df['GATK']=='Inherited','GATK'] = 0
    inh_df.loc[inh_df['GATK']=='loConfDeNovo','GATK'] = 1
    inh_df.loc[inh_df['GATK']=='hiConfDeNovo','GATK'] = 2
    inh_df['GATK'] = inh_df['GATK'].astype(int)
    inh_df['DeNovo'] = (inh_df['GATK']==2).astype(int)
    inh_df['Father'] = ((inh_df['Pattern'].isin(['F','MF','NF','MU','UF','UU'])) | (inh_df['Pattern'].isin(['MN','NN']) & (inh_df['GATK']!=2))).astype(int)
    inh_df['Mother'] = ((inh_df['Pattern'].isin(['M','MF','MN','MU','UF','UU'])) | (inh_df['Pattern'].isin(['NF','NN']) & (inh_df['GATK']!=2))).astype(int)
    inh_df['Unknown'] = ((inh_df['DeNovo']==0) & (inh_df['Father']==0) & (inh_df['Mother']==0)).astype(int)

    # Merge the two
    rank_df = rank_df.loc[rank_df.index.isin(inh_df.index)].copy()
    added_features = ['DeNovo','Father','Mother','Unknown','geneEnsId']
    rank_df[added_features] = inh_df[added_features]
    rank_df.loc[rank_df['geneEnsId'].isna(), 'geneEnsId'] = [str(i) for i in range(np.sum(rank_df['geneEnsId'].isna()))]

    # Gene level variant summary
    genes = list(set((','.join(rank_df['geneEnsId'].tolist())).split(',')))  ### Changed from set to list ###
    gene_vars_df = pd.DataFrame({#'Gene':genes,
                            'Gene_vars_Father':0,
                            'Gene_vars_Mother':0,
                            'Gene_vars_DeNovo':0,
                            'Gene_vars_Unknown':0}, index=genes)

    genes_HMImpact = list(set((','.join(rank_df.loc[rank_df['IMPACT.from.Tier']>=3,'geneEnsId'].tolist())).split(',')))

    for g in genes_HMImpact:
        gene_df = rank_df.loc[rank_df['geneEnsId'].str.contains(g)].copy()
        gene_df = gene_df.loc[gene_df['IMPACT.from.Tier']>=3].copy()
        gene_vars_df.loc[g, 'Gene_vars_Father'] = gene_df['Father'].sum()
        gene_vars_df.loc[g, 'Gene_vars_Mother'] = gene_df['Mother'].sum()
        gene_vars_df.loc[g, 'Gene_vars_DeNovo'] = gene_df['DeNovo'].sum()
        gene_vars_df.loc[g, 'Gene_vars_Unknown'] = gene_df['Unknown'].sum()

    rank_df['varId'] = rank_df.index
    rank_df_new = pd.DataFrame(np.ones((0,2)), columns=['varId','geneEnsId'])

    for i in range(rank_df.shape[0]):
            n_lines = len(rank_df.loc[rank_df.index[i],'geneEnsId'].split(','))
            rank_df_new_tmp = rank_df.loc[np.repeat(rank_df.index[i],n_lines),['varId','geneEnsId']]
            rank_df_new_tmp['geneEnsId'] = rank_df.loc[rank_df.index[i],'geneEnsId'].split(',')

            rank_df_new = pd.concat([rank_df_new, rank_df_new_tmp], axis=0)
            rank_df_new.index = np.arange(rank_df_new.shape[0])
            
    rank_df_new['Gene_vars_Father'] = 0
    rank_df_new['Gene_vars_Mother'] = 0
    rank_df_new['Gene_vars_DeNovo'] = 0
    rank_df_new['Gene_vars_Unknown'] = 0
    rank_df_new.loc[:,['Gene_vars_Father','Gene_vars_Mother','Gene_vars_DeNovo','Gene_vars_Unknown']] = gene_vars_df.loc[rank_df_new['geneEnsId'].tolist(), 
                                                                  gene_vars_df.columns.tolist()].to_numpy()
    rank_df_new = rank_df_new.groupby(by='varId')[['Gene_vars_Father','Gene_vars_Mother','Gene_vars_DeNovo','Gene_vars_Unknown']].max()
    rank_df_new = rank_df_new.loc[rank_df.index]

    assert np.all(rank_df.index == rank_df_new.index)

    rank_df['Gene_vars_Father'] = rank_df_new['Gene_vars_Father']
    rank_df['Gene_vars_Mother'] = rank_df_new['Gene_vars_Mother']
    rank_df['Gene_vars_DeNovo'] = rank_df_new['Gene_vars_DeNovo']
    rank_df['Gene_vars_Unknown'] = rank_df_new['Gene_vars_Unknown']

    rank_df['No.Risk.Var'] = rank_df['No.Var.H'] + rank_df['No.Var.M']
    rank_df['chrom'] = rank_df.index.str.split('_').str[0]
    rank_df['chrom'] = rank_df['chrom'].astype(int)

    rank_df['Homozygous_Recessive'] = ((rank_df['zyg']==2) & (rank_df['chrom']!=23)).astype(int)
    rank_df['De_Novo'] = ((rank_df['DeNovo']==1) & (rank_df['zyg']==1)).astype(int)
    rank_df['X_Linked'] = (rank_df['chrom']==23).astype(int)

    rank_df['Compound_heterozygous'] = (
                                        (rank_df['zyg']==1) &\
                                        (rank_df['No.Risk.Var']>0) &\
                                        (((rank_df['Father']==1) & ((rank_df['Gene_vars_Mother']+rank_df['Gene_vars_DeNovo'])>0)) |\
                                         ((rank_df['Mother']==1) & ((rank_df['Gene_vars_Father']+rank_df['Gene_vars_DeNovo'])>0)))
                                       ).astype(int)


    rank_df['Inherited_dominant'] = (
                                     (rank_df['zyg']==1) &\
                                     (rank_df['DeNovo']==0) &\
                                     (rank_df['X_Linked']==0) &\
                                     (rank_df['Compound_heterozygous']==0)
                                    ).astype(int)


    rank_df['Homozygous_Recessive_Risk'] = rank_df['Homozygous_Recessive'].copy()
    rank_df.loc[rank_df['IMPACT.from.Tier']<3,'Homozygous_Recessive_Risk'] = 0

    rank_df['De_Novo_Risk'] = rank_df['De_Novo'].copy()
    rank_df.loc[rank_df['IMPACT.from.Tier']<3,'De_Novo_Risk'] = 0

    rank_df['X_Linked_Risk'] = rank_df['X_Linked'].copy()
    rank_df.loc[rank_df['IMPACT.from.Tier']<3,'X_Linked_Risk'] = 0

    rank_df['Compound_heterozygous_Risk'] = rank_df['Compound_heterozygous'].copy()
    rank_df.loc[rank_df['IMPACT.from.Tier']<3,'Compound_heterozygous_Risk'] = 0

    rank_df['Inherited_dominant_Risk'] = rank_df['Inherited_dominant'].copy()
    rank_df.loc[rank_df['IMPACT.from.Tier']<3,'Inherited_dominant_Risk'] = 0

    return rank_df #.to_csv('matrix_inh/%s.txt'%ID)
    










