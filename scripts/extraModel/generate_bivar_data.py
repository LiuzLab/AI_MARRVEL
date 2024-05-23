import pandas as pd
import os
import numpy as np
from glob import glob
from tqdm import tqdm, trange
import logging
from multiprocessing import Pool
from os.path import exists
import shutil

def process_sample(data_folder, sample_id,
                   default_pred, labeling=False, n_thread = 10):

    # recessive_folder = f'{data_folder}/recessive_matrix'
    # if not os.path.exists(recessive_folder):
    #     os.mkdir(recessive_folder)

    recessive_folder = f'{data_folder}/recessive_matrix/'
    if not os.path.exists(recessive_folder):
        os.mkdir(recessive_folder)

    tmp_folder = f'{data_folder}/recessive_matrix/tmp/'
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)

    ## read feature matrix for single var
    feature_fn = f"/out/final_matrix/{sample_id}.csv"
    #feature_df = []
    #for feature_fn in feature_fns:
    if 'csv' in feature_fn:
        feature_df = pd.read_csv(feature_fn, index_col=0, sep=",")
    else:
        feature_df = pd.read_csv(feature_fn, index_col=0, sep="\t")
        #feature_df.append(df)
    #feature_df = pd.concat(feature_df)
    feature_df = feature_df.loc[~feature_df.index.duplicated(keep='first')]   

    ## Group variants by Ens IDs
    expanded_fn = f"/out/final_matrix_expanded/{sample_id}.expanded.csv.gz"
    gene_dict = {}
    #for expanded_fn in expanded_fns:
    if 'csv' in expanded_fn:
        df = pd.read_csv(expanded_fn, sep=",", index_col=0, compression='infer')
    else:
        df = pd.read_csv(expanded_fn, sep="\t", index_col=0, compression='infer')
        
    if df.shape[0] > 100000:
        df = df.loc[df['IMPACT.from.Tier'] > 1]

    for varId, gene in zip(df.index.tolist(), df.geneEnsId.tolist()):
        if gene.startswith('ENSG'):
            if gene not in gene_dict:
                gene_dict[gene] = {varId}
            else:
                gene_dict[gene].add(varId)
    
    ## create parameters and run bivar feature generation in parallel
    params = [{'feature_df': feature_df,
                'default_pred': default_pred,
                'labeling': labeling,
                'out_folder': f'{tmp_folder}',
                'gene': gene,
               'varIDs': list(gene_dict[gene])} for gene in gene_dict
               ]
    print(f"Now starting to generate recessive feature matrix for {len(gene_dict)} genes, {feature_df.shape[0]} variants using {n_thread} threads.")
    p = Pool(processes=n_thread)
    
    with tqdm(total=len(params)) as pbar:
        for result in p.imap_unordered(process_gene, params):
            pbar.update()
    p.close()
    p.join()
    print("Recessive features for each gene finished, now putting together...")

    bivar_feature_mats = []
    fns = glob(tmp_folder + "/*.csv")
    for fn in tqdm(fns, total=len(fns)):
        df = pd.read_csv(fn, index_col=0)
        bivar_feature_mats.append(df)
    if bivar_feature_mats:
        bivar_feature_mats = pd.concat(bivar_feature_mats)
        bivar_feature_mats.to_csv(f"{recessive_folder}/{sample_id}.csv")
        print("Recessive features saved, now removing tmp files...")
        print("Recessive features generation finished...")
    else:
        print("No recessive variant pair found. End and remove temp files")
    
    shutil.rmtree(f'{tmp_folder}')

def process_gene(param):
    
    feature_names = 'diffuse_Phrank_STRING,hgmdSymptomScore,omimSymMatchFlag,hgmdSymMatchFlag,clinVarSymMatchFlag,omimGeneFound,omimVarFound,hgmdGeneFound,hgmdVarFound,clinVarVarFound,clinVarGeneFound,clinvarNumP,clinvarNumLP,clinvarNumLB,clinvarNumB,dgvVarFound,decipherVarFound,curationScoreHGMD,curationScoreOMIM,curationScoreClinVar,conservationScoreDGV,omimSymptomSimScore,hgmdSymptomSimScore,GERPpp_RS,gnomadAF,gnomadAFg,LRT_score,LRT_Omega,phyloP100way_vertebrate,gnomadGeneZscore,gnomadGenePLI,gnomadGeneOELof,gnomadGeneOELofUpper,IMPACT,CADD_phred,CADD_PHRED,DANN_score,REVEL_score,fathmm_MKL_coding_score,conservationScoreGnomad,conservationScoreOELof,Polyphen2_HDIV_score,Polyphen2_HVAR_score,SIFT_score,zyg,FATHMM_score,M_CAP_score,MutationAssessor_score,ESP6500_AA_AF,ESP6500_EA_AF,hom,hgmd_rs,spliceAImax,nc_ClinVar_Exp,nc_HGMD_Exp,nc_isPLP,nc_isBLB,c_isPLP,c_isBLB,nc_CLNREVSTAT,c_CLNREVSTAT,nc_RANKSCORE,c_RANKSCORE,CLASS,phrank,isB/LB,isP/LP,cons_transcript_ablation,cons_splice_acceptor_variant,cons_splice_donor_variant,cons_stop_gained,cons_frameshift_variant,cons_stop_lost,cons_start_lost,cons_transcript_amplification,cons_inframe_insertion,cons_inframe_deletion,cons_missense_variant,cons_protein_altering_variant,cons_splice_region_variant,cons_splice_donor_5th_base_variant,cons_splice_donor_region_variant,c_ClinVar_Exp_Del_to_Missense,c_ClinVar_Exp_Different_pChange,c_ClinVar_Exp_Same_pChange,c_HGMD_Exp_Del_to_Missense,c_HGMD_Exp_Different_pChange,c_HGMD_Exp_Same_pChange,c_HGMD_Exp_Stop_Loss,c_HGMD_Exp_Start_Loss,IMPACT.from.Tier,TierAD,TierAR,TierAR.adj,No.Var.HM,No.Var.H,No.Var.M,No.Var.L,AD.matched,AR.matched,recessive,dominant,simple_repeat'
    feature_names = feature_names.split(",")
    nFeatures = len(feature_names)
    bivar_feature_mats = []

    gene = param['gene']
    varIDs, labeling = param['varIDs'], param['labeling']

    
    feature_df, out_folder = param['feature_df'], param['out_folder']

    seen_pairs = set() 
    
    allVars = feature_df.index.tolist()
    varIDs = [v for v in varIDs if v in allVars]

    if len(varIDs) == 0:
        return

    if len(varIDs) > 6:
        defaultPred = param['default_pred'].copy()
        defaultPred = defaultPred.loc[varIDs,:].sort_values("predict", ascending=False)
        varIDs = defaultPred.index.tolist()[:6]

    gene_feats = feature_df.loc[varIDs, feature_names].copy()
    gene_feats = gene_feats.values
    assert gene_feats.shape[0] == len(varIDs)
    # print(gene_feats, varIDs)

    def generate_feature_vector(i, j):
        fmat = np.zeros((2, nFeatures * 2 + 1))
        dist = np.absolute(int(varIDs[i].split("-")[1])-int(varIDs[j].split("-")[1]))
        
        fmat[0] = np.append(np.append(gene_feats[i], gene_feats[j]), [dist])
        fmat[1] = np.append(np.append(gene_feats[j], gene_feats[i]), [dist])

        f1_features = [f'{f}_1' for f in feature_names]
        f2_features = [f'{f}_2' for f in feature_names]
        
        return pd.DataFrame(data=fmat, columns=f1_features + f2_features + ['var_dist'],
                            index=[f'{varIDs[i]}_{varIDs[j]}', f'{varIDs[j]}_{varIDs[i]}'])
        

    for i in range(len(varIDs)):
        varii_feature_matrix = generate_feature_vector(i,i)            
        var1 = varIDs[i]
        
        if f'{var1}_{var1}' not in seen_pairs:               
            if feature_df.loc[var1, 'zyg'] == 2:

                if labeling:                        
                    if feature_df.loc[var1, 'is_strong'] == 1:
                        is_causal_wo_OMIM = 1
                    else:
                        is_causal_wo_OMIM = 0
                    varii_feature_matrix['is_causal'] = is_causal_wo_OMIM
                
                bivar_feature_mats.append(varii_feature_matrix)
                seen_pairs.add(f'{var1}_{var1}')
        
        if i < len(varIDs) - 1:
            for j in range(i+1, len(varIDs)):
                var2 = varIDs[j]
                if (f'{var1}_{var2}' not in seen_pairs) and (f'{var2}_{var1}' not in seen_pairs) and (var1 != var2):
                    varij_feature_matrix = generate_feature_vector(i, j)
                    
                    if labeling:
                        is_causal_wo_OMIM = 0

                        if (feature_df.loc[[var1, var2], 'is_strong'] == 1).all():
                            if (feature_df.loc[[var1, var2], 'zyg'] == 1).all():
                                is_causal_wo_OMIM = 1    

                        varij_feature_matrix['is_causal'] = is_causal_wo_OMIM

                    bivar_feature_mats.append(varij_feature_matrix)

                    seen_pairs.add(f'{var1}_{var2}')
                    seen_pairs.add(f'{var2}_{var1}')
                else:
                    pass

    if len(bivar_feature_mats) > 0:          
        bivar_feature_mats = pd.concat(bivar_feature_mats, axis=0)
        bivar_feature_mats = bivar_feature_mats.drop_duplicates()
        bivar_feature_mats.to_csv(f"{out_folder}/{gene}.csv")    
    
    return