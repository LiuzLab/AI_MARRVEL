from random import sample
import numpy as np
import pandas as pd

__author__ = "Chaozhong Liu"
__email__ = "chaozhol@bcm.edu"

###feature engineering and incorporation of tier features 
###keeps NA fields

def feature_engineering(score_file, tier_file):
    variable_name = ['varId_dash', 'hgmdSymptomScore',
                 'omimSymMatchFlag', 'hgmdSymMatchFlag', 'clinVarSymMatchFlag', 
                 'omimGeneFound', 'omimVarFound', 'hgmdGeneFound', 'hgmdVarFound', 
                 'clinVarVarFound', 'clinVarGeneFound', 
                 'clinvarTotalNumVars','clinvarNumP', 'clinvarNumLP', 'clinvarNumLB', 'clinvarNumB',
                 'dgvVarFound','decipherVarFound', 
                 'curationScoreHGMD', 'curationScoreOMIM', 'curationScoreClinVar', 
                 'conservationScoreDGV', 'omimSymptomSimScore', 'hgmdSymptomSimScore', 'clin_code', 
                 'GERPpp_RS', 'gnomadAF', 'gnomadAFg', 'LRT_score', 'LRT_Omega', 'phyloP100way_vertebrate', 
                 'gnomadGeneZscore', 'gnomadGenePLI', 'gnomadGeneOELof', 'gnomadGeneOELofUpper', 'IMPACT', 'CADD_phred',
                 'CADD_PHRED', 'DANN_score', 'REVEL_score', 'fathmm_MKL_coding_score', 'conservationScoreGnomad', 
                 'conservationScoreOELof', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 'SIFT_score', 'zyg',
                 'FATHMM_score', 'M_CAP_score', 'MutationAssessor_score', 'ESP6500_AA_AF',
                 'ESP6500_EA_AF', 'hom', 'hgmd_rs', 'spliceAImax', 'clin_PLP_perc', 'Consequence',
                 'nc_ClinVar_Exp', 'c_ClinVar_Exp', 'c_HGMD_Exp', 'nc_HGMD_Exp',
                 'nc_isPLP', 'nc_isBLB', 'c_isPLP', 'c_isBLB',
                 'nc_CLNREVSTAT', 'c_CLNREVSTAT',
                 'nc_RANKSCORE', 'c_RANKSCORE', 'CLASS', "phrank"]


    
    #sample_name = score_file.split('/')[-1].split('_')[1]
    #patient = pd.read_csv(score_file,sep='\t')
    patient = score_file

    patient = patient[variable_name]
    patient = patient.fillna('-')
    indel_index = [len(i.split('-')[-1]) != 1 or len(i.split('-')[-2]) != 1 for i in patient['varId_dash'].to_list()]

    patient.loc[patient['phrank'] == '-','phrank'] = np.NaN
    patient['phrank'] = patient['phrank'].astype('float64')

    patient.loc[patient['hgmdSymptomScore'] == '-','hgmdSymptomScore'] = np.NaN
    patient['hgmdSymptomScore'] = patient['hgmdSymptomScore'].astype('float64')

    patient.loc[patient['omimSymMatchFlag'] == '-','omimSymMatchFlag'] = np.NaN
    patient['omimSymMatchFlag'] = patient['omimSymMatchFlag'].astype('float64')

    patient.loc[patient['hgmdSymMatchFlag'] == '-','hgmdSymMatchFlag'] = np.NaN
    patient['hgmdSymMatchFlag'] = patient['hgmdSymMatchFlag'].astype('float64')

    patient.loc[patient['clinVarSymMatchFlag'] == '-','clinVarSymMatchFlag'] = np.NaN
    patient['clinVarSymMatchFlag'] = patient['clinVarSymMatchFlag'].astype('float64')


    patient['clinvarNumLP'] = (patient['clinvarNumLP']+patient['clinvarNumP']) / patient['clinvarTotalNumVars']
    patient['clinvarNumP'] = patient['clinvarNumP'] / patient['clinvarTotalNumVars']
    patient['clinvarNumLB'] = (patient['clinvarNumLB']+patient['clinvarNumB']) / patient['clinvarTotalNumVars']
    patient['clinvarNumB'] = patient['clinvarNumB'] / patient['clinvarTotalNumVars']
    patient = patient.drop(columns=['clinvarTotalNumVars'])
    variable_name.remove('clinvarTotalNumVars')


    patient.loc[patient['curationScoreHGMD'] == 'Low','curationScoreHGMD'] = 1
    patient.loc[patient['curationScoreHGMD'] == 'Medium','curationScoreHGMD'] = 2
    patient.loc[patient['curationScoreHGMD'] == 'High','curationScoreHGMD'] = 3
    patient['curationScoreHGMD'] = patient['curationScoreHGMD'].astype('float64')

    patient.loc[patient['curationScoreOMIM'] == 'Low','curationScoreOMIM'] = 1
    patient.loc[patient['curationScoreOMIM'] == 'Medium','curationScoreOMIM'] = 2
    patient.loc[patient['curationScoreOMIM'] == 'High','curationScoreOMIM'] = 3
    patient['curationScoreOMIM'] = patient['curationScoreOMIM'].astype('float64')

    patient.loc[patient['curationScoreClinVar'] == 'Low','curationScoreClinVar'] = 1
    patient.loc[patient['curationScoreClinVar'] == 'Medium','curationScoreClinVar'] = 2
    patient.loc[patient['curationScoreClinVar'] == 'High','curationScoreClinVar'] = 3
    patient['curationScoreClinVar'] = patient['curationScoreClinVar'].astype('float64')

    patient.loc[patient['conservationScoreDGV'] == 'Low','conservationScoreDGV'] = 1
    patient.loc[patient['conservationScoreDGV'] == 'Medium','conservationScoreDGV'] = 2
    patient.loc[patient['conservationScoreDGV'] == 'High','conservationScoreDGV'] = 3
    patient['conservationScoreDGV'] = patient['conservationScoreDGV'].astype('float64')


    patient.loc[patient['omimSymptomSimScore'] == '-','omimSymptomSimScore'] = np.NaN
    patient['omimSymptomSimScore'] = patient['omimSymptomSimScore'].astype('float64')
    
    patient.loc[patient['hgmdSymptomSimScore'] == '-','hgmdSymptomSimScore'] = np.NaN
    patient['hgmdSymptomSimScore'] = patient['hgmdSymptomSimScore'].astype('float64')


    patient['isB/LB'] = 0
    patient['isP/LP'] = 0
    patient.loc[patient['clin_code'].str.contains('Benign'),'isB/LB'] = 1
    patient.loc[patient['clin_code'].str.contains('Likely_benign'),'isB/LB'] = 1
    patient.loc[patient['clin_code'].str.contains('Conflicting_interpretations_of_pathogenicity'),'isB/LB'] = 0
    
    patient.loc[patient['clin_code'].str.contains('Likely_pathogenic'),'isP/LP'] = 1
    patient.loc[patient['clin_code'].str.contains('Pathogenic'),'isP/LP'] = 1
    patient.loc[patient['clin_code'].str.contains('Conflicting_interpretations_of_pathogenicity'),'isP/LP'] = 0

    patient.loc[patient['clin_PLP_perc']!='-','isP/LP'] = patient.loc[patient['clin_PLP_perc']!='-','clin_PLP_perc'].astype('float64')

    variable_name.append('isB/LB')
    variable_name.append('isP/LP')
    variable_name.remove('clin_PLP_perc')
    variable_name.remove('clin_code')


    patient.loc[patient['GERPpp_RS'] == '-', 'GERPpp_RS'] = np.NaN
    patient['GERPpp_RS'] = patient['GERPpp_RS'].astype('float64')

    
    patient.loc[patient['gnomadAFg'] == '-', 'gnomadAFg'] = patient.loc[patient['gnomadAFg'] == '-', 'gnomadAF'] #0.0 #
    patient.loc[patient['gnomadAFg'] == '.', 'gnomadAFg'] = patient.loc[patient['gnomadAFg'] == '.', 'gnomadAF'] #
    patient.loc[patient['gnomadAF'] == '-', 'gnomadAF'] = patient.loc[patient['gnomadAF'] == '-', 'gnomadAFg'] #0.0 #
    patient.loc[patient['gnomadAF'] == '.', 'gnomadAF'] = patient.loc[patient['gnomadAF'] == '.', 'gnomadAFg']
    patient.loc[patient['gnomadAFg'] == '-', 'gnomadAFg'] = np.NaN
    patient.loc[patient['gnomadAFg'] == '.', 'gnomadAFg'] = np.NaN
    patient.loc[patient['gnomadAF'] == '-', 'gnomadAF'] = np.NaN
    patient.loc[patient['gnomadAF'] == '.', 'gnomadAF'] = np.NaN
    patient['gnomadAFg'] = patient['gnomadAFg'].astype('float64')
    patient['gnomadAF'] = patient['gnomadAF'].astype('float64')
    #patient.loc[patient['gnomadAFg'] == -100, 'gnomadAFg'] = patient['gnomadAFg'].describe()['max']


    patient.loc[patient['LRT_score'] == '-', 'LRT_score'] = np.NaN
    patient['LRT_score'] = patient['LRT_score'].astype('float64')
    
    
    patient.loc[patient['LRT_Omega'] == '-', 'LRT_Omega'] = np.NaN
    patient['LRT_Omega'] = patient['LRT_Omega'].astype('float64')
    
    
    patient.loc[patient['phyloP100way_vertebrate'] == '-', 'phyloP100way_vertebrate'] = np.NaN
    patient['phyloP100way_vertebrate'] = patient['phyloP100way_vertebrate'].astype('float64')


    patient.loc[patient['gnomadGeneZscore'] == '-', 'gnomadGeneZscore'] = np.NaN
    patient['gnomadGeneZscore'] = patient['gnomadGeneZscore'].astype('float64')
    
    patient.loc[patient['gnomadGenePLI'] == '-', 'gnomadGenePLI'] = np.NaN
    patient['gnomadGenePLI'] = patient['gnomadGenePLI'].astype('float64')
    
    patient.loc[patient['gnomadGeneOELof'] == '-', 'gnomadGeneOELof'] = np.NaN
    patient['gnomadGeneOELof'] = patient['gnomadGeneOELof'].astype('float64')
    
    patient.loc[patient['gnomadGeneOELofUpper'] == '-', 'gnomadGeneOELofUpper'] = np.NaN
    patient['gnomadGeneOELofUpper'] = patient['gnomadGeneOELofUpper'].astype('float64')


    patient.loc[patient['IMPACT'] == '-','IMPACT'] = np.NaN
    patient.loc[patient['IMPACT'] == 'MODIFIER','IMPACT'] = 1
    patient.loc[patient['IMPACT'] == 'LOW','IMPACT'] = 2
    patient.loc[patient['IMPACT'] == 'MODERATE','IMPACT'] = 3
    patient.loc[patient['IMPACT'] == 'HIGH','IMPACT'] = 4
    patient['IMPACT'] = patient['IMPACT'].astype('float64')


    #"CADD_phred"
    patient.loc[patient['CADD_phred'] == '-', 'CADD_phred'] = np.NaN
    patient['CADD_phred'] = patient['CADD_phred'].astype('float64')

    patient.loc[patient['CADD_PHRED'] == '-', 'CADD_PHRED'] = np.NaN
    patient['CADD_PHRED'] = patient['CADD_PHRED'].astype('float64')


    #DANN_score
    patient.loc[patient['DANN_score'] == '-', 'DANN_score'] = np.NaN
    patient['DANN_score'] = patient['DANN_score'].astype('float64')


    #REVEL_score
    patient.loc[patient['REVEL_score'] == '-', 'REVEL_score'] = np.NaN
    patient['REVEL_score'] = patient['REVEL_score'].astype('float64')

    
    #fathmm_MKL_coding_score
    patient.loc[patient['fathmm_MKL_coding_score'] == '-', 'fathmm_MKL_coding_score'] = np.NaN
    patient['fathmm_MKL_coding_score'] = patient['fathmm_MKL_coding_score'].astype('float64')


    #conservationScoreGnomad
    patient.loc[patient['conservationScoreGnomad'] == '-','conservationScoreGnomad'] = np.NaN
    patient.loc[patient['conservationScoreGnomad'] == 'Low','conservationScoreGnomad'] = 1
    patient.loc[patient['conservationScoreGnomad'] == 'High','conservationScoreGnomad'] = 2
    patient['conservationScoreGnomad'] = patient['conservationScoreGnomad'].astype('float64')
    
    #conservationScoreOELof
    patient.loc[patient['conservationScoreOELof'] == '-','conservationScoreOELof'] = np.NaN
    patient.loc[patient['conservationScoreOELof'] == 'Low','conservationScoreOELof'] = 1
    patient.loc[patient['conservationScoreOELof'] == 'High','conservationScoreOELof'] = 2
    patient['conservationScoreOELof'] = patient['conservationScoreOELof'].astype('float64')


    #Polyphen2_HDIV_score
    patient['Polyphen2_HDIV_score'] = patient['Polyphen2_HDIV_score'] #.str.split(',')
    patient.loc[patient['Polyphen2_HDIV_score'] == '-','Polyphen2_HDIV_score'] = -100.0
    for i in patient[patient['Polyphen2_HDIV_score'] != -100.0].index:
        score_list = patient.loc[i, 'Polyphen2_HDIV_score'].split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'Polyphen2_HDIV_score'] = np.NaN
        else:
            patient.loc[i, 'Polyphen2_HDIV_score'] = max(score_list)
    patient['Polyphen2_HDIV_score'] = patient['Polyphen2_HDIV_score'].astype('float64')
    patient.loc[patient['Polyphen2_HVAR_score'] == -100, 'Polyphen2_HVAR_score'] = np.NaN
    
    #Polyphen2_HVAR_score
    patient['Polyphen2_HVAR_score'] = patient['Polyphen2_HVAR_score'] #.str.split(',')
    patient.loc[patient['Polyphen2_HVAR_score'] == '-','Polyphen2_HVAR_score'] = -100.0
    for i in patient[patient['Polyphen2_HVAR_score'] != -100.0].index:
        score_list = patient.loc[i, 'Polyphen2_HVAR_score'].split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'Polyphen2_HVAR_score'] = np.NaN
        else:
            patient.loc[i, 'Polyphen2_HVAR_score'] = max(score_list)
    patient['Polyphen2_HVAR_score'] = patient['Polyphen2_HVAR_score'].astype('float64')
    patient.loc[patient['Polyphen2_HVAR_score'] == -100, 'Polyphen2_HVAR_score'] = np.NaN


    #SIFT_score
    patient['SIFT_score'] = patient['SIFT_score'] #.str.split(',')
    patient.loc[patient['SIFT_score'] == '-','SIFT_score'] = -100.0
    for i in patient[patient['SIFT_score'] != -100.0].index:
        score_list = patient.loc[i, 'SIFT_score'].split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'SIFT_score'] = np.NaN
        else:
            patient.loc[i, 'SIFT_score'] = min(score_list)
    patient['SIFT_score'] = patient['SIFT_score'].astype('float64')
    patient.loc[patient['SIFT_score'] == -100, 'SIFT_score'] = np.NaN


    patient.loc[patient['zyg'] == 'HET','zyg'] = 1
    patient.loc[patient['zyg'] == 'HOM','zyg'] = 2
    patient.loc[patient['zyg'] == '-','zyg'] = np.NaN

    patient['zyg'] = patient['zyg'].astype('float64')


    #patient.loc[patient['GERPpp_NR'] == '-', 'GERPpp_NR'] = np.NaN
    #patient['GERPpp_NR'] = patient['GERPpp_NR'].astype('float64')


    #FATHMM_score
    patient['FATHMM_score'] = patient['FATHMM_score'] #.str.split(',')
    patient.loc[patient['FATHMM_score'] == '-','FATHMM_score'] = -100.0
    for i in patient[patient['FATHMM_score'] != -100.0].index:
        score_list = patient.loc[i, 'FATHMM_score'].split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'FATHMM_score'] = np.NaN
        else:
            patient.loc[i, 'FATHMM_score'] = min(score_list)
    patient['FATHMM_score'] = patient['FATHMM_score'].astype('float64')
    patient.loc[patient['FATHMM_score'] == -100, 'FATHMM_score'] = np.NaN


    #M_CAP_score
    patient.loc[patient['M_CAP_score'] == '-', 'M_CAP_score'] = np.NaN
    patient['M_CAP_score'] = patient['M_CAP_score'].astype('float64')


    #MutationAssessor_score
    patient['MutationAssessor_score'] = patient['MutationAssessor_score'] #.str.split(',')
    patient.loc[patient['MutationAssessor_score'] == '-','MutationAssessor_score'] = -100.0
    for i in patient[patient['MutationAssessor_score'] != -100.0].index:
        score_list = patient.loc[i, 'MutationAssessor_score'].split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'MutationAssessor_score'] = np.NaN
        else:
            patient.loc[i, 'MutationAssessor_score'] = max(score_list)
    patient['MutationAssessor_score'] = patient['MutationAssessor_score'].astype('float64')
    patient.loc[patient['MutationAssessor_score'] == -100, 'MutationAssessor_score'] = np.NaN



    #ESP6500_AA_AF
    patient.loc[patient['ESP6500_AA_AF'] == '-', 'ESP6500_AA_AF'] = np.NaN
    patient['ESP6500_AA_AF'] = patient['ESP6500_AA_AF'].astype('float64')


    #ESP6500_AA_AF
    patient.loc[patient['ESP6500_EA_AF'] == '-', 'ESP6500_EA_AF'] = np.NaN
    patient['ESP6500_EA_AF'] = patient['ESP6500_EA_AF'].astype('float64')


    #hom
    patient.loc[patient['hom'] == '-', 'hom'] = np.NaN
    patient['hom'] = patient['hom'].astype('float64')

    patient['hgmd_rs'] = patient['hgmd_rs'].apply(lambda x: x.split(",")[0])
    patient.loc[patient['hgmd_rs'] == '-','hgmd_rs'] = np.NaN
    patient['hgmd_rs'] = patient['hgmd_rs'].astype('float64')


    #spliceAImax
    patient.loc[patient['spliceAImax'] == '-', 'spliceAImax'] = np.NaN
    patient['spliceAImax'] = patient['spliceAImax'].astype('float64')
    patient.loc[pd.isna(patient['spliceAImax']), 'spliceAImax'] = patient.loc[~pd.isna(patient['spliceAImax']), 'spliceAImax'].describe()['min']


    #Consequence
    consequence = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 
               'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification', 'inframe_insertion', 
               'inframe_deletion', 'missense_variant', 'protein_altering_variant', 'splice_region_variant',
               'splice_donor_5th_base_variant', 'splice_donor_region_variant']

    for cons in consequence:
        patient['cons_%s'%cons] = patient['Consequence'].str.contains(cons).astype('int')
        variable_name.append('cons_%s'%cons)

    variable_name.remove('Consequence')


    #nc_isPLP
    patient.loc[patient['nc_isPLP']!=True, 'nc_isPLP'] = 0
    patient.loc[patient['nc_isPLP']==True, 'nc_isPLP'] = 1
    patient['nc_isPLP'] = patient['nc_isPLP'].astype('float64')

    #nc_isBLB
    patient.loc[patient['nc_isBLB']!=True, 'nc_isBLB'] = 0
    patient.loc[patient['nc_isBLB']==True, 'nc_isBLB'] = 1
    patient['nc_isBLB'] = patient['nc_isBLB'].astype('float64')

    #c_isPLP
    patient.loc[patient['c_isPLP']!=True, 'c_isPLP'] = 0
    patient.loc[patient['c_isPLP']==True, 'c_isPLP'] = 1
    patient['c_isPLP'] = patient['c_isPLP'].astype('float64')

    #c_isBLB
    patient.loc[patient['c_isBLB']!=True, 'c_isBLB'] = 0
    patient.loc[patient['c_isBLB']==True, 'c_isBLB'] = 1
    patient['c_isBLB'] = patient['c_isBLB'].astype('float64')


    #nc_ClinVar_Exp
    patient.loc[patient['nc_ClinVar_Exp']!='nonCoding','nc_ClinVar_Exp'] = 0
    patient.loc[patient['nc_ClinVar_Exp']=='nonCoding','nc_ClinVar_Exp'] = 1
    patient['nc_ClinVar_Exp'] = patient['nc_ClinVar_Exp'].astype('float64')

    #nc_HGMD_Exp
    patient.loc[patient['nc_HGMD_Exp']!='nonCoding','nc_HGMD_Exp'] = 0
    patient.loc[patient['nc_HGMD_Exp']=='nonCoding','nc_HGMD_Exp'] = 1
    patient['nc_HGMD_Exp'] = patient['nc_HGMD_Exp'].astype('float64')

    #c_ClinVar_Exp
    c_clin_exp = ['Del_to_Missense', 'Different_pChange', 'Same_pChange']

    for exp in c_clin_exp:

        patient['c_ClinVar_Exp_%s'%exp] = patient['c_ClinVar_Exp'].str.contains(exp).astype('int')
        variable_name.append('c_ClinVar_Exp_%s'%exp)

    variable_name.remove('c_ClinVar_Exp')

    #c_HGMD_Exp
    c_hgmd_exp = ['Del_to_Missense', 'Different_pChange', 'Same_pChange', 'Stop_Loss', 'Start_Loss']

    for exp in c_hgmd_exp:

        patient['c_HGMD_Exp_%s'%exp] = patient['c_HGMD_Exp'].str.contains(exp).astype('int')
        variable_name.append('c_HGMD_Exp_%s'%exp)

    variable_name.remove('c_HGMD_Exp')


    #nc_CLNREVSTAT
    CLN_STAT = ['no_assertion_provided', 'no_assertion_criteria_provided','no_assertion_for_the_individual_variant',
                'criteria_provided,_single_submitter', 'criteria_provided,_conflicting_interpretations',
                'criteria_provided,_multiple_submitters,_no_conflicts', 'reviewed_by_expert_panel', 'practice_guideline']

    CLN_STAT_SCORE = [0,0,0,1,1,2,3,4]
    patient.loc[patient['nc_CLNREVSTAT'] == '-', 'nc_CLNREVSTAT'] = np.NaN
    patient.loc[patient['c_CLNREVSTAT'] == '-', 'c_CLNREVSTAT'] = np.NaN

    for i in range(len(CLN_STAT)):
        patient.loc[patient['nc_CLNREVSTAT'] == CLN_STAT[i], 'nc_CLNREVSTAT'] = CLN_STAT_SCORE[i]
        patient.loc[patient['c_CLNREVSTAT'] == CLN_STAT[i], 'c_CLNREVSTAT'] = CLN_STAT_SCORE[i]

    patient['nc_CLNREVSTAT'] = patient['nc_CLNREVSTAT'].astype('float64')
    patient['c_CLNREVSTAT'] = patient['c_CLNREVSTAT'].astype('float64')


    #nc_RANKSCORE
    patient.loc[patient['nc_RANKSCORE'] == '-', 'nc_RANKSCORE'] = np.NaN
    patient['nc_RANKSCORE'] = patient['nc_RANKSCORE'].astype('float64')

    #nc_RANKSCORE
    patient.loc[patient['c_RANKSCORE'] == '-', 'c_RANKSCORE'] = np.NaN
    patient['c_RANKSCORE'] = patient['c_RANKSCORE'].astype('float64')

    #CLASS
    patient.loc[patient['CLASS'] == '-', 'CLASS'] = np.NaN
    patient.loc[patient['CLASS'] == 'DM?', 'CLASS'] = 1
    patient.loc[patient['CLASS'] == 'DM', 'CLASS'] = 2
    patient['CLASS'] = patient['CLASS'].astype('float64')


    patient.loc[:,['gnomadAF','gnomadAFg','gnomadGeneOELof','gnomadGeneOELofUpper',
    'SIFT_score', 'FATHMM_score', 'ESP6500_AA_AF', 'ESP6500_EA_AF']] = -patient.loc[:,['gnomadAF','gnomadAFg','gnomadGeneOELof','gnomadGeneOELofUpper','SIFT_score', 'FATHMM_score', 'ESP6500_AA_AF', 'ESP6500_EA_AF']]
    patient = patient.groupby(['varId_dash'], sort=False)[variable_name[1:]].max()
    patient.loc[:,['gnomadAF','gnomadAFg','gnomadGeneOELof','gnomadGeneOELofUpper',
    'SIFT_score', 'FATHMM_score', 'ESP6500_AA_AF', 'ESP6500_EA_AF']] = -patient.loc[:,['gnomadAF','gnomadAFg','gnomadGeneOELof','gnomadGeneOELofUpper','SIFT_score', 'FATHMM_score', 'ESP6500_AA_AF', 'ESP6500_EA_AF']]


    #adding in tier features
    #tier = pd.read_csv(tier_file, sep='\t')
    tier=tier_file
    tier_vars = ['IMPACT.from.Tier','TierAD','TierAR','TierAR.adj',
                 'No.Var.HM','No.Var.H','No.Var.M','No.Var.L',
                 'AD.matched','AR.matched',
                 'OMIM.Inheritance.recessive', 'OMIM.Inheritance.dominant',
                 'OMIM.Inheritance.both', 'OMIM.Inheritance.others']

    tier['OMIM.Inheritance.recessive'] = 0
    tier['OMIM.Inheritance.dominant'] = 0
    tier['OMIM.Inheritance.both'] = 0
    tier['OMIM.Inheritance.others'] = 0
    tier.loc[tier['OMIM.Inheritance'] == 'recessive', 'OMIM.Inheritance.recessive'] = 1
    tier.loc[tier['OMIM.Inheritance'] == 'dominant', 'OMIM.Inheritance.dominant'] = 1
    tier.loc[tier['OMIM.Inheritance'] == 'both', 'OMIM.Inheritance.both'] = 1
    tier.loc[tier['OMIM.Inheritance'] == 'others', 'OMIM.Inheritance.others'] = 1

    tier.loc[:,['TierAD','TierAR','TierAR.adj']] = -tier.loc[:,['TierAD','TierAR','TierAR.adj']]
    tier = tier.groupby(['Uploaded_variation'], sort=False)[tier_vars].max()
    tier.loc[:,['TierAD','TierAR','TierAR.adj']] = -tier.loc[:,['TierAD','TierAR','TierAR.adj']]
    
    patient = pd.concat([patient, tier], axis=1)


    return patient




    

