from random import sample
import numpy as np
import pandas as pd

__author__ = "Chaozhong Liu"
__email__ = "chaozhol@bcm.edu"

###Feature engineering and incorporation of tier module outputs. 
###fills NA w/ 0

def getValFromStr(valStr: str, select: str = 'min'):
    select_method = {'min':min, 'max':max}
    vals = valStr.split(',')
    if '-' in vals or '.' in vals:
        return np.NaN
    else:
        vals = [float(i) for i in vals]
        return select_method[select](vals)

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
    feature_stat = pd.read_csv('/run/data_dependencies/annotate/feature_stats.csv', index_col=0)

    patient = score_file.copy()

    patient = patient[variable_name]
    patient = patient.fillna('-')
    indel_index = [len(i.split('-')[-1]) != 1 or len(i.split('-')[-2]) != 1 for i in patient['varId_dash'].to_list()]

    patient.loc[patient['phrank'] == '-','phrank'] = 0
    patient['phrank'] = patient['phrank'].astype('float64')

    patient.loc[patient['hgmdSymptomScore'] == '-','hgmdSymptomScore'] = 0
    patient['hgmdSymptomScore'] = patient['hgmdSymptomScore'].astype('float64')

    patient.loc[patient['omimSymMatchFlag'] == '-','omimSymMatchFlag'] = 0
    patient['omimSymMatchFlag'] = patient['omimSymMatchFlag'].astype('float64')

    patient.loc[patient['hgmdSymMatchFlag'] == '-','hgmdSymMatchFlag'] = 0
    patient['hgmdSymMatchFlag'] = patient['hgmdSymMatchFlag'].astype('float64')

    patient.loc[patient['clinVarSymMatchFlag'] == '-','clinVarSymMatchFlag'] = 0
    patient['clinVarSymMatchFlag'] = patient['clinVarSymMatchFlag'].astype('float64')


    patient['clinvarNumLP'] = (patient['clinvarNumLP']+patient['clinvarNumP']) / patient['clinvarTotalNumVars']
    patient['clinvarNumLP'] = patient['clinvarNumLP'].fillna(0.0)
    patient['clinvarNumP'] = patient['clinvarNumP'] / patient['clinvarTotalNumVars']
    patient['clinvarNumP'] = patient['clinvarNumP'].fillna(0.0)
    patient['clinvarNumLB'] = (patient['clinvarNumLB']+patient['clinvarNumB']) / patient['clinvarTotalNumVars']
    patient['clinvarNumLB'] = patient['clinvarNumLB'].fillna(1.0)
    patient['clinvarNumB'] = patient['clinvarNumB'] / patient['clinvarTotalNumVars']
    patient['clinvarNumB'] = patient['clinvarNumB'].fillna(1.0)
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


    patient.loc[patient['omimSymptomSimScore'] == '-','omimSymptomSimScore'] = 0.0
    patient['omimSymptomSimScore'] = patient['omimSymptomSimScore'].astype('float64')
    
    patient.loc[patient['hgmdSymptomSimScore'] == '-','hgmdSymptomSimScore'] = 0.0
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
    if np.all(pd.isna(patient['GERPpp_RS'])):
        # Single variant query
        patient.loc[(pd.isna(patient['GERPpp_RS'])) & np.array(indel_index), 'GERPpp_RS'] = feature_stat.loc['GERPpp_RS', 'mean']
        patient.loc[pd.isna(patient['GERPpp_RS']), 'GERPpp_RS'] = feature_stat.loc['GERPpp_RS', 'mean']
        patient['GERPpp_RS'] = patient['GERPpp_RS'].astype('float64')
    else:
        # Patient whole VCF
        patient['GERPpp_RS'] = patient['GERPpp_RS'].astype('float64')
        patient.loc[(pd.isna(patient['GERPpp_RS'])) & np.array(indel_index), 'GERPpp_RS'] = patient['GERPpp_RS'].describe()['mean']
        patient.loc[pd.isna(patient['GERPpp_RS']), 'GERPpp_RS'] = patient['GERPpp_RS'].describe()['mean']


    gnomadAFVal = patient['gnomadAF'].copy()
    gnomadAFVal = np.array([ getValFromStr(str(i), 'min') for i in gnomadAFVal ])
    gnomadAFVal = gnomadAFVal.astype(float)
    patient['gnomadAF'] = gnomadAFVal

    gnomadAFgVal = patient['gnomadAFg'].copy()
    gnomadAFgVal = np.array([ getValFromStr(str(i), 'min') for i in gnomadAFgVal ])
    gnomadAFgVal = gnomadAFgVal.astype(float)
    patient['gnomadAFg'] = gnomadAFgVal
    patient.loc[patient['gnomadAFg'].isna(), 'gnomadAFg'] = patient.loc[patient['gnomadAFg'].isna(), 'gnomadAF']
    patient.loc[patient['gnomadAF'].isna(), 'gnomadAF'] = patient.loc[patient['gnomadAF'].isna(), 'gnomadAFg']
    patient['gnomadAF'] = patient['gnomadAF'].fillna(0.0)
    patient['gnomadAFg'] = patient['gnomadAFg'].fillna(0.0)

    #patient.loc[patient['gnomadAFg'] == '-', 'gnomadAFg'] = patient.loc[patient['gnomadAFg'] == '-', 'gnomadAF'] #0.0 #
    #patient.loc[patient['gnomadAFg'] == '.', 'gnomadAFg'] = patient.loc[patient['gnomadAFg'] == '.', 'gnomadAF'] #
    #patient.loc[patient['gnomadAF'] == '-', 'gnomadAF'] = patient.loc[patient['gnomadAF'] == '-', 'gnomadAFg'] #0.0 #
    #patient.loc[patient['gnomadAF'] == '.', 'gnomadAF'] = patient.loc[patient['gnomadAF'] == '.', 'gnomadAFg']
    #patient.loc[patient['gnomadAFg'] == '-', 'gnomadAFg'] = 0.0
    #patient.loc[patient['gnomadAFg'] == '.', 'gnomadAFg'] = 0.0
    #patient.loc[patient['gnomadAF'] == '-', 'gnomadAF'] = 0.0
    #patient.loc[patient['gnomadAF'] == '.', 'gnomadAF'] = 0.0
    #patient['gnomadAFg'] = patient['gnomadAFg'].astype('float64')
    #patient['gnomadAF'] = patient['gnomadAF'].astype('float64')
    #patient.loc[patient['gnomadAFg'] == -100, 'gnomadAFg'] = patient['gnomadAFg'].describe()['max']


    patient.loc[patient['LRT_score'] == '-', 'LRT_score'] = np.NaN
    if np.all(pd.isna(patient['LRT_score'])):
        # Single variant query
        patient.loc[(pd.isna(patient['LRT_score'])) & np.array(indel_index), 'LRT_score'] = feature_stat.loc['LRT_score', 'mean']
        patient.loc[pd.isna(patient['LRT_score']), 'LRT_score'] = feature_stat.loc['LRT_score', 'mean']
        patient['LRT_score'] = patient['LRT_score'].astype('float64')
    else:
        # Patient whole VCF
        patient['LRT_score'] = patient['LRT_score'].astype('float64')
        patient.loc[(pd.isna(patient['LRT_score'])) & np.array(indel_index), 'LRT_score'] = patient['LRT_score'].describe()['mean']
        patient.loc[pd.isna(patient['LRT_score']), 'LRT_score'] = patient['LRT_score'].describe()['mean']
    
    
    patient.loc[patient['LRT_Omega'] == '-', 'LRT_Omega'] = np.NaN
    if np.all(pd.isna(patient['LRT_Omega'])):
        # Single variant query
        patient.loc[(pd.isna(patient['LRT_Omega'])) & np.array(indel_index), 'LRT_Omega'] = feature_stat.loc['LRT_Omega', 'mean']
        patient.loc[pd.isna(patient['LRT_Omega']), 'LRT_Omega'] = feature_stat.loc['LRT_Omega', 'mean']
        patient['LRT_Omega'] = patient['LRT_Omega'].astype('float64')
    else:
        # Patient whole VCF
        patient['LRT_Omega'] = patient['LRT_Omega'].astype('float64')
        patient.loc[(pd.isna(patient['LRT_Omega'])) & np.array(indel_index), 'LRT_Omega'] = patient['LRT_Omega'].describe()['mean']
        patient.loc[pd.isna(patient['LRT_Omega']), 'LRT_Omega'] = patient['LRT_Omega'].describe()['mean']
    
    
    patient.loc[patient['phyloP100way_vertebrate'] == '-', 'phyloP100way_vertebrate'] = np.NaN
    if np.all(pd.isna(patient['phyloP100way_vertebrate'])):
        # Single variant query
        patient.loc[(pd.isna(patient['phyloP100way_vertebrate'])) & np.array(indel_index), 'phyloP100way_vertebrate'] = feature_stat.loc['phyloP100way_vertebrate', 'mean']
        patient.loc[pd.isna(patient['phyloP100way_vertebrate']), 'phyloP100way_vertebrate'] = feature_stat.loc['phyloP100way_vertebrate', 'mean']
        patient['phyloP100way_vertebrate'] = patient['phyloP100way_vertebrate'].astype('float64')
    else:
        # Patient whole VCF
        patient['phyloP100way_vertebrate'] = patient['phyloP100way_vertebrate'].astype('float64')
        patient.loc[(pd.isna(patient['phyloP100way_vertebrate'])) & np.array(indel_index), 'phyloP100way_vertebrate'] = patient['phyloP100way_vertebrate'].describe()['mean']
        patient.loc[pd.isna(patient['phyloP100way_vertebrate']), 'phyloP100way_vertebrate'] = patient['phyloP100way_vertebrate'].describe()['mean']


    patient.loc[patient['gnomadGeneZscore'] == '-', 'gnomadGeneZscore'] = np.NaN
    if np.all(pd.isna(patient['gnomadGeneZscore'])):
        # Single variant query
        patient.loc[pd.isna(patient['gnomadGeneZscore']), 'gnomadGeneZscore'] = feature_stat.loc['gnomadGeneZscore', 'mean']
        patient['gnomadGeneZscore'] = patient['gnomadGeneZscore'].astype('float64')
    else:
        # Patient whole VCF
        patient['gnomadGeneZscore'] = patient['gnomadGeneZscore'].astype('float64')
        patient.loc[pd.isna(patient['gnomadGeneZscore']), 'gnomadGeneZscore'] = patient['gnomadGeneZscore'].describe()['mean']


    patient.loc[patient['gnomadGenePLI'] == '-', 'gnomadGenePLI'] = np.NaN
    if np.all(pd.isna(patient['gnomadGenePLI'])):
        # Single variant query
        patient.loc[pd.isna(patient['gnomadGenePLI']), 'gnomadGenePLI'] = feature_stat.loc['gnomadGenePLI', 'mean']
        patient['gnomadGenePLI'] = patient['gnomadGenePLI'].astype('float64')
    else:
        # Patient whole VCF
        patient['gnomadGenePLI'] = patient['gnomadGenePLI'].astype('float64')
        patient.loc[pd.isna(patient['gnomadGenePLI']), 'gnomadGenePLI'] = patient['gnomadGenePLI'].describe()['mean']

    
    patient.loc[patient['gnomadGeneOELof'] == '-', 'gnomadGeneOELof'] = np.NaN
    if np.all(pd.isna(patient['gnomadGeneOELof'])):
        # Single variant query
        patient.loc[pd.isna(patient['gnomadGeneOELof']), 'gnomadGeneOELof'] = feature_stat.loc['gnomadGeneOELof', 'mean']
        patient['gnomadGeneOELof'] = patient['gnomadGeneOELof'].astype('float64')
    else:
        # Patient whole VCF
        patient['gnomadGeneOELof'] = patient['gnomadGeneOELof'].astype('float64')
        patient.loc[pd.isna(patient['gnomadGeneOELof']), 'gnomadGeneOELof'] = patient['gnomadGeneOELof'].describe()['mean']
    
    patient.loc[patient['gnomadGeneOELofUpper'] == '-', 'gnomadGeneOELofUpper'] = np.NaN
    if np.all(pd.isna(patient['gnomadGeneOELofUpper'])):
        # Single variant query
        patient.loc[pd.isna(patient['gnomadGeneOELofUpper']), 'gnomadGeneOELofUpper'] = feature_stat.loc['gnomadGeneOELofUpper', 'mean']
        patient['gnomadGeneOELofUpper'] = patient['gnomadGeneOELofUpper'].astype('float64')
    else:
        # Patient whole VCF
        patient['gnomadGeneOELofUpper'] = patient['gnomadGeneOELofUpper'].astype('float64')
        patient.loc[pd.isna(patient['gnomadGeneOELofUpper']), 'gnomadGeneOELofUpper'] = patient['gnomadGeneOELofUpper'].describe()['mean']


    patient.loc[patient['IMPACT'] == '-','IMPACT'] = 0
    patient.loc[patient['IMPACT'] == 'MODIFIER','IMPACT'] = 1
    patient.loc[patient['IMPACT'] == 'LOW','IMPACT'] = 2
    patient.loc[patient['IMPACT'] == 'MODERATE','IMPACT'] = 3
    patient.loc[patient['IMPACT'] == 'HIGH','IMPACT'] = 4
    patient['IMPACT'] = patient['IMPACT'].astype('float64')


    #"CADD_phred"
    patient.loc[patient['CADD_phred'] == '-', 'CADD_phred'] = np.NaN
    if np.all(pd.isna(patient['CADD_phred'])):
        # Single variant query
        patient.loc[(pd.isna(patient['CADD_phred'])) & np.array(indel_index), 'CADD_phred'] = feature_stat.loc['CADD_phred', 'mean']
        patient.loc[pd.isna(patient['CADD_phred']), 'CADD_phred'] = feature_stat.loc['CADD_phred', 'mean']
        patient['CADD_phred'] = patient['CADD_phred'].astype('float64')
    else:
        # Patient whole VCF
        patient['CADD_phred'] = patient['CADD_phred'].astype('float64')
        patient.loc[(pd.isna(patient['CADD_phred'])) & np.array(indel_index), 'CADD_phred'] = patient['CADD_phred'].describe()['mean']
        patient.loc[pd.isna(patient['CADD_phred']), 'CADD_phred'] = patient['CADD_phred'].describe()['mean']


    patient.loc[patient['CADD_PHRED'] == '-', 'CADD_PHRED'] = np.NaN
    if np.all(pd.isna(patient['CADD_PHRED'])):
        # Single variant query
        patient.loc[(pd.isna(patient['CADD_PHRED'])) & np.array(indel_index), 'CADD_PHRED'] = feature_stat.loc['CADD_PHRED', 'mean']
        patient.loc[pd.isna(patient['CADD_PHRED']), 'CADD_PHRED'] = feature_stat.loc['CADD_PHRED', 'mean']
        patient['CADD_PHRED'] = patient['CADD_PHRED'].astype('float64')
    else:
        # Patient whole VCF
        patient['CADD_PHRED'] = patient['CADD_PHRED'].astype('float64')
        patient.loc[(pd.isna(patient['CADD_PHRED'])) & np.array(indel_index), 'CADD_PHRED'] = patient['CADD_PHRED'].describe()['mean']
        patient.loc[pd.isna(patient['CADD_PHRED']), 'CADD_PHRED'] = patient['CADD_PHRED'].describe()['mean']


    #DANN_score
    patient.loc[patient['DANN_score'] == '-', 'DANN_score'] = np.NaN
    if np.all(pd.isna(patient['DANN_score'])):
        # Single variant query
        patient.loc[(pd.isna(patient['DANN_score'])) & np.array(indel_index), 'DANN_score'] = feature_stat.loc['DANN_score', 'mean']
        patient.loc[pd.isna(patient['DANN_score']), 'DANN_score'] = feature_stat.loc['DANN_score', 'mean']
        patient['DANN_score'] = patient['DANN_score'].astype('float64')
    else:
        # Patient whole VCF
        patient['DANN_score'] = patient['DANN_score'].astype('float64')
        patient.loc[(pd.isna(patient['DANN_score'])) & np.array(indel_index), 'DANN_score'] = patient['DANN_score'].describe()['mean']
        patient.loc[pd.isna(patient['DANN_score']), 'DANN_score'] = patient['DANN_score'].describe()['mean']


    #REVEL_score
    patient.loc[patient['REVEL_score'] == '-', 'REVEL_score'] = np.NaN
    for i in patient[~patient['REVEL_score'].isna()].index:
        score_list = str(patient.loc[i, 'REVEL_score']).split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'REVEL_score'] = np.NaN
        else:
            patient.loc[i, 'REVEL_score'] = max(score_list)

    if np.all(pd.isna(patient['REVEL_score'])):
        # Single variant query
        patient.loc[(pd.isna(patient['REVEL_score'])) & np.array(indel_index), 'REVEL_score'] = feature_stat.loc['REVEL_score', 'mean']
        patient.loc[pd.isna(patient['REVEL_score']), 'REVEL_score'] = feature_stat.loc['REVEL_score', 'mean']
        patient['REVEL_score'] = patient['REVEL_score'].astype('float64')
    else:
        # Patient whole VCF
        patient['REVEL_score'] = patient['REVEL_score'].astype('float64')
        patient.loc[(patient['REVEL_score'].isna()) & np.array(indel_index), 'REVEL_score'] = patient['REVEL_score'].describe()['mean']
        patient.loc[pd.isna(patient['REVEL_score']), 'REVEL_score'] = patient['REVEL_score'].describe()['mean']


    
    #fathmm_MKL_coding_score
    patient.loc[patient['fathmm_MKL_coding_score'] == '-', 'fathmm_MKL_coding_score'] = np.NaN
    if np.all(pd.isna(patient['fathmm_MKL_coding_score'])):
        # Single variant query
        patient.loc[(pd.isna(patient['fathmm_MKL_coding_score'])) & np.array(indel_index), 'fathmm_MKL_coding_score'] = feature_stat.loc['fathmm_MKL_coding_score', 'mean']
        patient.loc[pd.isna(patient['fathmm_MKL_coding_score']), 'fathmm_MKL_coding_score'] = feature_stat.loc['fathmm_MKL_coding_score', 'mean']
        patient['fathmm_MKL_coding_score'] = patient['fathmm_MKL_coding_score'].astype('float64')
    else:
        # Patient whole VCF
        patient['fathmm_MKL_coding_score'] = patient['fathmm_MKL_coding_score'].astype('float64')
        patient.loc[(pd.isna(patient['fathmm_MKL_coding_score'])) & np.array(indel_index), 'fathmm_MKL_coding_score'] = patient['fathmm_MKL_coding_score'].describe()['mean']
        patient.loc[pd.isna(patient['fathmm_MKL_coding_score']), 'fathmm_MKL_coding_score'] = patient['fathmm_MKL_coding_score'].describe()['mean']
    

    #conservationScoreGnomad
    patient.loc[patient['conservationScoreGnomad'] == '-','conservationScoreGnomad'] = 1
    patient.loc[patient['conservationScoreGnomad'] == 'Low','conservationScoreGnomad'] = 1
    patient.loc[patient['conservationScoreGnomad'] == 'High','conservationScoreGnomad'] = 2
    patient['conservationScoreGnomad'] = patient['conservationScoreGnomad'].astype('float64')
    
    #conservationScoreOELof
    patient.loc[patient['conservationScoreOELof'] == '-','conservationScoreOELof'] = 1
    patient.loc[patient['conservationScoreOELof'] == 'Low','conservationScoreOELof'] = 1
    patient.loc[patient['conservationScoreOELof'] == 'High','conservationScoreOELof'] = 2
    patient['conservationScoreOELof'] = patient['conservationScoreOELof'].astype('float64')
    

    #Polyphen2_HDIV_score
    patient['Polyphen2_HDIV_score'] = patient['Polyphen2_HDIV_score']
    patient.loc[patient['Polyphen2_HDIV_score'] == '-','Polyphen2_HDIV_score'] = np.NaN
    for i in patient[~patient['Polyphen2_HDIV_score'].isna()].index:
        score_list = str(patient.loc[i, 'Polyphen2_HDIV_score']).split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'Polyphen2_HDIV_score'] = np.NaN
        else:
            patient.loc[i, 'Polyphen2_HDIV_score'] = max(score_list)

    if np.all(pd.isna(patient['Polyphen2_HDIV_score'])):
        # Single variant query
        patient.loc[(pd.isna(patient['Polyphen2_HDIV_score'])) & np.array(indel_index), 'Polyphen2_HDIV_score'] = feature_stat.loc['Polyphen2_HDIV_score', '50%']
        patient.loc[pd.isna(patient['Polyphen2_HDIV_score']), 'Polyphen2_HDIV_score'] = feature_stat.loc['Polyphen2_HDIV_score', '50%']
        patient['Polyphen2_HDIV_score'] = patient['Polyphen2_HDIV_score'].astype('float64')
    else:
        # Patient whole VCF
        patient['Polyphen2_HDIV_score'] = patient['Polyphen2_HDIV_score'].astype('float64')
        patient.loc[(pd.isna(patient['Polyphen2_HDIV_score'])) & np.array(indel_index), 'Polyphen2_HDIV_score'] = patient['Polyphen2_HDIV_score'].describe()['50%']
        patient.loc[pd.isna(patient['Polyphen2_HDIV_score']), 'Polyphen2_HDIV_score'] = patient['Polyphen2_HDIV_score'].describe()['50%']

    
    #Polyphen2_HVAR_score
    patient['Polyphen2_HVAR_score'] = patient['Polyphen2_HVAR_score']
    patient.loc[patient['Polyphen2_HVAR_score'] == '-','Polyphen2_HVAR_score'] = np.NaN
    for i in patient[~patient['Polyphen2_HVAR_score'].isna()].index:
        score_list = str(patient.loc[i, 'Polyphen2_HVAR_score']).split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'Polyphen2_HVAR_score'] = np.NaN
        else:
            patient.loc[i, 'Polyphen2_HVAR_score'] = max(score_list)

    if np.all(pd.isna(patient['Polyphen2_HVAR_score'])):
        # Single variant query
        patient.loc[(pd.isna(patient['Polyphen2_HVAR_score'])) & np.array(indel_index), 'Polyphen2_HVAR_score'] = feature_stat.loc['Polyphen2_HVAR_score', '50%']
        patient.loc[pd.isna(patient['Polyphen2_HVAR_score']), 'Polyphen2_HVAR_score'] = feature_stat.loc['Polyphen2_HVAR_score', '50%']
        patient['Polyphen2_HVAR_score'] = patient['Polyphen2_HVAR_score'].astype('float64')
    else:
        # Patient whole VCF
        patient['Polyphen2_HVAR_score'] = patient['Polyphen2_HVAR_score'].astype('float64')
        patient.loc[(pd.isna(patient['Polyphen2_HVAR_score'])) & np.array(indel_index), 'Polyphen2_HVAR_score'] = patient['Polyphen2_HVAR_score'].describe()['50%']
        patient.loc[pd.isna(patient['Polyphen2_HVAR_score']), 'Polyphen2_HVAR_score'] = patient['Polyphen2_HVAR_score'].describe()['50%']


    #SIFT_score
    patient['SIFT_score'] = patient['SIFT_score']
    patient.loc[patient['SIFT_score'] == '-','SIFT_score'] = np.NaN
    for i in patient[~patient['SIFT_score'].isna()].index:
        score_list = str(patient.loc[i, 'SIFT_score']).split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'SIFT_score'] = np.NaN
        else:
            patient.loc[i, 'SIFT_score'] = min(score_list)

    if np.all(pd.isna(patient['SIFT_score'])):
        # Single variant query
        patient.loc[(pd.isna(patient['SIFT_score'])) & np.array(indel_index), 'SIFT_score'] = feature_stat.loc['SIFT_score', '50%']
        patient.loc[pd.isna(patient['SIFT_score']), 'SIFT_score'] = feature_stat.loc['SIFT_score', '50%']
        patient['SIFT_score'] = patient['SIFT_score'].astype('float64')
    else:
        # Patient whole VCF
        patient['SIFT_score'] = patient['SIFT_score'].astype('float64')
        patient.loc[(pd.isna(patient['SIFT_score'])) & np.array(indel_index), 'SIFT_score'] = patient['SIFT_score'].describe()['50%']
        patient.loc[pd.isna(patient['SIFT_score']), 'SIFT_score'] = patient['SIFT_score'].describe()['50%']


    patient.loc[patient['zyg'] == 'HET','zyg'] = 1
    patient.loc[patient['zyg'] == 'HOM','zyg'] = 2
    patient.loc[patient['zyg'] == '-','zyg'] = 0

    patient['zyg'] = patient['zyg'].astype('float64')


    #patient.loc[patient['GERPpp_NR'] == '-', 'GERPpp_NR'] = np.NaN
    #patient['GERPpp_NR'] = patient['GERPpp_NR'].astype('float64')
    #patient.loc[(pd.isna(patient['GERPpp_NR'])) & np.array(indel_index), 'GERPpp_NR'] = patient.loc[patient['GERPpp_NR'] != -100, 'GERPpp_NR'].describe()['mean']
    #patient.loc[pd.isna(patient['GERPpp_NR']), 'GERPpp_NR'] = patient['GERPpp_NR'].describe()['min']


    #FATHMM_score
    patient['FATHMM_score'] = patient['FATHMM_score']
    patient.loc[patient['FATHMM_score'] == '-','FATHMM_score'] = np.NaN
    for i in patient[~patient['FATHMM_score'].isna()].index:
        score_list = str(patient.loc[i, 'FATHMM_score']).split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'FATHMM_score'] = np.NaN
        else:
            patient.loc[i, 'FATHMM_score'] = min(score_list)

    if np.all(pd.isna(patient['FATHMM_score'])):
        # Single variant query
        patient.loc[(pd.isna(patient['FATHMM_score'])) & np.array(indel_index), 'FATHMM_score'] = feature_stat.loc['FATHMM_score', '50%']
        patient.loc[pd.isna(patient['FATHMM_score']), 'FATHMM_score'] = feature_stat.loc['FATHMM_score', '50%']
        patient['FATHMM_score'] = patient['FATHMM_score'].astype('float64')
    else:
        # Patient whole VCF
        patient['FATHMM_score'] = patient['FATHMM_score'].astype('float64')
        patient.loc[(pd.isna(patient['FATHMM_score'])) & np.array(indel_index), 'FATHMM_score'] = patient['FATHMM_score'].describe()['50%']
        patient.loc[pd.isna(patient['FATHMM_score']), 'FATHMM_score'] = patient['FATHMM_score'].describe()['50%']


    #M_CAP_score
    patient.loc[patient['M_CAP_score'] == '-', 'M_CAP_score'] = np.NaN
    if np.all(pd.isna(patient['M_CAP_score'])):
        # Single variant query
        patient.loc[(pd.isna(patient['M_CAP_score'])) & np.array(indel_index), 'M_CAP_score'] = feature_stat.loc['M_CAP_score', 'mean']
        patient.loc[pd.isna(patient['M_CAP_score']), 'M_CAP_score'] = feature_stat.loc['M_CAP_score', 'mean']
        patient['M_CAP_score'] = patient['M_CAP_score'].astype('float64')
    else:
        # Patient whole VCF
        patient['M_CAP_score'] = patient['M_CAP_score'].astype('float64')
        patient.loc[(pd.isna(patient['M_CAP_score'])) & np.array(indel_index), 'M_CAP_score'] = patient['M_CAP_score'].describe()['mean']
        patient.loc[pd.isna(patient['M_CAP_score']), 'M_CAP_score'] = patient['M_CAP_score'].describe()['mean']


    #MutationAssessor_score
    patient['MutationAssessor_score'] = patient['MutationAssessor_score']
    patient.loc[patient['MutationAssessor_score'] == '-','MutationAssessor_score'] = np.NaN
    for i in patient[~patient['MutationAssessor_score'].isna()].index:
        score_list = str(patient.loc[i, 'MutationAssessor_score']).split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'MutationAssessor_score'] = np.NaN
        else:
            patient.loc[i, 'MutationAssessor_score'] = max(score_list)

    if np.all(pd.isna(patient['MutationAssessor_score'])):
        # Single variant query
        patient.loc[(pd.isna(patient['MutationAssessor_score'])) & np.array(indel_index), 'MutationAssessor_score'] = feature_stat.loc['MutationAssessor_score', '50%']
        patient.loc[pd.isna(patient['MutationAssessor_score']), 'MutationAssessor_score'] = feature_stat.loc['MutationAssessor_score', '50%']
        patient['MutationAssessor_score'] = patient['MutationAssessor_score'].astype('float64')
    else:
        # Patient whole VCF
        patient['MutationAssessor_score'] = patient['MutationAssessor_score'].astype('float64')
        patient.loc[(pd.isna(patient['MutationAssessor_score'])) & np.array(indel_index), 'MutationAssessor_score'] = patient['MutationAssessor_score'].describe()['50%']
        patient.loc[pd.isna(patient['MutationAssessor_score']), 'MutationAssessor_score'] = patient['MutationAssessor_score'].describe()['50%']


    #ESP6500_AA_AF
    patient.loc[patient['ESP6500_AA_AF'] == '-', 'ESP6500_AA_AF'] = 0.0
    patient['ESP6500_AA_AF'] = patient['ESP6500_AA_AF'].astype('float64')


    #ESP6500_AA_AF
    patient.loc[patient['ESP6500_EA_AF'] == '-', 'ESP6500_EA_AF'] = 0.0
    patient['ESP6500_EA_AF'] = patient['ESP6500_EA_AF'].astype('float64')


    #hom
    patient.loc[patient['hom'] == '-','hom'] = np.NaN
    for i in patient[~patient['hom'].isna()].index:
        score_list = str(patient.loc[i, 'hom']).split(',')
        score_list = [float(i) for i in score_list if i!='-' and i!='.']
        if score_list == []:
            patient.loc[i, 'hom'] = np.NaN
        else:
            patient.loc[i, 'hom'] = max(score_list)
    patient['hom'] = patient['hom'].fillna(0)
    patient['hom'] = patient['hom'].astype('float64')

    patient['hgmd_rs'] = patient['hgmd_rs'].apply(lambda x: x.split(",")[0] if type(x) == str else str(x))

    patient.loc[patient['hgmd_rs'] == '-','hgmd_rs'] = 0

    patient['hgmd_rs'] = patient['hgmd_rs'].astype('float64')


    #spliceAImax
    patient.loc[patient['spliceAImax'] == '-', 'spliceAImax'] = np.NaN
    if np.all(pd.isna(patient['spliceAImax'])):
        # Single variant query
        patient.loc[pd.isna(patient['spliceAImax']), 'spliceAImax'] = feature_stat.loc['spliceAImax', 'mean']
        patient['spliceAImax'] = patient['spliceAImax'].astype('float64')
    else:
        # Patient whole VCF
        patient['spliceAImax'] = patient['spliceAImax'].astype('float64')
        patient.loc[pd.isna(patient['spliceAImax']), 'spliceAImax'] = patient.loc[~pd.isna(patient['spliceAImax']), 'spliceAImax'].describe()['mean']


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
    CLN_STAT = ['-', 'no_assertion_provided', 'no_assertion_criteria_provided','no_assertion_for_the_individual_variant',
                'criteria_provided,_single_submitter', 'criteria_provided,_conflicting_interpretations',
                'criteria_provided,_multiple_submitters,_no_conflicts', 'reviewed_by_expert_panel', 'practice_guideline']

    CLN_STAT_SCORE = [0,0,0,0,1,1,2,3,4]

    for i in range(len(CLN_STAT)):
        patient.loc[patient['nc_CLNREVSTAT'] == CLN_STAT[i], 'nc_CLNREVSTAT'] = CLN_STAT_SCORE[i]
        patient.loc[patient['c_CLNREVSTAT'] == CLN_STAT[i], 'c_CLNREVSTAT'] = CLN_STAT_SCORE[i]

    patient['nc_CLNREVSTAT'] = patient['nc_CLNREVSTAT'].astype('float64')
    patient['c_CLNREVSTAT'] = patient['c_CLNREVSTAT'].astype('float64')


    #nc_RANKSCORE
    patient.loc[patient['nc_RANKSCORE'] == '-', 'nc_RANKSCORE'] = 0
    patient['nc_RANKSCORE'] = patient['nc_RANKSCORE'].astype('float64')

    #nc_RANKSCORE
    patient.loc[patient['c_RANKSCORE'] == '-', 'c_RANKSCORE'] = 0
    patient['c_RANKSCORE'] = patient['c_RANKSCORE'].astype('float64')

    #CLASS
    patient.loc[patient['CLASS'] == '-', 'CLASS'] = 0
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
    tier=tier_file.copy()
    tier_vars = ['IMPACT.from.Tier','TierAD','TierAR','TierAR.adj',
                 'No.Var.HM','No.Var.H','No.Var.M','No.Var.L',
                 'AD.matched','AR.matched',
                 'recessive', 'dominant'] #,
                 #'OMIM.Inheritance.both', 'OMIM.Inheritance.others']
    """ CL: deleted
    # recessive and dominant is enough
    tier['OMIM.Inheritance.recessive'] = 0
    tier['OMIM.Inheritance.dominant'] = 0
    tier['OMIM.Inheritance.both'] = 0
    tier['OMIM.Inheritance.others'] = 0
    tier.loc[tier['OMIM.Inheritance'] == 'recessive', 'OMIM.Inheritance.recessive'] = 1
    tier.loc[tier['OMIM.Inheritance'] == 'dominant', 'OMIM.Inheritance.dominant'] = 1
    tier.loc[tier['OMIM.Inheritance'] == 'both', 'OMIM.Inheritance.both'] = 1
    tier.loc[tier['OMIM.Inheritance'] == 'others', 'OMIM.Inheritance.others'] = 1
    """

    tier.loc[:,['TierAD','TierAR','TierAR.adj']] = -tier.loc[:,['TierAD','TierAR','TierAR.adj']]
    tier = tier.groupby(['Uploaded_variation'], sort=False)[tier_vars].max()
    tier.loc[:,['TierAD','TierAR','TierAR.adj']] = -tier.loc[:,['TierAD','TierAR','TierAR.adj']]
    
    patient = pd.concat([patient, tier], axis=1)
    patient.loc[patient['IMPACT.from.Tier'].isna(), 'IMPACT.from.Tier'] = 1
    patient.loc[patient['TierAD'].isna(), 'TierAD'] = 4
    patient.loc[patient['TierAR'].isna(), 'TierAR'] = 4
    patient.loc[patient['TierAR.adj'].isna(), 'TierAR.adj'] = 4

    for var in tier_vars[4:8]:
        patient.loc[patient[var].isna(), var] = 0

    patient.loc[patient['AD.matched'].isna(), 'AD.matched'] = 0
    patient.loc[patient['AR.matched'].isna(), 'AR.matched'] = 0

    patient.loc[patient['recessive'].isna(), 'recessive'] = 0
    patient.loc[patient['dominant'].isna(), 'dominant'] = 0
    #patient.loc[patient['OMIM.Inheritance.both'].isna(), 'OMIM.Inheritance.both'] = 0 ##CL: deleted
    #patient.loc[patient['OMIM.Inheritance.others'].isna(), 'OMIM.Inheritance.others'] = 0 ##CL: deleted

    return patient







