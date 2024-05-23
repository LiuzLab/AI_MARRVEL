#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np




def load_raw_matrix(score):
    #score_file = '/out/rami-test/r1_%s_scores.txt'%(ID) #'/out/rami-test/r1_%s_scores.txt'%(ID)
    #score = pd.read_csv(score_file, sep='\t')
    
    raw_features = ["chrom", "pos", "ref", 
                    "alt",'varId','varId_dash','zyg','geneSymbol','geneEnsId','gnomadAF','gnomadAFg','CADD_phred',
                    'CADD_PHRED','GERPpp_RS','GERPpp_NR','DANN_score','FATHMM_pred','FATHMM_score',
                    'Polyphen2_HDIV_score','Polyphen2_HVAR_score','REVEL_score','SIFT_score',
                    'fathmm_MKL_coding_score','LRT_score','LRT_Omega','phyloP100way_vertebrate',
                    'M_CAP_score','MutationAssessor_score','MutationTaster_score','ESP6500_AA_AF',
                    'ESP6500_EA_AF','symptomName','omimSymptomSimScore','omimSymMatchFlag','hgmdSymptomScore',
                    'hgmdSymptomSimScore','hgmdSymMatchFlag','clinVarSymMatchFlag','gnomadGeneZscore',
                    'gnomadGenePLI','gnomadGeneOELof','gnomadGeneOELofUpper','IMPACT','Consequence',
                    'omimGeneFound','omimVarFound','hgmdGeneFound','hgmdVarFound',
                    'clinVarVarFound','clinVarGeneFound','clinvarTotalNumVars',
                    'clinvarNumP','clinvarNumLP','clinvarNumLB','clinvarNumB','clinvarSignDesc',
                    'clinvarCondition','dgvVarFound','decipherVarFound','curationScoreHGMD','curationScoreOMIM',
                    'curationScoreClinVar','conservationScoreDGV','conservationScoreGnomad','conservationScoreOELof',
                    'hom','hgmd_rs','clin_dict','clin_PLP','clin_PLP_perc','spliceAImax','clin_code','hgmd_id',
                    'hgmd_CLASS','rsId', 'HGVSc', 'HGVSp', 'phenoList', 'phenoInhList']
    
    return score.loc[:,raw_features].copy()




def omimCurate(vardf):
    '''
    When related value changed,

    change the following value:
        omimSymptomSimScore
        omimSymMatchFlag
        curationScoreOMIM

    according to:
        omimGeneFound
        omimVarFound

    '''
    print('\nIn OMIM curation now...')

    # change omimSymptomSimScore according to wether gene found in OMIM database to avoid info leaking
    zero_score_bool = vardf['omimGeneFound']==0
    vardf.loc[zero_score_bool, 'omimVarFound'] = 0
    vardf.loc[zero_score_bool, 'omimSymptomSimScore'] = 0.0

    # change omimSymMatchFlag accordingly
    flag_0_bool = vardf['omimSymptomSimScore'] < 0.2
    vardf.loc[flag_0_bool, 'omimSymMatchFlag'] = 0
    vardf.loc[~flag_0_bool, 'omimSymMatchFlag'] = 1

    # change final curation score
    high_bool = (vardf['omimVarFound']==1) & (vardf['omimSymMatchFlag']==1)
    medi_bool = ((vardf['omimVarFound']==1) & (vardf['omimSymMatchFlag']==0)) | ((vardf['omimVarFound']==0) & (vardf['omimGeneFound']==1) & (vardf['omimSymMatchFlag']==1))

    vardf['curationScoreOMIM'] = 'Low'
    vardf.loc[high_bool, 'curationScoreOMIM'] = 'High'
    vardf.loc[medi_bool, 'curationScoreOMIM'] = 'Medium'

    return vardf




def hgmdCurate(vardf):
    '''
    When related value changed,

    change the following value:
        hgmdSymptomScore
        hgmdSymptomSimScore
        hgmdSymMatchFlag
        curationScoreHGMD

    according to:
        hgmdGeneFound
        hgmdVarFound

    '''
    print('\nIn HGMD curation now...')

    # change hgmdSymptomSimScore according to wether gene found in HGMD database to avoid info leaking
    zero_score_bool = vardf['hgmdGeneFound']==0
    vardf.loc[zero_score_bool, 'hgmdVarFound'] = 0
    vardf.loc[zero_score_bool, 'hgmdSymptomSimScore'] = 0.0

    # change hgmdSymptomScore according to wether variant found in HGMD database to avoid info leaking
    zero_score_bool = vardf['hgmdVarFound']==0
    vardf.loc[zero_score_bool, 'hgmdSymptomScore'] = 0.0

    # change hgmdSymMatchFlag accordingly
    score_numeric = vardf['hgmdSymptomSimScore'].copy()
    score_numeric[score_numeric=='-'] = np.NaN
    score_numeric = score_numeric.astype(float)
    flag_0_bool = score_numeric < 0.2
    vardf.loc[flag_0_bool, 'hgmdSymMatchFlag'] = 0
    vardf.loc[~flag_0_bool, 'hgmdSymMatchFlag'] = 1

    # change final curation score
    high_bool = (vardf['hgmdVarFound']==1) & (vardf['hgmdSymMatchFlag']==1)
    medi_bool = ((vardf['hgmdVarFound']==1) & (vardf['hgmdSymMatchFlag']==0)) | ((vardf['hgmdVarFound']==0) & (vardf['hgmdGeneFound']==1) & (vardf['hgmdSymMatchFlag']==1))

    vardf['curationScoreHGMD'] = 'Low'
    vardf.loc[high_bool, 'curationScoreHGMD'] = 'High'
    vardf.loc[medi_bool, 'curationScoreHGMD'] = 'Medium'

    return vardf





def clinvarCurate(vardf):
    '''
    When related value changed,

    change the following value:
        clinVarSymMatchFlag
        curationScoreClinVar

    according to:
        clinVarVarFound
        clinVarGeneFound
        clinvarSignDesc

    '''
    print('\nIn ClinVar curation now...')

    # change clinVarSymMatchFlag
    vardf.loc[vardf['clinVarGeneFound']==0, 'clinVarVarFound'] = 0
    one_bool = (vardf['clinVarSymMatchFlag']==1) & ((vardf['clinVarVarFound']==1) | (vardf['clinVarGeneFound']==1))
    vardf.loc[~one_bool, 'clinVarSymMatchFlag'] = 0

    # change final curation score

    pathList=['Pathogenic','Likely pathogenic','Pathogenic, Affects','Pathogenic/Likely pathogenic, other',
                'Pathogenic/Likely pathogenic','Pathogenic/Likely pathogenic, drug response',
                'Pathogenic/Likely pathogenic, risk factor','Likely pathogenic, drug response','Likely pathogenic, risk factor','Likely pathogenic, association',
                'Likely pathogenic, other','Pathogenic, association, protective','Pathogenic, Affects','Pathogenic, association','Pathogenic, other',
                'Pathogenic, protective','Pathogenic, protective, risk factor','Pathogenic, risk factor',
                'Pathogenic/Likely pathogenic, other','Pathogenic/Likely pathogenic, risk factor']
    
    #pathList=['Pathogenic','Likely pathogenic', 'Pathogenic/Likely pathogenic']
    benignList=['Benign', 'Likely benign', 'Benign/Likely benign', 'Benign, association', 'Benign, drug response',
               'Benign, other', 'Benign, protective', ' Benign/Likely benign, Affects', 'Benign/Likely benign, association',
               'Benign/Likely benign, drug response', 'Benign/Likely benign, drug response, risk factor',
               'Benign/Likely benign, other', 'Benign/Likely benign, protective', 'Benign/Likely benign, protective, risk factor',
               'Benign/Likely benign, risk factor', 'Likely benign, drug response, other','Likely benign, other',
               'Likely benign, other, risk factor', 'Likely benign, risk factor']
                                             
    #right logic
    #path_bool = vardf['clinvarSignDesc'].str.contains('Pathogenic', case=False) & (~vardf['clinvarSignDesc'].str.contains('Conflicting', case=False))
    #bngn_bool = vardf['clinvarSignDesc'].str.contains('Benign', case=False) & (~vardf['clinvarSignDesc'].str.contains('Conflicting', case=False))
    # wrong but current logic
    path_bool = vardf['clinvarSignDesc'].isin(pathList)
    bngn_bool = vardf['clinvarSignDesc'].isin(benignList)

    high_bool = (vardf['hgmdVarFound']==1) & path_bool & (vardf['hgmdSymMatchFlag']==1)
    medi_bool = ((vardf['hgmdVarFound']==1) & path_bool & (vardf['hgmdSymMatchFlag']==0)) | \
                ((vardf['hgmdVarFound']==1) & (~path_bool) & (~bngn_bool) & (vardf['hgmdSymMatchFlag']==1) & (vardf['hgmdSymMatchFlag']==1)) | \
                ((vardf['hgmdVarFound']==0) & (vardf['hgmdGeneFound']==1) & (vardf['clinVarSymMatchFlag']==1))


    vardf['curationScoreHGMD'] = 'Low'
    vardf.loc[high_bool, 'curationScoreHGMD'] = 'High'
    vardf.loc[medi_bool, 'curationScoreHGMD'] = 'Medium'

    return vardf




def conservationCurate(vardf):
    '''
    When related value changed,

    change the following value:
        conservationScoreDGV
        conservationScoreGnomad
        conservationScoreOELof

    according to:
        gnomadAF
        gnomadAFg
        dgvVarFound
        gnomadGeneOELofUpper

    '''

    # gnomad
    gnomadAFVal = vardf['gnomadAF'].copy()
    gnomadAFVal = np.array([ getValFromStr(str(i), 'min') for i in gnomadAFVal ])
    #gnomadAFVal[gnomadAFVal=='-'] = np.NaN
    gnomadAFVal = gnomadAFVal.astype(float)

    gnomadAFgVal = vardf['gnomadAFg'].copy()
    gnomadAFgVal = np.array([ getValFromStr(str(i), 'min') for i in gnomadAFgVal ])
    #gnomadAFgVal[gnomadAFgVal=='-'] = np.NaN
    gnomadAFgVal = gnomadAFgVal.astype(float)

    low_bool = (gnomadAFVal>=0.01) | (gnomadAFgVal>=0.01)
    vardf.loc[low_bool, 'conservationScoreGnomad'] = 'Low'
    vardf.loc[~low_bool, 'conservationScoreGnomad'] = 'High'


    # conservationScoreDGV
    low_bool = (vardf['conservationScoreDGV'] == 'Low') & (vardf['dgvVarFound']==1)
    vardf['conservationScoreDGV'] = 'High'
    vardf.loc[low_bool,'conservationScoreDGV'] = 'Low'

    # conservationScoreOELof
    gnomadGeneOELofUpperVal = vardf['gnomadGeneOELofUpper'].copy()
    gnomadGeneOELofUpperVal[gnomadGeneOELofUpperVal=='-'] = np.NaN
    gnomadGeneOELofUpperVal = gnomadGeneOELofUpperVal.astype(float)

    high_bool = gnomadGeneOELofUpperVal < 0.35
    vardf['conservationScoreOELof'] = 'Low'
    vardf.loc[high_bool, 'conservationScoreOELof'] = 'High'

    return vardf

def getValFromStr(valStr: str, select: str = 'min'):
    """
    Function to convert string to float,
    and takes care of situation when multiple values exist

    Some test example:
    >> a = '-'
    >> b = '-,0.012'
    >> c = '0.0003'
    >> d = '6.38284e-05'
    >> e = '6.38284e-05,3.19142e-05'
    >> f = '-,-'
    >> getValFromStr(a)
    '-'
    >> getValFromStr(b)
    '-'
    >> getValFromStr(c)
    0.0003
    >> getValFromStr(d)
    6.38284e-05
    >> getValFromStr(e)
    3.19142e-05
    >> getValFromStr(f)
    '-'
    """
    select_method = {'min':min, 'max':max}
    vals = valStr.split(',')
    if '-' in vals:
        return np.NaN
    else:
        vals = [float(i) for i in vals]
        return select_method[select](vals)




