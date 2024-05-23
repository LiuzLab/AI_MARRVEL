## UPDATE
#     Changes made by Chaozhong noted as #CL

from utils_1 import  *
import numpy as np
import time
import re

def getDecipherUsingMarrvelFlatFile(varObj, decipherDf):
    '''
    Get Decipher annotation using marrvel flat file
    Params:
    varObj:variant object
    decipherDf:dataframe of input flat file
    Return:
    Decipher list
    '''
    #print('in DECIPHER call')
    decipherDictList=[]
    decipherDeletionObsList=[]
    decipherStudyList=[]
    decipherVarFound=0
    deletionObs='-'
    #get the varaint object info from varObj
    chromVal=int(varObj.chrom)
    posVal=int(varObj.pos)
    startVal=int(varObj.start)
    stopVal=int(varObj.stop)

    # CL 03-14-2023: changed column names to be compatible with hg38
    #vals=decipherDf[ ( decipherDf['hg19Chr'] == chromVal ) & ( decipherDf['hg19Start']==startVal ) & (decipherDf['hg19Stop']==stopVal) ]
    vals=decipherDf[ ( decipherDf['Chr'] == chromVal ) & ( decipherDf['Start']==startVal ) & (decipherDf['Stop']==stopVal) ]
    numRows=len(vals.index)

    if numRows>0:
        decipherVarFound=1
        deletionObs=vals.iloc[0]['deletion.obs']
        decipherDeletionObsList.append(deletionObs)

    #print('\tchrom:', chromVal,'posVal:', posVal,'start:', startVal,'stopVal:', stopVal)
    #print('\tdecipherVarFound:',decipherVarFound,'decipherDeletionObs:', deletionObs)
    retList=[decipherDictList,decipherDeletionObsList,decipherStudyList, decipherVarFound]
    return retList

def getDGVUsingMarrvelFlatFile(varObj, dgvDf):
    '''
    Get DGV annotation using marrvel flat file
    Params:
    varObj:variant object
    dgvDf:dataframe of input flat file
    Return:
    DGV list
    '''
    #print('in DGV call')
    dgvDictList=[]
    typeList=[]
    subtypeList=[]
    dgvVarFound=0
    dgvType='-'
    dgvSubtype='-'
    
    chromVal=int(varObj.chrom)
    posVal=int(varObj.pos)
    startVal=int(varObj.start)
    stopVal=int(varObj.stop)

    #CL 03-14-2023: changed column names to be compatible with hg38
    #vals=dgvDf[ ( dgvDf['hg19Chr'] == chromVal ) & ( dgvDf['hg19Start']<=startVal ) & (dgvDf['hg19Stop']>=stopVal) ]
    vals=dgvDf[ ( dgvDf['Chr'] == chromVal ) & ( dgvDf['Start']<=startVal ) & (dgvDf['Stop']>=stopVal) ]
    numRows=len(vals.index)

    if numRows>0:
        dgvVarFound=1
        #print('\tnumrows:',numRows)
        #print('\t type of vals:', type(vals))
        #print('\tvals:', vals)
        dgvType=vals.iloc[0]['type']
        dgvSubtype=vals.iloc[0]['subType']

    #print('\tchrom:', chromVal,'posVal:', posVal,'start:', startVal,'stopVal:', stopVal)
    #print('\tdgvVarFound:',dgvVarFound,'dgvType:', dgvType, 'dgvsubtype:', dgvSubtype)
    typeList.append(dgvType)
    subtypeList.append(dgvSubtype)
    retList=[dgvDictList,typeList,subtypeList,dgvVarFound]
    return retList

def getGnomadGeneMetricFlatFile(varObj, gnomadMetricsGeneDf):
    '''
    Get GnomadGeneMetric from flat file
    Params:
    varObj:variant object
    gnomadMetricsGeneDf:dataframe of input flat file
    Return:
    gnomadMetricsGeneDf list
    '''
    #print('in gnomadGeneMetricFlatFile call')
    #check if gene in file
    vals=gnomadMetricsGeneDf[gnomadMetricsGeneDf['gene'] == varObj.geneSymbol]
    numRows=len(vals.index)

    if numRows>0:#pLI, oe_lof, oe_lof_upper,mis_z
        gnomadGeneZscore=vals.iloc[0]['mis_z']
        gnomadGenePLI=vals.iloc[0]['pLI']
        gnomadGeneOELof=vals.iloc[0]['oe_lof']
        gnomadGeneOELofUpper=vals.iloc[0]['oe_lof_upper']
    else:
        #get the values
        gnomadGeneZscore='-'
        gnomadGenePLI='-'
        gnomadGeneOELof='-'
        gnomadGeneOELofUpper='-'

    retList=[gnomadGeneZscore, gnomadGenePLI, gnomadGeneOELof, gnomadGeneOELofUpper]
    #print('\t gene symbol:', varObj.geneSymbol, 'zScore:', gnomadGeneZscore, 'gnomadGenePLI:',gnomadGenePLI,'gnomadGeneOELof:', gnomadGeneOELof,'gnomadGeneOELofUpper:', gnomadGeneOELofUpper)
    #print('\t gene symbol:', varObj.geneSymbol, 'zScore:', gnomadGeneZscore)
    return retList

def getConservationScore(varObj, diseaseInh):
    '''
    Get conservation score
    Params:
    varObj:variant object
    diseaseInh: the disease inheritance
    Return:
    return three conservation scores: conservationScoreGnomad, conservationScoreDGV, conservationScoreOELof
    '''
    #print('in conservation score call')
    #create conservation scores
    #set up thresholds
    gnomadGFreqCut=0.01#genome g
    gnomadEFreqCut=0.01#exome e
    # GERP ranges
    gerpPPCut=4
    #LRT cut
    #phylo cut

    #gnomad
    conservationScoreGnomad='-'
    gnomadAFVal = getValFromStr(str(varObj.gnomadAF), 'min')
    gnomadAFgVal = getValFromStr(str(varObj.gnomadAFg), 'min')
    if gnomadAFVal!='-' and gnomadAFgVal!='-':
        #gnomadAFVal=float(varObj.gnomadAF)
        #gnomadAFgVal=float(varObj.gnomadAFg)
        if gnomadAFVal<0.01 and gnomadAFgVal<0.01:
            conservationScoreGnomad='High'
        else:
            conservationScoreGnomad='Low'
    else:
        conservationScoreGnomad='-'

    #DGV
    conservationScoreDGV='-'
    if 'deletion' in varObj.dgvSubtypeList or 'loss' in varObj.dgvSubtypeList:
        conservationScoreDGV='Low'
    else:
        conservationScoreDGV='High'

    #gene O/E score
    conservationScoreOELof='-'
    if varObj.gnomadGeneOELofUpper!='-':
        gnomadGeneOELofUpperVal=float(varObj.gnomadGeneOELofUpper)
        if gnomadGeneOELofUpperVal<0.35:
            conservationScoreOELof='High'
        else:
            conservationScoreOELof='Low'

    #return 
    retList=[conservationScoreGnomad, conservationScoreDGV, conservationScoreOELof]
    return retList

def getValFromStr(valStr: str, select: str = 'min'):
    """
    Function to convert string to float,
    and takes care of situation when multiple values exist
    """
    select_method = {'min':min, 'max':max}
    vals = valStr.split(',')
    if '-' in vals:
        return '-'
    else:
        vals = [float(i) for i in vals]
        return select_method[select](vals)


    

