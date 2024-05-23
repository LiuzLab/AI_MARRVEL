## UPDATE
#     Changes made by Chaozhong noted as #CL

from utils_1 import  *
from utils_for_marrvel_flatfile_module_2 import *
import numpy as np
import time
import re

def getOMIMUsingMarrvelFlatFile(varObj, omimGeneList, omimAlleleList):
    '''
    function to get OMIM info using marrrvel flat file
    Params:a varaint object read from VEP annotation
    
    Returns:
    List of OMIM annotations
    '''
    #print('In OMIM call:')
    #print('\tvar:', varObj.varId_dash, 'rsId:', varObj.rsId)
    #check if rsId has a comma  in it
    inputSnpList=[]
    if ',' in varObj.rsId:
        inputSnpList = varObj.rsId.split(',')
    else:
        inputSnpList=varObj.rsId
    #print('\tinputSnpList:', inputSnpList)
    varFound=0
    geneFound=0
    omimDict={}
    omimGeneDict={}
    omimAlleleDict={}
    phenoList=[]
    phenoInhList=[]
    phenoMimList=[]
    #check gene
    #keys: dict_keys(['phenotypes', 'allelicVariants', 'mimNumber', 'status', 'title', 'description', 'geneEntrezId', 'geneSymbol'])
    if any(d['geneSymbol'] == varObj.geneSymbol for d in omimGeneList):
        #print('\tgene:', varObj.geneSymbol, 'found')
        geneFound=1
        omimGeneDict = next((sub for sub in omimGeneList if sub['geneSymbol'] == varObj.geneSymbol), None)
        #print('\tgeneSymbol from dict:', omimGeneDict['geneSymbol'])
        #print('\talleleic:', type(omimGeneDict['allelicVariants']), ' len:', len(omimGeneDict['allelicVariants']) )
        #check if snp found
        snpList=[]
        for a in omimGeneDict['allelicVariants']:
            #print('a:', a)
            #print('type:', type(a))
            if 'dbSnps' in a:
                snpList.append(a['dbSnps'])
        #print('\tsnpList:', snpList)
        #print('\tlen snpList:', len(snpList))
        #check if input snpID matches the OMIM one
        set1 = set(inputSnpList) 
        set2 = set(snpList) 
        if set1.intersection(set2): 
            varFound=1 
        else: 
            varFound=0
           
        #get disease info from OMIM
        #print('\tphenotypes:', type(omimGeneDict['phenotypes']), ' len:', len(omimGeneDict['phenotypes']) )
        for a in omimGeneDict['phenotypes']:
          #print('type:', type(a))
          pheno=a['phenotype']
          if 'phenotypeMimNumber' in a:
            phenoMim=a['phenotypeMimNumber']
          else:
            phenoMim='-'
          if 'phenotypeInheritance' in a:
            phenoInh=a['phenotypeInheritance']
          else:
            phenoInh='-'
          phenoList.append(pheno)
          phenoInhList.append(phenoInh)
          phenoMimList.append(str(phenoMim))
          #print('phenotype:', pheno,phenoMim,phenoInh)
        
    #print('\tvarFound:', varFound)
    #print('\tphenoList:', phenoList)
    #print('\tphenoInhList:', phenoInhList)
    #print('\tphenoMimList:', phenoMimList)

    retList=[varFound, geneFound, omimDict, omimGeneDict, omimAlleleDict, phenoList, phenoInhList, phenoMimList]
    return retList

def getClinVarUsingMarrvelFlatFile(varObj,clinvarAlleleDf,clinvarGeneDf):
    '''
    function to get clinvar info using marrrvel flat file
    Params:a varaint object read from VEP annotation
    
    Returns:
    List of clinvar annotations
    '''
    #print('in clinvar using flatfile')
    #print('\tvar:',varObj.varId_dash)
    varFound=0
    varDict={}
    geneFound=0
    geneDict={}
    clinvarTotalNumVars=0
    clinvarNumP=0
    clinvarNumLP=0
    clinvarNumLB=0
    clinvarNumB=0
    ## CL: no use but kept
    clinvarTitle=''
    clinvarSignDesc=''
    clinvarCondition=''

    #print('in function type clinavrDf:', type(clinvarAlleleDf) )
    chromVal=int(varObj.chrom)
    posVal=int(varObj.pos)

    """ #CL: Old implementation
    #using inetger column not index
    if 1:
        vals=clinvarAlleleDf.loc[(clinvarAlleleDf['chr']==chromVal ) &
                             (clinvarAlleleDf['start'] == posVal ) &
                             (clinvarAlleleDf['stop'] == posVal )]
        numRows=len(vals.index)
        print('columns numRows:', numRows)
    #using index
    if 0:
        idVal=str(chromVal)+'_'+str(posVal)+'_'+str(posVal)
        #vals=clinvarAlleleDf.loc[(clinvarAlleleDf['id1']==idVal)]
        try:
            vals=clinvarAlleleDf.loc[idVal]
            numRows=len(vals.index)
            print('index numRows:', numRows)
        except:
            numRows=0

    print('\tcheck var numRows:', numRows)
    if numRows > 0:
        varFound=1
        print('\tclinvar var found')
        print('\t clinvar vals:', vals)
        clinvarTitle=vals.iloc[0]['title']
        clinvarSignDesc=vals.iloc[0]['significance.description']
        clinvarCondition=vals.iloc[0]['condition']
        print('\ttitle:', clinvarTitle, 'signDes:', clinvarSignDesc, 'condition:', clinvarCondition)
    else:
        varFound=0
    """

    #CL: check if var annotated in clinvar using VEP clinvar.vcf.gz custom annotation
    if varObj.clinvar_AlleleID != '-':
        varFound = 1
        #print('\tclinvar var found')
    else:
        varFound = 0



    #check if gene annotated in clinvar using clinvarGeneDf
    geneSymbol=varObj.geneSymbol
    #print('\tgene symbol:', geneSymbol)
    try:
        vals=clinvarGeneDf.loc[geneSymbol]
        numRows=len(vals.index)
    except:
        numRows=0
    #print('\tcheck gene numRows:', numRows)
    if numRows>0:
        #print('\tclinvar gene found')
        geneFound=1
        clinvarTotalNumVars=vals['totalClinvarVars']
        clinvarNumP=vals['P']
        clinvarNumLP=vals['LP']
        clinvarNumLB=vals['LB']
        clinvarNumB=vals['B']
    #return
    retList=[varFound,varDict, geneFound, geneDict, clinvarTotalNumVars,clinvarNumP,clinvarNumLP,clinvarNumLB,
             clinvarNumB, clinvarTitle, clinvarSignDesc,clinvarCondition]
    return retList

def getHGMDUsingFlatFile(varObj,hgmdDf):
    '''
    function to get HGMD from local flat file
    Params:
    varObj:a varaint object read from VEP annotation
    hgmdDf: HGMD data frame read from local file (CL: now it refers to hgmdHPOScoreDf in main.py)
    
    Returns:
    List of HGMD annotations

    Update by CL:
    '''
    #print('\nin HGMD')
    #print('\tvar:', varObj.varId_dash, 'var-gene:', varObj.geneSymbol)
    HGMDDict={}
    hgmdGeneFound=0
    hgmdVarFound=0
    hgmdVarPhenIdList=[]
    hgmdVarHPOIdList=[]
    hgmdVarHPOStrList=[]
    chromVal=varObj.chrom
    posVal=int(varObj.pos)
    startVal=int(varObj.start)
    stopVal=int(varObj.stop)
    #print('\tpos type:',type(varObj.start),'chrom:', type(chromVal) )

    """
    #using int columns
    if 1:
        vals=hgmdDf[ ( hgmdDf['chromosome'] == chromVal ) & ( hgmdDf['startCoord']==startVal ) & (hgmdDf['endCoord']==stopVal) ]
        numRows=len(vals.index)
    #using index
    if 0:
        idVal=str(chromVal)+'_'+str(startVal)+'_'+str(stopVal)
        try:
            vals=hgmdVarDf.loc[idVal]
            numRows=len(vals.index)
            print('index numRows:', numRows)
        except:
            numRows=0

        print('\tvar numRows:', numRows)


    if numRows>0:
        hgmdVarFound=1
        print('\tnumrows:',numRows)
        print('\tvals:', vals)
        if 'phen_id' in vals:
            hgmdVarPhenIdList.extend(vals['phen_id'].tolist())
        if 'hpo_id' in vals:
            hgmdVarHPOIdList.extend(vals['hpo_id'].tolist())
        if 'hpo_str' in vals:
            hgmdVarHPOStrList.extend(vals['hpo_str'].tolist())
        print('\tvals:', vals)
    """
    # CL: check VarFound
    if varObj.hgmd_id != '-':
        hgmdVarFound = 1
    else:
        hgmdVarFound = 0

    """
    #check gene
    #vals=hgmdGeneDf[(hgmdGeneDf['gene']==varObj.geneSymbol)]
    print('\t1 HGMD geneSymbol:', varObj.geneSymbol)
    try:
        print('\t2 HGMD geneSymbol:', varObj.geneSymbol)
        vals=hgmdDf.loc[varObj.geneSymbol]
        #vals=hgmdDf[ ( hgmdDf['gene'] == varObj.geneSymbol ) ]
        numRows=len(vals.index)
    except:
        numRows=0

    print('\tHGMD gene found numRows:', numRows)
    if numRows>0:
        hgmdGeneFound=1
    """
    # CL: check geneFound
    if np.any(hgmdDf['gene_sym'].isin([varObj.geneSymbol])):
        hgmdGeneFound = 1
    else:
        hgmdGeneFound = 0


    #print('\thgmdVarFound:',hgmdVarFound,'hgmdGeneFound:',hgmdGeneFound,
    #      'hgmdVarPhenIdList:',hgmdVarPhenIdList,'hgmdVarHPOIdList:',hgmdVarHPOIdList,
    #      'hgmdVarHPOStrList:',hgmdVarHPOStrList)
    #return
    retList=[hgmdVarFound,hgmdGeneFound,hgmdVarPhenIdList,hgmdVarHPOIdList,hgmdVarHPOStrList]
    return retList

def getVarObjFromMarrvelFlatFile(varDf,genomeRef,clinvarGeneDf,clinvarAlleleDf, omimGeneList, omimAlleleList,hgmdDf):
    '''
    Function to get the marrvel annotaion if not done in parallel. Use parallel option better
    '''
    #dict_keys(['chr', 'pos', 'ref', 'alt', 'omim', 'clinvar', 'gtex', 'dgv',
    #'decipher', 'decipherDisease', 'gnomad', 'dbnsfp'])
    numVars=len(varDf.index)
    print('type:', type(varDf), 'len/numVars:', numVars)
    varObjList = [ Variant() for i in range(numVars) ]

    debugFlag=0

    i=0
    for index, row in varDf.iterrows():
        print('type of row:', type(row))
        transcriptId=row.Feature
        #s=row.Uploaded_variation.split('_') '1_10204_-/T'
        s=row[0].split('_')
        #print('s:', s)
        chrom=s[0]
        pos=int(s[1])
        tmp=s[2]
        s=tmp.split('/')
        ref=s[0]
        alt=s[1]
        #get the start and stop from second column like '1:10203-10204'
        if '-' in row[1]:
            s=row[1].split(':')
            tmp=s[1]
            s=tmp.split('-')
            print('s:',s)
            start=int(s[0])
            stop=int(s[1])
        else:
            #start and stop the same
            s=row[1].split(':')
            start=int(s[1])
            stop=int(s[1])


        print('chrom:', chrom,'pos:',pos,'ref:',ref,'alt:',alt,'start:',start,'stop:',stop)
        #change chrom X and Y and MT to numbers
        if chrom=='X':
            chrom=23
        elif chrom=='Y':
            chrom=24
        elif chrom=='MT':
            chrom=25
        chrom=int(chrom)
        #if it is hg38 get its hg19 coordinates
        if genomeRef=='hg38':
            varObjList[i].hg38Chrom=chrom
            varObjList[i].hg38Pos=pos
            retList=gethg19LocFromHg38(chrom, pos)
            retList=[newChrom, newPos]
            varObjList[i].hg19Chrom=retList[0]
            varObjList[i].hg19Pos=retList[1]
            chrom=retList[0]
            pos=retList[1]
            #get the start
            retList=gethg19LocFromHg38(chrom, start)
            varObjList[i].start=int(retList[1])
            #get the stop
            retList=gethg19LocFromHg38(chrom, stop)
            varObjList[i].stop=int(retList[1])
        else:
            varObjList[i].hg19Chrom=chrom
            varObjList[i].hg19Pos=pos
            varObjList[i].chrom=chrom
            varObjList[i].pos=pos
            varObjList[i].start=start
            varObjList[i].stop=stop


        geneSymbol=row.SYMBOL
        print('gene:', geneSymbol)
        varObjList[i].geneSymbol=geneSymbol
        varObjList[i].CADD_phred=row.CADD_phred
        varObjList[i].CADD_PHRED=row.CADD_PHRED
        #assign
        varObjList[i].ref=ref
        varObjList[i].alt=alt
        varObjList[i].varId_dash='-'.join([str(chrom),str(pos),ref,alt])
        varId='_'.join([str(chrom),str(pos),ref,alt,transcriptId])
        varObjList[i].varId=varId
        varObjList[i].zyg=row.ZYG
        varObjList[i].geneEnsId=row.Gene
        varObjList[i].rsId=row.Existing_variation
        varObjList[i].GERPpp_RS=row.GERPpp_RS
        varObjList[i].featureType=row.Feature_type
        varObjList[i].gnomadAF=row.gnomAD_AF
        varObjList[i].gnomadAFg=row.gnomADg_AF
        varObjList[i].CLIN_SIG=row.CLIN_SIG
        varObjList[i].HGVSc=row.HGVSc
        varObjList[i].HGVSp=row.HGVSp
        #dbnsfp attributes 
        varObjList[i].GERPpp_NR=row.GERPpp_NR
        varObjList[i].DANN_score=row.DANN_score
        varObjList[i].FATHMM_pred=row.FATHMM_pred
        varObjList[i].FATHMM_score=row.FATHMM_score
        varObjList[i].GTEx_V8_gene=row.GTEx_V8_gene
        varObjList[i].GTEx_V8_tissue=row.GTEx_V8_tissue
        varObjList[i].Polyphen2_HDIV_score=row.Polyphen2_HDIV_score
        varObjList[i].Polyphen2_HVAR_score=row.Polyphen2_HVAR_score
        varObjList[i].REVEL_score=row.REVEL_score
        varObjList[i].SIFT_score=row.SIFT_score
        varObjList[i].clinvar_clnsig=row.clinvar_clnsig
        varObjList[i].fathmm_MKL_coding_score=row.fathmm_MKL_coding_score
        varObjList[i].LRT_score=row.LRT_score
        varObjList[i].LRT_Omega=row.LRT_Omega
        varObjList[i].phyloP100way_vertebrate=row.phyloP100way_vertebrate
        varObjList[i].M_CAP_score=row.M_CAP_score
        varObjList[i].MutationAssessor_score=row.MutationAssessor_score
        varObjList[i].MutationTaster_score=row.MutationTaster_score
        varObjList[i].ESP6500_AA_AC=row.ESP6500_AA_AC
        varObjList[i].ESP6500_AA_AF=row.ESP6500_AA_AF
        varObjList[i].ESP6500_EA_AC=row.ESP6500_EA_AC
        varObjList[i].ESP6500_EA_AF=row.ESP6500_EA_AF
        #dbnsfp
        
        #for module one will need OMIM, HGMD, clinvar
        #get OMIM
        if 1:
            print('\nReading OMIM')
            #varObjList[i].omimList=jsonDict['omim']
            #retList=[varFound, geneFound, omimDict, omimGeneDict, omimAlleleDict]
            omimRet=getOMIMUsingMarrvelFlatFile(varObjList[i], omimGeneList, omimAlleleList)
            varObjList[i].omimVarFound=omimRet[0]
            varObjList[i].omimGeneFound=omimRet[1]
            varObjList[i].omimDict=omimRet[2]
            varObjList[i].omimGeneDict=omimRet[3]
            varObjList[i].omimAlleleDict=omimRet[4]
            varObjList[i].phenoList=omimRet[5]
            varObjList[i].phenoInhList=omimRet[6]
            varObjList[i].phenoMimList=omimRet[7]
            print('OMIM res:')
            print('\tgeneFound:',varObjList[i].omimGeneFound,'varFound:',varObjList[i].omimVarFound )

        #get clinvar
        if 1:
            print('\nReading clinVar')
            #clinVarList=jsonDict['clinvar']
            #varObjList[i].clinVarList=clinVarList
            #retList=[varFound,varDict, geneFound, geneDict, clinvarTotalNumVars,clinvarNumP,clinvarNumLP,clinvarNumLB,clinvarNumB]
            clinVarRet=getClinVarUsingMarrvelFlatFile(varObjList[i], clinvarAlleleDf, clinvarGeneDf)
            varObjList[i].clinVarVarFound=clinVarRet[0]
            varObjList[i].clinVarVarDict=clinVarRet[1]
            varObjList[i].clinVarGeneFound=clinVarRet[2]
            varObjList[i].clinVarGeneDict=clinVarRet[3]
            varObjList[i].clinvarTotalNumVars=clinVarRet[4]
            varObjList[i].clinvarNumP=clinVarRet[5]
            varObjList[i].clinvarNumLP=clinVarRet[6]
            varObjList[i].clinvarNumLB=clinVarRet[7]
            varObjList[i].clinvarNumB=clinVarRet[8]
            varObjList[i].clinvarTitle= clinVarRet[9]
            varObjList[i].clinvarSignDesc=clinVarRet[10]
            varObjList[i].clinvarCondition=clinVarRet[11]
            print('clinVar res:')
            if debugFlag==1:
                print('\tgeneFound::',varObjList[i].clinVarGeneFound,'varFound:',varObjList[i].clinVarVarFound)
                print('\tnumVars:',varObjList[i].clinvarTotalNumVars,'numPathologic:',varObjList[i].clinvarNumP,'numBenign:',varObjList[i].clinvarNumB)
                print('\tsignDesc:', varObjList[i].clinvarSignDesc)

        #get HGMD annotation
        if 1:
            print('\nReading HGMD')
            hgmdRet=getHGMDUsingFlatFile(varObjList[i],hgmdDf)
            #hgmdVarFound,hgmdGeneFound,hgmdVarPhenIdList,hgmdVarHPOIdList,hgmdVarHPOStrList
            varObjList[i].hgmdVarFound=hgmdRet[0]
            varObjList[i].hgmdGeneFound=hgmdRet[1]
            varObjList[i].hgmdVarPhenIdList=hgmdRet[2]
            varObjList[i].hgmdVarHPOIdList=hgmdRet[3]
            varObjList[i].hgmdVarHPOStrList=hgmdRet[4]
            print('HGMD results:')
            print('\thgmdVarFound:',varObjList[i].hgmdVarFound,'hgmdGeneFound:',varObjList[i].hgmdGeneFound,
                  'hgmdVarPhenIdList:',varObjList[i].hgmdVarPhenIdList,'hgmdVarHPOIdList:',
                  varObjList[i].hgmdVarHPOIdList,
                  'hgmdVarHPOStrList:',varObjList[i].hgmdVarHPOStrList)

        #get DGV annotation
        if 1:
            print('Reading DGV')
            dgvList=jsonDict['dgv']
            #print('dgvList:', dgvList)
            dgvRet=getDGVUsingMarrvelAnnotScript(varObjList[i],dgvList)
            varObjList[i].dgvDictList=dgvRet[0]
            varObjList[i].dgvTypeList=dgvRet[1]
            varObjList[i].dgvSubtypeList=dgvRet[2]
        #get DECIPHER annotation
        if 1:
            print('Reading DECIPHER')
            decipherList=jsonDict['decipher']
            decipherRet=getDecipherUsingMarrvelAnnotScript(varObjList[i], decipherList)
            #
            varObjList[i].decipherDictList=decipherRet[0]
            varObjList[i].decipherDeletionObsList=decipherRet[1]
            varObjList[i].decipherStudyList=decipherRet[2]

        #break for debug
        i=i+1
        if debugFlag==1:
            if i==2:
                break


    return varObjList

def getVarObjListFromDict(varDictList, numVars):
    '''
    Make a var obj list after paralle processing 
    the list has three levels since from pool get a first level list which is equal to the number of cores used
    '''
    numProc=len(varDictList)
    print('numProc:', numProc)
    #initialize
    varObjList = [ Variant() for i in range(numVars) ]
    i=0
    for procList in varDictList:
        for varDict in procList:
            #print('i:', i)
            #for varDict in dictList:
            varObjList[i].chrom=varDict['chrom']
            varObjList[i].pos=varDict['pos']
            varObjList[i].start=varDict['start']
            varObjList[i].stop=varDict['stop']
            varObjList[i].ref=varDict['ref']
            varObjList[i].alt=varDict['alt']
            varObjList[i].hg19Chrom=varDict['hg19Chrom']
            varObjList[i].hg19Pos=varDict['hg19Pos']
            varObjList[i].varId=varDict['varId']
            varObjList[i].varId_dash=varDict['varId_dash']

            #varObjList[i].patientID=varDict['patientID']
            varObjList[i].zyg=varDict['ZYG']
            varObjList[i].geneSymbol=varDict['geneSymbol']
            varObjList[i].CADD_phred=varDict['CADD_phred']
            varObjList[i].CADD_PHRED=varDict['CADD_PHRED']
            varObjList[i].geneEnsId=varDict['Gene']
            varObjList[i].rsId=varDict['Existing_variation']
            varObjList[i].GERPpp_RS=varDict['GERPpp_RS']
            varObjList[i].featureType=varDict['Feature_type']
            varObjList[i].gnomadAF=varDict['gnomadAF']
            varObjList[i].gnomadAFg=varDict['gnomadAFg']
            varObjList[i].CLIN_SIG=varDict['CLIN_SIG']
            varObjList[i].HGVSc=varDict['HGVSc']
            varObjList[i].HGVSp=varDict['HGVSp']

            varObjList[i].VARIANT_CLASS=varDict['VARIANT_CLASS']
            varObjList[i].Feature=varDict['Feature']
            varObjList[i].hom=varDict['hom']
            varObjList[i].hgmd_rs=varDict['hgmd_rs']

            varObjList[i].clin_dict=varDict['clin_dict']
            varObjList[i].clin_PLP=varDict['clin_PLP']
            varObjList[i].clin_PLP_perc=varDict['clin_PLP_perc']
            varObjList[i].spliceAI=varDict['spliceAI']
            varObjList[i].spliceAImax=varDict['spliceAImax']


            #dbnsfp attributes 
            varObjList[i].GERPpp_NR=varDict['GERPpp_NR']
            varObjList[i].DANN_score=varDict['DANN_score']
            varObjList[i].FATHMM_pred=varDict['FATHMM_pred']
            varObjList[i].FATHMM_score=varDict['FATHMM_score']
            varObjList[i].GTEx_V8_gene=varDict['GTEx_V8_gene']
            varObjList[i].GTEx_V8_tissue=varDict['GTEx_V8_tissue']
            varObjList[i].Polyphen2_HDIV_score=varDict['Polyphen2_HDIV_score']
            varObjList[i].Polyphen2_HVAR_score=varDict['Polyphen2_HVAR_score']
            varObjList[i].REVEL_score=varDict['REVEL_score']
            varObjList[i].SIFT_score=varDict['SIFT_score']
            #varObjList[i].clinvar_clnsig=varDict['clinvar_clnsig']

            varObjList[i].clinvar_AlleleID = varDict['clinvar_AlleleID'] # Clinvar allele ID from clinvar.vcf.gz
            varObjList[i].clinvar_clnsig = varDict['clinvar_clnsig'] #CL: Clinvar SIG from clinvar.vcf.gz
            varObjList[i].clinvar_CLNREVSTAT = varDict['clinvar_CLNREVSTAT'] #CL: Clinvar STAT from clinvar.vcf.gz, for interface only
            varObjList[i].clinvar_CLNSIGCONF = varDict['clinvar_CLNSIGCONF'] #CL: Clinvar STAT from clinvar.vcf.gz
            varObjList[i].clin_code = varDict['clin_code'] #CL: feature for ai


            varObjList[i].fathmm_MKL_coding_score=varDict['fathmm_MKL_coding_score']
            varObjList[i].LRT_score=varDict['LRT_score']
            varObjList[i].LRT_Omega=varDict['LRT_Omega']
            varObjList[i].phyloP100way_vertebrate=varDict['phyloP100way_vertebrate']
            varObjList[i].M_CAP_score=varDict['M_CAP_score']
            varObjList[i].MutationAssessor_score=varDict['MutationAssessor_score']
            varObjList[i].MutationTaster_score=varDict['MutationTaster_score']
            varObjList[i].ESP6500_AA_AC=varDict['ESP6500_AA_AC']
            varObjList[i].ESP6500_AA_AF=varDict['ESP6500_AA_AF']
            varObjList[i].ESP6500_EA_AC=varDict['ESP6500_EA_AC']
            varObjList[i].ESP6500_EA_AF=varDict['ESP6500_EA_AF']
            #conservation
            varObjList[i].LRT_Omega=varDict['LRT_Omega']
            varObjList[i].LRT_score=varDict['LRT_score']
            varObjList[i].phyloP100way_vertebrate=varDict['phyloP100way_vertebrate']
            varObjList[i].IMPACT=varDict['IMPACT']
            varObjList[i].Consequence=varDict['Consequence']
            #DGV
            varObjList[i].dgvDictList=varDict['dgvDictList']
            varObjList[i].dgvTypeList=varDict['dgvTypeList']
            varObjList[i].dgvSubtypeList=varDict['dgvSubtypeList']
            varObjList[i].dgvVarFound=varDict['dgvVarFound']
            #Decipher
            varObjList[i].decipherDictList=varDict['decipherDictList']
            varObjList[i].decipherDeletionObsList=varDict['decipherDeletionObsList']
            varObjList[i].decipherStudyList=varDict['decipherStudyList']
            varObjList[i].decipherVarFound=varDict['decipherVarFound']
            #gnomadG
            varObjList[i].gnomadGeneZscore=varDict['gnomadGeneZscore']
            varObjList[i].gnomadGenePLI=varDict['gnomadGenePLI']
            varObjList[i].gnomadGeneOELof=varDict['gnomadGeneOELof']
            varObjList[i].gnomadGeneOELofUpper=varDict['gnomadGeneOELofUpper']

            varObjList[i].omimVarFound=varDict['omimVarFound']
            varObjList[i].omimGeneFound=varDict['omimGeneFound']
            varObjList[i].omimDict=varDict['omimDict']
            varObjList[i].omimGeneDict=varDict['omimGeneDict']
            varObjList[i].omimAlleleDict=varDict['omimAlleleDict']
            varObjList[i].phenoList=varDict['phenoList']
            varObjList[i].phenoInhList=varDict['phenoInhList']
            varObjList[i].phenoMimList=varDict['phenoMimList']

            varObjList[i].clinVarVarFound=varDict['clinVarVarFound']
            varObjList[i].clinVarVarDict=varDict['clinVarVarDict']
            varObjList[i].clinVarGeneFound=varDict['clinVarGeneFound']
            varObjList[i].clinVarGeneDict=varDict['clinVarGeneDict']
            varObjList[i].clinvarTotalNumVars=varDict['clinvarTotalNumVars']
            varObjList[i].clinvarNumP=varDict['clinvarNumP']
            varObjList[i].clinvarNumLP=varDict['clinvarNumLP']
            varObjList[i].clinvarNumLB=varDict['clinvarNumLB']
            varObjList[i].clinvarNumB=varDict['clinvarNumB']
            varObjList[i].clinvarTitle= varDict['clinvarTitle']
            varObjList[i].clinvarSignDesc=varDict['clinvarSignDesc']
            varObjList[i].clinvarCondition=varDict['clinvarCondition']

            varObjList[i].hgmdVarFound=varDict['hgmdVarFound']
            varObjList[i].hgmdGeneFound=varDict['hgmdGeneFound']
            varObjList[i].hgmdVarPhenIdList=varDict['hgmdVarPhenIdList']
            varObjList[i].hgmdVarHPOIdList=varDict['hgmdVarHPOIdList']
            varObjList[i].hgmdVarHPOStrList=varDict['hgmdVarHPOStrList']

            varObjList[i].hgmd_id = varDict['hgmd_id'] # CL added
            varObjList[i].hgmd_symbol = varDict['hgmd_symbol'] # CL added
            varObjList[i].hgmd_PHEN = varDict['hgmd_PHEN'] # CL added
            varObjList[i].hgmd_CLASS = varDict['hgmd_CLASS'] # CL added

            #symptom
            varObjList[i].SymptomMatched=varDict['SymptomMatched']
            varObjList[i].symptomScore=varDict['symptomScore']
            varObjList[i].symptomName=varDict['symptomName']
            varObjList[i].omimSymptomSimScore=varDict['omimSymptomSimScore']
            varObjList[i].omimSymMatchFlag=varDict['omimSymMatchFlag']
            varObjList[i].hgmdSymptomScore=varDict['hgmdSymptomScore']
            varObjList[i].hgmdSymptomSimScore=varDict['hgmdSymptomSimScore']
            varObjList[i].hgmdSymMatchFlag=varDict['hgmdSymMatchFlag']
            varObjList[i].clinVarSymMatchFlag=varDict['clinVarSymMatchFlag']
            #increment
            i=i+1
            #print('i index:', i)
    #return
    return varObjList

def getAnnotateInfo_2(rowList,genomeRef,clinvarGeneDf,clinvarAlleleDf, omimGeneList, omimAlleleList,hgmdDf, moduleList, dgvDf, decipherDf, gnomadMetricsGeneDf):
    '''
    Function to read the input annoatated varaint file and call the marrvel flat files to get annotations such as OMIM, DGV, Decipher
    Parameters:
    rowList: the input variants as a lits
    other inputs are the flat file objects read locally instaed of an API
    
    Returns:
    A dict of annoatated varaints
    '''
    debugFlag=0
    #print('\nlen numVars:', len(rowList), 'typerow00:', type(rowList[0][0]), 'typerow[0]:',
    #      type(rowList[0]) )
    #print('len(row[0][0]):', len(rowList[0][0]))
    #print('moduleList:', moduleList)
    varObjList=[]
    varDictList=[]
    i=0
    #loop thru the rows
    for row in rowList:
        # CL 03-14-2023: commented all printing lines
        #print('type of row:', type(row))
        varObj=Variant()
        transcriptId=row.Feature
        #s=row.Uploaded_variation.split('_') '1_10204_-/T' 1_1588250_T_A
        ####row[0]: 21_11039079_C/A
        ####s: ['21', '11039079', 'C/A']
        #print('row[0]:', row[0])

        #two ways of input of first column either 1_1588250_T_A OR 21_11039079_C/A, so use the option flag
        optFlag=0
        if(row[0].find('/')!=-1):
            optFlag=1
        
        if optFlag==0:
            s=row[0].split('_')
            #print('s:', s)
            chrom=s[0]
            pos=int(s[1])
            ref=s[2]
            alt=s[3]
        elif optFlag==1:
            s=row[0].split('_')
            #print('s:', s)
            chrom=s[0]
            pos=int(s[1])
            s=s[2].split('/')
            ref=s[0]
            alt=s[1]
       
        #get the start and stop from second column like '1:10203-10204'
        if '-' in row[1]:
            s=row[1].split(':')
            tmp=s[1]
            s=tmp.split('-')
            #print('s:',s)
            start=int(s[0])
            stop=int(s[1])
        else:
            #start and stop the same
            s=row[1].split(':')
            start=int(s[1])
            stop=int(s[1])

        #print('chrom:', chrom,'pos:',pos,'ref:',ref,'alt:',alt,'start:',start,'stop:',stop)
        #change chrom X and Y and MT to numbers
        if chrom=='X':
            chrom=23
        elif chrom=='Y':
            chrom=24
        elif chrom=='MT':
            chrom=25
        elif re.search(r'GL', chrom):
            chrom=26
        chrom=int(chrom)

        #if it is hg38 get its hg19 coordinates
        # CL 03-14-2023: we have separate database for hg19 and hg38,
        #                we don't need to use LiftOver which is inaccurate
        #                related codes commented and modified
        if genomeRef=='hg38':
            varObj.hg38Chrom=chrom
            varObj.hg38Pos=pos
            varObj.chrom=chrom
            varObj.pos=pos
            varObj.start=start
            varObj.stop=stop

            '''
            retList=gethg19LocFromHg38(chrom, pos)#called from the utils_1.py
            # retList=[newChrom, newPos]
            varObj.hg19Chrom=retList[0]
            varObj.hg19Pos=retList[1]
            varObj.chrom=retList[0]
            varObj.pos=retList[1]
            #get the start
            retList=gethg19LocFromHg38(chrom, start)
            varObj.start=int(retList[1])
            #get the stop
            retList=gethg19LocFromHg38(chrom, stop)
            varObj.stop=int(retList[1])
            '''
        else:
            varObj.hg19Chrom=chrom
            varObj.hg19Pos=pos
            varObj.chrom=chrom
            varObj.pos=pos
            varObj.start=start
            varObj.stop=stop


        geneSymbol=row.SYMBOL
        #print('gene:', geneSymbol)
        varObj.geneSymbol=geneSymbol
        varObj.CADD_phred=row.CADD_phred
        varObj.CADD_PHRED=row.CADD_PHRED

        #assign
        varObj.ref=ref
        varObj.alt=alt
        varObj.varId_dash='-'.join([str(chrom),str(start),ref,alt])
        #print('varId dash:', varObj.varId_dash)
        varId='_'.join([str(chrom),str(pos),ref,alt,transcriptId])
        varObj.varId=varId
        if 'ZYG' in row:
            varObj.zyg=row.ZYG
        varObj.geneEnsId=row.Gene
        varObj.rsId=row.Existing_variation
        varObj.GERPpp_RS=row.GERPpp_RS
        varObj.featureType=row.Feature_type
        varObj.gnomadAF=row.gnomAD_AF
        varObj.gnomadAFg=row.gnomADg_AF
        varObj.CLIN_SIG=row.CLIN_SIG #CL: useless but kept for now
        varObj.LRT_Omega=row.LRT_Omega
        varObj.LRT_score=row.LRT_score
        varObj.phyloP100way_vertebrate=row.phyloP100way_vertebrate
        varObj.IMPACT=row.IMPACT
        varObj.Consequence=row.Consequence
        varObj.HGVSc=row.HGVSc
        varObj.HGVSp=row.HGVSp
        #dbnsfp attributes 
        varObj.GERPpp_NR=row.GERPpp_NR
        varObj.DANN_score=row.DANN_score
        varObj.FATHMM_pred=row.FATHMM_pred
        varObj.FATHMM_score=row.FATHMM_score
        varObj.GTEx_V8_gene=row.GTEx_V8_gene
        varObj.GTEx_V8_tissue=row.GTEx_V8_tissue
        varObj.Polyphen2_HDIV_score=row.Polyphen2_HDIV_score
        varObj.Polyphen2_HVAR_score=row.Polyphen2_HVAR_score
        varObj.REVEL_score=row.REVEL_score
        varObj.SIFT_score=row.SIFT_score

        varObj.clinvar_AlleleID = row.clinvar # Clinvar allele ID from clinvar.vcf.gz
        varObj.clinvar_clnsig= row.clinvar_CLNSIG #CL: Clinvar SIG from clinvar.vcf.gz
        #varObj.clinvar_clnsig = row.clinvar_clnsig #CL: Clinvar SIG from VEP, deleted
        varObj.clinvar_CLNREVSTAT= row.clinvar_CLNREVSTAT #CL: Clinvar STAT from clinvar.vcf.gz, for interface only
        varObj.clinvar_CLNSIGCONF= row.clinvar_CLNSIGCONF #CL: Clinvar SIGCONF from clinvar.vcf.gz
        varObj.clin_code = row.clinvar_CLNSIG #CL: feature name for ai

        varObj.fathmm_MKL_coding_score=row.fathmm_MKL_coding_score
        varObj.LRT_score=row.LRT_score
        varObj.LRT_Omega=row.LRT_Omega
        varObj.phyloP100way_vertebrate=row.phyloP100way_vertebrate
        varObj.M_CAP_score=row.M_CAP_score
        varObj.MutationAssessor_score=row.MutationAssessor_score
        varObj.MutationTaster_score=row.MutationTaster_score
        varObj.ESP6500_AA_AC=row.ESP6500_AA_AC
        varObj.ESP6500_AA_AF=row.ESP6500_AA_AF
        varObj.ESP6500_EA_AC=row.ESP6500_EA_AC
        varObj.ESP6500_EA_AF=row.ESP6500_EA_AF

        varObj.VARIANT_CLASS=row.VARIANT_CLASS
        varObj.Feature=row.Feature
        varObj.hom=row.gnomADg_controls_nhomalt
        varObj.hgmd_id=row.hgmd # CL added
        varObj.hgmd_symbol=row.hgmd_GENE # CL added
        varObj.hgmd_rs=row.hgmd_RANKSCORE
        varObj.hgmd_PHEN=row.hgmd_PHEN # CL added
        varObj.hgmd_CLASS=row.hgmd_CLASS # CL added

        if row.clinvar_CLNSIGCONF !="-":
            clin_dict=dict()
            for ro in row.clinvar_CLNSIGCONF.split("|_"):
                temp=ro.split("(")
                clin_dict[temp[0]]=int(temp[1][0])
            PLP_sum=clin_dict.get("Pathogenic", 0)+clin_dict.get("Likely_pathogenic", 0)
            varObj.clin_dict=clin_dict
            varObj.clin_PLP=PLP_sum
            varObj.clin_PLP_perc=PLP_sum/sum(clin_dict.values())
        else:
            if "benign" in row.clinvar_clnsig.lower():
                varObj.clin_PLP_perc=0
            elif "pathogenic" in row.clinvar_clnsig.lower():
                varObj.clin_PLP_perc=1
            else:
                varObj.clin_PLP_perc="-"
            varObj.clin_PLP="-"
            varObj.clin_dict="-"


        if row.SpliceAI_pred !="-":
            varObj.spliceAI=row.SpliceAI_pred
            temp=row.SpliceAI_pred.split("|")
            varObj.spliceAImax=max(float(temp[1]),float(temp[2]),float(temp[3]),float(temp[4]))
        else:
            varObj.spliceAI="-"
            varObj.spliceAImax="-"

        #get dgv
        if 'conserve' in moduleList:
            #print('\nGetting DGV')
            dgvRet=getDGVUsingMarrvelFlatFile(varObj, dgvDf)
            varObj.dgvDictList=dgvRet[0]
            varObj.dgvTypeList=dgvRet[1]
            varObj.dgvSubtypeList=dgvRet[2]
            varObj.dgvVarFound=dgvRet[3]

        #get decipher
        if 'conserve' in moduleList:
            #print('\nGetting DECIPHER')
            decipherRet=getDecipherUsingMarrvelFlatFile(varObj, decipherDf)
            #[decipherDictList,decipherDeletionObsList,decipherStudyList, decipherVarFound]
            varObj.decipherDictList=decipherRet[0]
            varObj.decipherDeletionObsList=decipherRet[1]
            varObj.decipherStudyList=decipherRet[2]
            varObj.decipherVarFound=decipherRet[3]

        #get gnomad gene metrics from gnomad file
        if 'conserve' in moduleList:
            #print('\nGetting gnomad gene metrics')
            gnomadGeneMetricRet=getGnomadGeneMetricFlatFile(varObj, gnomadMetricsGeneDf)
            #[decipherDictList,decipherDeletionObsList,decipherStudyList, decipherVarFound]
            varObj.gnomadGeneZscore=gnomadGeneMetricRet[0]
            varObj.gnomadGenePLI=gnomadGeneMetricRet[1]
            varObj.gnomadGeneOELof=gnomadGeneMetricRet[2]#O/E lof
            varObj.gnomadGeneOELofUpper=gnomadGeneMetricRet[3]#O/E lof upper

        #for module one will need OMIM, HGMD, clinvar
        #get OMIM
        if 'curate' in moduleList:
            #print('\nGetting OMIM')
            #varObj.omimList=jsonDict['omim']
            #retList=[varFound, geneFound, omimDict, omimGeneDict, omimAlleleDict]
            omimRet=getOMIMUsingMarrvelFlatFile(varObj, omimGeneList, omimAlleleList)
            varObj.omimVarFound=omimRet[0]
            varObj.omimGeneFound=omimRet[1]
            varObj.omimDict=omimRet[2]
            varObj.omimGeneDict=omimRet[3]
            varObj.omimAlleleDict=omimRet[4]
            varObj.phenoList=omimRet[5]
            varObj.phenoInhList=omimRet[6]
            varObj.phenoMimList=omimRet[7]
            #print('OMIM res:')
            #print('\tgeneFound:',varObj.omimGeneFound,'varFound:',varObj.omimVarFound )

        #get clinvar
        if 'curate' in moduleList:
            #print('\nReading clinVar')
            clinVarRet=getClinVarUsingMarrvelFlatFile(varObj, clinvarAlleleDf, clinvarGeneDf)
            varObj.clinVarVarFound=clinVarRet[0]
            varObj.clinVarVarDict=clinVarRet[1]
            varObj.clinVarGeneFound=clinVarRet[2]
            varObj.clinVarGeneDict=clinVarRet[3]
            varObj.clinvarTotalNumVars=clinVarRet[4]
            varObj.clinvarNumP=clinVarRet[5]
            varObj.clinvarNumLP=clinVarRet[6]
            varObj.clinvarNumLB=clinVarRet[7]
            varObj.clinvarNumB=clinVarRet[8]
            varObj.clinvarTitle= clinVarRet[9]
            varObj.clinvarSignDesc= row.clinvar_CLNSIG #clinVarRet[10] #CL: changed to clinvar.vcf.gz annotation
            varObj.clinvarCondition=clinVarRet[11]
            #print('clinVar res:')
            '''
            if debugFlag==1:
                print('\tgeneFound::',varObj.clinVarGeneFound,'varFound:',varObj.clinVarVarFound)
                print('\tnumVars:',varObj.clinvarTotalNumVars,'numPathologic:',varObj.clinvarNumP,'numBenign:',varObj.clinvarNumB)
                print('\tsignDesc:', varObj.clinvarSignDesc)
            '''

        #get HGMD
        if 'curate' in moduleList:
            #print('\nReading HGMD')
            hgmdRet=getHGMDUsingFlatFile(varObj,hgmdDf)
            #hgmdVarFound,hgmdGeneFound,hgmdVarPhenIdList,hgmdVarHPOIdList,hgmdVarHPOStrList
            varObj.hgmdVarFound=hgmdRet[0]
            varObj.hgmdGeneFound=hgmdRet[1]
            varObj.hgmdVarPhenIdList=hgmdRet[2]
            varObj.hgmdVarHPOIdList=hgmdRet[3]
            varObj.hgmdVarHPOStrList=hgmdRet[4]
            #print('HGMD results:')
            #print('\thgmdVarFound:',varObj.hgmdVarFound,'hgmdGeneFound:',varObj.hgmdGeneFound,
            #      'hgmdVarPhenIdList:',varObj.hgmdVarPhenIdList,'hgmdVarHPOIdList:',
            #      varObj.hgmdVarHPOIdList,
            #      'hgmdVarHPOStrList:',varObj.hgmdVarHPOStrList)

        #make a dict from the varObj. the varObj has one element for each row of the input VEP annotated file
        varDict={
            'hg19Chrom':varObj.hg19Chrom,
            'hg19Pos':varObj.hg19Pos,
            'chrom':varObj.chrom,
            'pos':varObj.pos,
            'start':varObj.start,
            'stop':varObj.stop,
            'geneSymbol':varObj.geneSymbol,
            'CADD_phred':varObj.CADD_phred,
            'CADD_PHRED':varObj.CADD_PHRED,
            'ref':varObj.ref,
            'alt':varObj.alt,
            'varId':varObj.varId,
            'ZYG':varObj.zyg,
            'HGVSc':varObj.HGVSc,
            'HGVSp':varObj.HGVSp,
            'Gene':varObj.geneEnsId,
            'Existing_variation':varObj.rsId,
            'GERPpp_RS':varObj.GERPpp_RS,
            'Feature_type':varObj.featureType,
            'gnomadAF':varObj.gnomadAF,
            'gnomadAFg':varObj.gnomadAFg,
            'CLIN_SIG':varObj.CLIN_SIG,
            'LRT_Omega':varObj.LRT_Omega,
            'LRT_score':varObj.LRT_score,
            'phyloP100way_vertebrate':varObj.phyloP100way_vertebrate,
            #dbnsfp attributes 
            'GERPpp_NR':varObj.GERPpp_NR,
            'DANN_score':varObj.DANN_score,
            'FATHMM_pred':varObj.FATHMM_pred,
            'FATHMM_score':varObj.FATHMM_score,
            'GTEx_V8_gene':varObj.GTEx_V8_gene,
            'GTEx_V8_tissue':varObj.GTEx_V8_tissue,
            'Polyphen2_HDIV_score':varObj.Polyphen2_HDIV_score,
            'Polyphen2_HVAR_score':varObj.Polyphen2_HVAR_score,
            'REVEL_score':varObj.REVEL_score,
            'SIFT_score':varObj.SIFT_score,

            'clinvar_AlleleID': varObj.clinvar_AlleleID, # Clinvar allele ID from clinvar.vcf.gz
            'clinvar_clnsig': varObj.clinvar_clnsig, #CL: Clinvar SIG from clinvar.vcf.gz
            'clinvar_CLNREVSTAT': varObj.clinvar_CLNREVSTAT, #CL: Clinvar STAT from clinvar.vcf.gz, for interface only
            'clinvar_CLNSIGCONF': varObj.clinvar_CLNSIGCONF, #CL: Clinvar SIGCONF from clinvar.vcf.gz
            'clin_code': varObj.clin_code, #CL: feature for ai

            'fathmm_MKL_coding_score':varObj.fathmm_MKL_coding_score,
            'LRT_score':varObj.LRT_score,
            'LRT_Omega':varObj.LRT_Omega,
            'phyloP100way_vertebrate':varObj.phyloP100way_vertebrate,
            'M_CAP_score':varObj.M_CAP_score,
            'MutationAssessor_score':varObj.MutationAssessor_score,
            'MutationTaster_score':varObj.MutationTaster_score,
            'ESP6500_AA_AC':varObj.ESP6500_AA_AC,
            'ESP6500_AA_AF':varObj.ESP6500_AA_AF,
            'ESP6500_EA_AC':varObj.ESP6500_EA_AC,
            'ESP6500_EA_AF':varObj.ESP6500_EA_AF,
            #dbnsfp
            'gnomadGeneZscore':varObj.gnomadGeneZscore,
            'gnomadGenePLI':varObj.gnomadGenePLI,
            'gnomadGeneOELof':varObj.gnomadGeneOELof,#O/E lof
            'gnomadGeneOELofUpper':varObj.gnomadGeneOELofUpper,#O/E lof upper,
            'IMPACT':varObj.IMPACT,
            'Consequence':varObj.Consequence,
            'omimVarFound':varObj.omimVarFound,
            'omimGeneFound':varObj.omimGeneFound,
            'omimDict':varObj.omimDict,
            'omimGeneDict':varObj.omimGeneDict,
            'omimAlleleDict':varObj.omimAlleleDict,
            'phenoList':varObj.phenoList,
            'phenoInhList':varObj.phenoInhList,
            'phenoMimList':varObj.phenoMimList,
            'clinVarVarFound':varObj.clinVarVarFound,
            'clinVarVarDict':varObj.clinVarVarDict,
            'clinVarGeneFound':varObj.clinVarGeneFound,
            'clinVarGeneDict':varObj.clinVarGeneDict,
            'clinvarTotalNumVars':varObj.clinvarTotalNumVars,
            'clinvarNumP':varObj.clinvarNumP,
            'clinvarNumLP':varObj.clinvarNumLP,
            'clinvarNumLB':varObj.clinvarNumLB,
            'clinvarNumB':varObj.clinvarNumB,
            'clinvarTitle':varObj.clinvarTitle,
            'clinvarSignDesc':varObj.clinvarSignDesc,
            'clinvarCondition':varObj.clinvarCondition,
            'hgmdVarFound':varObj.hgmdVarFound,
            'hgmdGeneFound':varObj.hgmdGeneFound,
            'hgmdVarPhenIdList':varObj.hgmdVarPhenIdList,
            'hgmdVarHPOIdList':varObj.hgmdVarHPOIdList,
            'hgmdVarHPOStrList':varObj.hgmdVarHPOStrList,
            'varId_dash':varObj.varId_dash,
            'dgvDictList':varObj.dgvDictList,
            'dgvTypeList':varObj.dgvTypeList,
            'dgvSubtypeList':varObj.dgvSubtypeList,
            'dgvVarFound':varObj.dgvVarFound,
            'decipherDictList':varObj.decipherDictList,
            'decipherDeletionObsList':varObj.decipherDeletionObsList,
            'decipherStudyList':varObj.decipherStudyList,
            'decipherVarFound':varObj.decipherVarFound,
            'gnomadGeneZscore':varObj.gnomadGeneZscore,
            'gnomadGenePLI':varObj.gnomadGenePLI,
            'gnomadGeneOELof':varObj.gnomadGeneOELof,
            'gnomadGeneOELofUpper':varObj.gnomadGeneOELofUpper,
            #symptom
            'SymptomMatched':varObj.SymptomMatched,
            'symptomScore':varObj.symptomScore,
            'symptomName':varObj.symptomName,
            'omimSymptomSimScore':varObj.omimSymptomSimScore,
            'omimSymMatchFlag':varObj.omimSymMatchFlag,
            'hgmdSymptomScore':varObj.hgmdSymptomScore,
            'hgmdSymptomSimScore':varObj.hgmdSymptomSimScore,
            'hgmdSymMatchFlag':varObj.hgmdSymMatchFlag,
            'clinVarSymMatchFlag':varObj.clinVarSymMatchFlag,
            'VARIANT_CLASS':varObj.VARIANT_CLASS,
            'Feature':varObj.Feature,
            'hom':varObj.hom,
            'hgmd_rs':varObj.hgmd_rs,

            'hgmd_id':varObj.hgmd_id, # CL added
            'hgmd_symbol':varObj.hgmd_symbol, # CL added
            'hgmd_PHEN':varObj.hgmd_PHEN, # CL added
            'hgmd_CLASS':varObj.hgmd_CLASS, # CL added

            'clin_dict':varObj.clin_dict,
            'clin_PLP':varObj.clin_PLP,
            'clin_PLP_perc':varObj.clin_PLP_perc,
            'spliceAI':varObj.spliceAI,
            'spliceAImax':varObj.spliceAImax

        }
        varDictList.append(varDict)
        #print('finsihed making dict')
        #increment index
        i=i+1
    
    #return object
    return varDictList

