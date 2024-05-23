## UPDATE
#    Changes made by Chaozhong noted as #CL
#        1. added clinvar feature and changed database to clinvar.vcf.gz in VEP annotation
#        2. chang the HGMD database file to VEP annotation and similarity screo files
# hg19
import time
import math
import os
import argparse
import re
import pandas as pd
import numpy as np
import json
import sys
import requests
import subprocess
#local methods
from utils_1 import  *
#from utils_for_marrvel_annotate_script import *
from utils_for_symMatch import *
from utils_for_marrvel_flatfile import *
#from utils_for_marrvel_flatfile_module_3 import *
import multiprocessing
from functools import partial
from marrvel_score_recalc import *


#############################################################
# Database list:
#   (common for all)  gnomad.v2.1.1.lof_metrics.by_gene.txt
#   (hg38 pos needed) gene_clinvar.csv
#   (common for all)  gene_omim.json
#   (common for all)  omim_alleric_variants.json
#   (hg38 ver needed) dgv.csv (/houston_30t/chaozhong/aimarrvel_pipeline/database/dgv/dgv.38.csv)
#   (hg38 pos needed) decipher.csv
#
#############################################################




def checkUserArgs(inArgs):
    #check the user arguments and if anything necessary is missing
    pass


def main():
    #example of running the script:
    #python3 main.py -outPrefix ./test/r1 \
    #-patientID UDN630665 -patientHPOsimiOMIM  /houston_10t/dongxue/MARRVEL_AI/BG/HPO2Gene/out_May12021/UDN/HPOsimi_UDN630665_simi_0.tsv \
    #-patientHPOsimiHGMD  /six_tera/chaozhong/module_1/out/UDN/HPOsimi_UDNUDN630665_simi_0.tsv \
    #-varFile /houston_20t/rami/vep_PassVarQual_filt_UDN_Apr2521/540200-UDN630665-P_HGYTFBCXY-2-ID01-vep_gnomADg_gnomADe_splice_clin_hgmd_filt.txt \
    #-inFileType vepAnnotTab -patientFile testPatientFile.txt -patientFileType one -genomeRef hg19 -diseaseInh AD -modules curate,conserve
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-varFile", "--varFile",  help="Proivde input variant file")
    parser.add_argument("-inFileType", "--inFileType",  help="Proivde type of input file:vcf, vepAnnotTab")
    parser.add_argument("-patientFile", "--patientFile",  help="Proivde HPO IDs")#currently not used
    parser.add_argument("-patientFileType", "--patientFileType",  help="Proivde type of file:one, two")#currently not used
    parser.add_argument("-patientID", "--patientID",  help="Proivde patientID")
    parser.add_argument("-patientHPOsimiOMIM", "--patientHPOsimiOMIM",  help="Proivde patient HPO similarity file-OMIM")
    parser.add_argument("-patientHPOsimiHGMD", "--patientHPOsimiHGMD",  help="Proivde patient HPO similarity file-HGMD")
    parser.add_argument("-diseaseInh", "--diseaseInh",  help="Proivde disease Inheritance:AD, AR, XD, XR")
    parser.add_argument("-modules", "--modules",  help="Select modules to run:all,curate,conserve,effectOnGene")
    parser.add_argument("-genomeRef", "--genomeRef",  help="Proivde genome ref: hg19, hg38")
    parser.add_argument("-outPrefix", "--outPrefix",  help="Proivde prefix for output files")
    args = parser.parse_args()
    #check the user args
    checkUserArgs(args)
    print('input file:', args.varFile)
    print('type of input file:', args.inFileType)
    print('outPrefix:', args.outPrefix)
    print('modules:', args.modules)
    moduleList=args.modules.split(',')
    print('modules list:', moduleList)
    
    #start time
    start_time = time.time()

    #need to write user options to a log file

    #read the gnomad gene metrics file
    fileName='/run/data_dependencies/annotate/anno_hg19/gnomad.v2.1.1.lof_metrics.by_gene.txt'
    #fileName='/database_test/gnomad.v2.1.1.lof_metrics.by_gene.txt'
    gnomadMetricsGeneDf=pd.read_csv(fileName, sep='\t')
    #print('gnomadMetricsGeneDf shape:', gnomadMetricsGeneDf.shape)

    #set up flags
    #debug flag; if 1 then can check/test the code
    debugFlag=0
    #print('debugFlag:', debugFlag)
    #use flat files to read the marrvel databases instead of using API
    flatFileFlag=1
    #use marrvel annotation script and API to get marrvel info. this is slower
    annotateFlag=0

    #the directory for marrvel database flat file
    #flatFilesDir='/home/dongxue/hdd_10t/MARRVEL_AI/resource/MARRVEL/'

    #initialization
    dgvDf=[]
    decipherDf=[]
    hgmdDf=[]
    omimHPOScoreDf=[]
    hgmdHPOScoreDf=[]
    clinvarGeneDf=[]
    clinvarAlleleDf=[]
    omimAlleleList=[]


    ## CL: It's no longer needed to extract info from HGMD_allmut.txt
    """
    #read HGMD files
    if 'curate' in moduleList:
        #hgmd dir
        fileName='/run/annotation/HGMD_allmut.txt'
        hgmdDf=pd.read_csv(fileName, sep='\t')
        print('hgmdDf shape:', hgmdDf.shape)
        hgmdDf=hgmdDf.fillna(0)
        hgmdDf['startCoord']=hgmdDf['startCoord'].astype(int)
        hgmdDf['endCoord']=hgmdDf['endCoord'].astype(int)
        hgmdDf['chromosome']=hgmdDf['chromosome'].replace('X',23  )
        hgmdDf['chromosome']=hgmdDf['chromosome'].replace('Y',24  )
        hgmdDf['chromosome']=hgmdDf['chromosome'].replace('MT',25  )
        hgmdDf['chromosome']=hgmdDf['chromosome'].replace('GL.*','26',regex=True )
        #make chr as int
        hgmdDf['chromosome']=hgmdDf['chromosome'].astype(int)
        #sort by gene name
        hgmdDf.sort_values("gene", inplace=True)
        #setting gene as index column
        hgmdDf.set_index("gene", inplace = True,drop=False)

        #test/check HGMD
        if debugFlag==1:
            print('in HGMD')
            #print('\tvar:', varObj.varId_dash, 'var-gene:', varObj.geneSymbol)
            geneSymbol='SAMD11'
            pos=865595
            chrom='1'
            print('\tchrom:', chrom, 'pos:',pos,'geneSymbol:', geneSymbol)
            HGMDDict={}
            hgmdGeneFound=0
            hgmdVarFound=0
            hgmdVarPhenIdList=[]
            hgmdVarHPOIdList=[]
            hgmdVarHPOStrList=[]
            #check variant
            print('\ttype pos:', type(hgmdVarDf['pos'][0]), 'type chrom:', type(hgmdVarDf['chrom'][0]))
            vals=hgmdVarDf[(hgmdVarDf['pos']==pos) & (hgmdVarDf['chrom'] == int(chrom))]
            numRows=len(vals.index)
            print('\tvar numRows:', numRows)
            if numRows>0:
                hgmdVarFound=1
                hgmdVarPhenIdList.extend(vals['phen_id'].tolist())
                hgmdVarHPOIdList.extend(vals['hpo_id'].tolist())
                hgmdVarHPOStrList.extend(vals['hpo_str'].tolist())

            #check gene
            vals=hgmdGeneDf[(hgmdGeneDf['gene']==geneSymbol)]
            numRows=len(vals.index)
            print('\tgene numRows:', numRows)
            if numRows>0:
                hgmdGeneFound=1
            print('\thgmdVarFound:',hgmdVarFound,'hgmdGeneFound:',hgmdGeneFound,
                  'hgmdVarPhenIdList:',hgmdVarPhenIdList,'hgmdVarHPOIdList:',hgmdVarHPOIdList,
                  'hgmdVarHPOStrList:',hgmdVarHPOStrList)

    """

    #read HPO files
    if 'curate' in moduleList:
        fileName=args.patientHPOsimiOMIM
        omimHPOScoreDf = pd.read_csv(fileName, sep='\t')
        print('patientHPOsimi-OMIM dimension:',omimHPOScoreDf.shape)

        fileName=args.patientHPOsimiHGMD
        hgmdHPOScoreDf = pd.read_csv(fileName, sep='\t')
        print('patientHPOsimi-HGMD dimension:',hgmdHPOScoreDf.shape)

    #read the clinvar genes file
    if 'curate' in moduleList:

        # Q: Why we use the same for hg19 and hg38?
        # A: 
        if args.genomeRef == 'hg38':
            fileName='/run/data_dependencies/annotate/anno_hg19/gene_clinvar.csv'
        else:
            fileName='/run/data_dependencies/annotate/anno_hg19/gene_clinvar.csv'
        #fileName='/database_test/gene_clinvar.csv'
        #print('\nfileName:', fileName)
        clinvarGeneDf=pd.read_csv(fileName, sep=',')
        #sort by gene name
        clinvarGeneDf.sort_values("symbol", inplace=True)
        clinvarGeneDf.set_index(['symbol'], inplace = True,drop=False)
        #print('clinvarGenes dim:', clinvarGeneDf.shape)#41754 x 10,, entrezId,symbol,chr,hg19Start,hg19Stop,totalClinvarVars,P,LP,LB,B
        #print(clinvarGeneDf.head())
        
        ### CL: this part was substituted by VEP custom annotation from clinvar.vcf.gz
        """
        #read the clinvar variant file
        fileName='/run/annotation/clinvar.csv'
        print('\nfileName:', fileName)
        clinvarAlleleDf=pd.read_csv(fileName, sep=',')
        print('clinvarAllele dim:', clinvarAlleleDf.shape)# (507196, 5)  title chr start stop significance.description

        #convert float to int
        clinvarAlleleDf = clinvarAlleleDf.fillna(0)
        clinvarAlleleDf['start'] = clinvarAlleleDf['start'].astype(int)
        clinvarAlleleDf['stop'] = clinvarAlleleDf['stop'].astype(int)
        #in chr column change X and Y to numbers
        clinvarAlleleDf['chr']=clinvarAlleleDf['chr'].replace('X',23)
        clinvarAlleleDf['chr']=clinvarAlleleDf['chr'].replace('Y',24)
        clinvarAlleleDf['chr']=clinvarAlleleDf['chr'].replace('MT',25)
        clinvarAlleleDf['chr']=clinvarAlleleDf['chr'].replace('GL.*','26',regex=True )
        #make chr as int
        clinvarAlleleDf['chr']=clinvarAlleleDf['chr'].astype(int)
        #make index
        clinvarAlleleDf['id1']=clinvarAlleleDf['chr'].astype(str)+ '_' +clinvarAlleleDf['start'].astype(str)+'_'+ clinvarAlleleDf['stop'].astype(str)
        # setting first name as index column
        clinvarAlleleDf.set_index("id1", inplace = True,drop=False)
        #make sure start and stop are int
        print('type chr:', type(clinvarAlleleDf['chr'][0]), 'start:', type(clinvarAlleleDf['start'][0]), 'stop:', type(clinvarAlleleDf['stop'][0]))
        """

    #read OMIM
    if 'curate' in moduleList:
        #read the OMIM gene file
        fileName='/run/data_dependencies/annotate/anno_hg19/gene_omim.json'
        #fileName='/database_test/gene_omim.json'
        #print('\nfileName:', fileName)
        with open(fileName) as f:
            omimGeneList=json.load(f)
        #print('type omimGeneList:', type(omimGeneList), 'len:', len(omimGeneList))
        #check the type of data and keys
        if debugFlag==1:
            for omimGeneDict in omimGeneList:
                print('type of omimGeneDict:', type(omimGeneDict))
                print('keys:', omimGeneDict.keys())
                for keyVal in omimGeneDict.keys():
                    print('keyVal:', keyVal)
                    print('\tsubkeys type:', type(omimGeneDict[keyVal]) )
                    if  isinstance(omimGeneDict[keyVal], list):
                        print('\n\t\tfound list')
                        print('\t\t type:',type(omimGeneDict[keyVal]))
                        #print('\t\t type:',omimGeneDict[keyVal])
                        #for row in omimGeneDict[keyVal]:
                        #    print('row keys:', row.keys())

                break

        #read the OMIM allele file
        fileName='/run/data_dependencies/annotate/anno_hg19/omim_alleric_variants.json'
        #fileName='/database_test/omim_alleric_variants.json'
        #print('\nfileName:', fileName)
        with open(fileName) as f:
            omimAlleleList=json.load(f)
        if debugFlag==1:
            print('type of omim allelic :', type(omimAlleleList), 'len:', len(omimAlleleList))
            for omimAlleleDict in omimAlleleList:
                print('type of omimAlleleDict:', type(omimAlleleDict))
                #check the keys
                print('keys:', omimAlleleDict.keys())
                for keyVal in omimAlleleDict.keys():
                    print('keyVal:', keyVal)
                    print('\tsubkeys type:', type(omimAlleleDict[keyVal]) )
                break
                
    #read DGV database
    if 'conserve' in moduleList:
        print('reading DGV flat file')
        if args.genomeRef == 'hg38':
            fileName='/run/data_dependencies/annotate/anno_hg38/dgv.csv'
        else:
            fileName='/run/data_dependencies/annotate/anno_hg19/dgv.csv'
        #fileName='/database_test/dgv.csv'
        dgvDf=pd.read_csv(fileName, sep=',')
        #print('dgvDf shape:', dgvDf.shape)
        #print('dgvHead:', dgvDf.head())
        #hg19Chr	hg19Start	hg19Stop
        dgvDf=dgvDf.fillna(0)

        #if args.genomeRef == 'hg19':
        dgvDf.columns = ['Chr', 'Start', 'Stop'] + dgvDf.columns.tolist()[3:]

        dgvDf['Start']=dgvDf['Start'].astype(int)
        dgvDf['Stop']=dgvDf['Stop'].astype(int)
        dgvDf['Chr']=dgvDf['Chr'].replace('X',23  )
        dgvDf['Chr']=dgvDf['Chr'].replace('Y',24  )
        dgvDf['Chr']=dgvDf['Chr'].replace('MT',25  )
        dgvDf['Chr']=dgvDf['Chr'].replace('GL.*',26,regex=True )

        #make chr as int
        #dgvDf = dgvDf.loc[dgvDf['Chr'].isin(list(range(27))),:].copy()
        dgvDf['Chr']=dgvDf['Chr'].astype(int)

        ''' CL: change column name to be compatible with hg38
        dgvDf['hg19Start']=dgvDf['hg19Start'].astype(int)
        dgvDf['hg19Stop']=dgvDf['hg19Stop'].astype(int)
        dgvDf['hg19Chr']=dgvDf['hg19Chr'].replace('X',23  )
        dgvDf['hg19Chr']=dgvDf['hg19Chr'].replace('Y',24  )
        dgvDf['hg19Chr']=dgvDf['hg19Chr'].replace('MT',25  )
        dgvDf['hg19Chr']=dgvDf['hg19Chr'].replace('GL.*','26',regex=True )
        #make chr as int
        dgvDf['hg19Chr']=dgvDf['hg19Chr'].astype(int)
        '''
        print('finsihed reading DGV')

    #read DECIPHER database
    if 'conserve' in moduleList:
        print('reading Decipher flat file')
        if args.genomeRef == 'hg38':
            fileName='/run/data_dependencies/annotate/anno_hg38/decipher.csv'
        else:
            fileName='/run/data_dependencies/annotate/anno_hg19/decipher.csv'
        #fileName='/database_test/decipher.csv'
        decipherDf=pd.read_csv(fileName, sep=',')
        #print('decipherDf shape:', decipherDf.shape)
        #hg19Chr,hg19Start,hg19Stop
        decipherDf=decipherDf.fillna(0)

        #if args.genomeRef == 'hg19':
        decipherDf.columns = ['Chr', 'Start', 'Stop'] + decipherDf.columns.tolist()[3:]

        decipherDf['Start']=decipherDf['Start'].astype(int)
        decipherDf['Stop']=decipherDf['Stop'].astype(int)
        decipherDf['Chr']=decipherDf['Chr'].replace('X',23  )
        decipherDf['Chr']=decipherDf['Chr'].replace('Y',24  )
        decipherDf['Chr']=decipherDf['Chr'].replace('MT',25  )
        decipherDf['Chr']=decipherDf['Chr'].replace('GL.*','26',regex=True )
        #make chr as int
        decipherDf['Chr']=decipherDf['Chr'].astype(int)

        ''' CL: change column name to be compatible with hg38
        decipherDf['hg19Start']=decipherDf['hg19Start'].astype(int)
        decipherDf['hg19Stop']=decipherDf['hg19Stop'].astype(int)
        decipherDf['hg19Chr']=decipherDf['hg19Chr'].replace('X',23  )
        decipherDf['hg19Chr']=decipherDf['hg19Chr'].replace('Y',24  )
        decipherDf['hg19Chr']=decipherDf['hg19Chr'].replace('MT',25  )
        decipherDf['hg19Chr']=decipherDf['hg19Chr'].replace('GL.*','26',regex=True )
        #make chr as int
        decipherDf['hg19Chr']=decipherDf['hg19Chr'].astype(int)
        '''
        print('finsihed reading DECIPHER')

    #read the input variant file
    if args.inFileType=='vepAnnotTab' and flatFileFlag==1:
        #numHeaderSkip=354#this depends on the type of annotation used
        numHeaderSkip=0
        with open(args.varFile,'r') as F:
            for line in F:
                if line.startswith('##'):
                    numHeaderSkip += 1
                else:
                    break
        print('input annoatated varFile:', args.varFile)
        t1 = time.time()
        varDf = pd.read_csv(args.varFile,sep='\t', skiprows=numHeaderSkip,error_bad_lines=False)
        #varDf=varDf[0:10]#do this if need to have a small test
        print('shape:', varDf.shape)
        t2 = time.time()
        inputReadTime = t2 - t1
        inputNumRows=len(varDf.index)

        #print('head:', varDf.head)
        #print('column names:', varDf.columns)
        #update the GERp column names since cant use ++ in the varaible name
        if 'GERP++_RS' in varDf.columns:
            print('found GERP++RS')
            varDf['GERPpp_RS'] = varDf['GERP++_RS']
        if 'GERP++_NR' in varDf.columns:
            print('found GERP++NR')
            varDf['GERPpp_NR'] = varDf['GERP++_NR']
        #update column names which have - in it
        if 'fathmm-MKL_coding_score' in varDf.columns:
            varDf['fathmm_MKL_coding_score']=varDf['fathmm-MKL_coding_score']
        if 'M-CAP_score' in varDf.columns:
            varDf['M_CAP_score']=varDf['M-CAP_score']
            

        print('start making variant objects')
        numCores=8
        p1 = time.time()
        #store rows into list for parallel function
        inputList =[]
        for index, rows in varDf.iterrows():
            inputList.append(rows)
        p2 = time.time()
        makeListTime = p2 - p1
        print('makeListTime :', makeListTime)
        print('type inputList:', type(inputList[0]), 'len inputList:', len(inputList))
        #group or chunks for parallel function
        groupSize=min((len(inputList) // numCores) + (len(inputList) %numCores), 1000)
        print('groupSize:', groupSize)
        groupList=list( getVarGroups(inputList, groupSize) )
        print('type groupList:', type(groupList), 'len:', len(groupList))

        #call parallel function, pass it the main function "getAnnotateInfo_2" from utils_for_marrvel_flatfile.py and any needed files. 
        #this function will get the marrvel annoatations from marrvel flat files
        parallelFlag=1
        if parallelFlag==1:
            pool = multiprocessing.Pool(processes=numCores)
            # CL: hgmdDf changed to hgmdHPOScoreDf
            func = partial(getAnnotateInfo_2, genomeRef=args.genomeRef,clinvarGeneDf=clinvarGeneDf, clinvarAlleleDf=clinvarAlleleDf,
                           omimGeneList=omimGeneList, omimAlleleList=omimAlleleList,hgmdDf=hgmdHPOScoreDf, moduleList=moduleList,
                           dgvDf=dgvDf, decipherDf=decipherDf, gnomadMetricsGeneDf=gnomadMetricsGeneDf)
            varDictList=pool.map(func,groupList)
            pool.close()
        else:
            pass

        print('finsihed getAnnotateInfo_2')
        #make a var obj list after paralle processing for each group, need a better more efficient way for doing this !!
        t1 = time.time()
        varObjList = getVarObjListFromDict(varDictList, len(varDf.index))
        t2 = time.time()
        makeVarObjListTime = t2 - t1
        #print('len varObjList:', len(varObjList))
        print('makeVarObjListTime 1:', makeVarObjListTime)

    #we now have the variant object list "varObjList" with annotations; loop thru it and get the scores
    for varObj in varObjList:
        #print('var:', varObj.varId_dash, 'var-gene:', varObj.geneSymbol, 'var-gene-pLI:', varObj.gnomadGenePLI)
        #the curate score is under the utils_1.py file
        if 'curate' in moduleList:
            #print('finding symptom matching')
            omimSymMatch(varObj, omimHPOScoreDf, args.inFileType)
            hgmdSymMatch(varObj, hgmdHPOScoreDf)
            clinVarSymMatch(varObj,args.inFileType)
            #OMIM and clinvar info
            retList=getCurationScore(varObj)#retList=[omimScore, hgmdScore, clinVarScore, totalScore]
            varObj.curationScoreOMIM=retList[0]
            varObj.curationScoreHGMD=retList[1]
            varObj.curationScoreClinVar=retList[2]
            varObj.curationScoreTotal=retList[3]
            #print('curation score total:', varObj.curationScoreTotal)
        #conserve score is under 
        if 'conserve' in moduleList:
            retList=getConservationScore(varObj, args.diseaseInh)
            #retList=[conservationScoreGnomad, conservationScoreDGV, conservationScoreOELof]
            varObj.conservationScoreGnomad=retList[0]
            varObj.conservationScoreDGV=retList[1]
            varObj.conservationScoreOELof=retList[2]
            #print('conservation score GERP++:',varObj.conservationScoreGerpPP,'DGV:',varObj.conservationScoreDGV,)
            #print('conserve score DGV:',varObj.conservationScoreDGV,'conserve score gnomad:', varObj.conservationScoreGnomad)


    #write the objects info scores to a file
    print('writing results')
    df = pd.DataFrame([t.__dict__ for t in varObjList ])
    print('shape of output file:',df.shape)
    #write
    fileName=args.outPrefix+'_'+args.patientID+'_scores.txt'
    print('out file name:', fileName)
    df.to_csv(fileName, sep='\t', index=False)

    end_time = time.time()
    process_time = end_time - start_time
    print('pipeline time:', process_time)
    #log file for times
    fileName=args.outPrefix+'_'+args.patientID+'_log.txt'
    f = open(fileName, "w")
    f.write("Process time:"+str(process_time)+ ' seconds and in mins:'+ str(process_time/60)+"\n")
    f.close()
    print('log file name:', fileName)
    print('input read time:', inputReadTime)
    print('input num rows:', inputNumRows)
    print('makeListTime from pd:', makeListTime)
    print('makeVarObjListTime 1:', makeVarObjListTime)
    print('m:', moduleList)

    print("Score re-calculation:")
    score = load_raw_matrix(df)
    score = omimCurate(score)
    score = hgmdCurate(score)
    score = clinvarCurate(score)
    score = conservationCurate(score)

    score.to_csv('/out/rami-test/%s_scores.csv'%(args.patientID), index=False)


    exit()



main()

