## UPDATE
#     Changes made by Chaozhong noted as #CL
#from liftover import get_lifter 
import pandas as pd
import numpy as np
import pprint
import json
import requests

#class for variant
class Variant():
    """Create a class to create an object for each variant
    Parameters:
        A class name for initiation call
    Output:
        An object of class Variant
    """
    def __init__(self):
        self.chrom='-'
        self.pos='-'
        self.start='-'
        self.stop='-'
        self.ref='-'
        self.alt='-'
        self.hg19Chrom='-'
        self.hg19Pos='-'
        self.hg38Chrom='-'
        self.hg38Pos='-'
        self.varId='-'
        self.varId_dash='-'#varinat ID separated by dash like '6-99365567-T-C'

        self.patientID = '-'
        self.zyg='-'
        self.geneSymbol='-'
        self.geneEnsId='-'
        self.rsId='-'
        self.featureType='-'
        self.gnomadAF='-'
        self.gnomadAFg='-'#genome gnomad score
        #dbnsfp attributes 
        self.CLIN_SIG='-'
        self.CADD_phred='-'#phred score from dbnsfp
        self.CADD_PHRED='-'
        self.GERPpp_RS='-'#GERP plus plus
        self.GERPpp_NR='-'
        self.DANN_score='-'
        self.FATHMM_pred='-'
        self.FATHMM_score='-'
        self.GTEx_V8_gene='-'
        self.GTEx_V8_tissue='-'
        self.Polyphen2_HDIV_score='-'
        self.Polyphen2_HVAR_score='-'
        self.REVEL_score='-'
        self.SIFT_score='-'
        
        self.clinvar_AlleleID = '-' #CL added
        self.clinvar_clnsig='-'
        self.clinvar_CLNREVSTAT =  '-' #CL added
        self.clinvar_CLNSIGCONF = '-' #CL added

        self.fathmm_MKL_coding_score='-'
        self.HGVSc='-'
        self.HGVSp='-' 
        #LRT
        self.LRT_score='-'
        self.LRT_Omega='-'
        #Phylo
        self.phyloP100way_vertebrate='-'
        self.M_CAP_score='-'
        self.MutationAssessor_score='-'
        self.MutationTaster_score='-'
        self.ESP6500_AA_AC='-'
        self.ESP6500_AA_AF='-'
        self.ESP6500_EA_AC='-'
        self.ESP6500_EA_AF='-'

        
        #symtom info
        self.SymptomMatched = {'omim':0, 'clinvar':0, 'hgmd':0} #list containing omimMatched, clinvarMatched, hgmdMatched
        self.symptomScore = {}#for OMIM
        self.symptomName = {}#for OMIM
        self.omimSymptomSimScore='-'
        self.omimSymMatchFlag= '-'
        self.hgmdSymptomScore='-'
        self.hgmdSymptomSimScore='-'
        self.hgmdSymMatchFlag='-'
        self.clinVarSymMatchFlag='-'
        
        #gnomad gene metrics from flat file
        self.gnomadGeneZscore='-'
        self.gnomadGenePLI='-'
        self.gnomadGeneOELof='-'#O/E lof
        self.gnomadGeneOELofUpper='-'#O/E lof upper
        #conseqeunce and impact
        self.IMPACT='-'
        self.Consequence='-'
        #omim
        self.omimDict={}
        self.omimGeneFound='-'
        self.omimVarFound='-'
        self.omimList=[]
        self.omimGeneDict={}
        self.omimAlleleDict={}
        self.phenoList=[]
        self.phenoInhList=[]
        self.phenoMimList=[]
        #HGMD
        self.HGMDDict={}
        self.hgmdGeneFound='-'
        self.hgmdVarFound='-'
        self.hgmdVarPhenIdList=[]
        self.hgmdVarHPOIdList=[]
        self.hgmdVarHPOStrList=[]
        #clinvar
        self.clinVarVarDict={}
        self.clinVarGeneDict={}
        self.clinVarVarFound='-'
        self.clinVarGeneFound='-'
        self.clinVarList=[]
        self.clinvarTotalNumVars='-'
        self.clinvarNumP='-'#number of clinvar variants pathogenic
        self.clinvarNumLP='-'#number of clinvar variants likely pathogenic
        self.clinvarNumLB='-'
        self.clinvarNumB='-'
        self.clinvarTitle='-'#title from the flat file
        self.clinvarSignDesc='-'#significance description from the flat file
        self.clinvarCondition='-'
        #dgv
        self.DGVDict={}
        self.dgvDictList=[]
        self.dgvTypeList=[]
        self.dgvSubtypeList=[]
        self.dgvVarFound='-'
        #Decipher
        self.DecipherDict={}
        self.decipherDictList=[]
        self.decipherDeletionObsList=[]
        self.decipherStudyList=[]
        self.decipherVarFound='-'
        #the module scores. The scores will be broken down so they can be added as features
        #curation score
        self.curationScoreTotal='-'
        self.curationScoreHGMD='-'
        self.curationScoreOMIM='-'
        self.curationScoreClinVar='-'
        #conservation scores
        self.conservationScoreTotal='-'
        self.conservationScorePhylop='-'
        self.conservationScoreLRT='-'
        self.conservationScoreGerpPP='-'
        self.conservationScoreCNV='-'
        self.conservationScoreDGV='-'
        self.conservationScoreDecipher='-'
        self.conservationScoreGnomad='-'
        self.conservationScoreGeneConstZ='-'
        self.conservationScoreGeneConstPLi='-'
        self.conservationScoreDomino='-'
        self.conservationScoreOELof='-'
        #module 3 scores
        self.effectOnGeneScore='-'
        self.geneDiseaseAssocScore='-'
        self.modelOrganismScore='-'


#function to get var objects from dataframe read from a VEP annoataed file
def getVarObjFromVEPAnnotDf(df, genomeRef):
    #create the objects
    varObjList = [ Variant() for i in range(len(df.index)) ]
    #create a list of var objects
    i=0
    for index, rows in df.iterrows():
        print('rows:', rows.Uploaded_variation, rows.Location, rows.Gene, rows.Feature_type)
        s=rows.Uploaded_variation.split('_')
        transcriptId=rows.Feature
        chrom=s[0]
        pos=s[1]
        #if it is hg38 get its hg19 coordinates
        if genomeRef=='hg38':
            pass
            varObjList[i].hg38Chrom=chrom
            varObjList[i].hg38pos=pos
            retList=gethg19LocFromHg38(chr, pos)
            retList=[newChrom, newPos]
            varObjList[i].hg19Chrom=retList[0]
            varObjList[i].hg19pos=retList[1]
            chrom=retList[0]
            pos=retList[1]
        else:
            varObjList[i].hg19Chrom=chrom
            varObjList[i].hg19pos=pos

        tmp=s[2]
        s=tmp.split('/')
        refAllele=s[0]
        altAllele=s[1]
        print(chr,pos, tmp, refAllele, altAllele)

        varObjList[i].chrom=chrom
        varObjList[i].pos=pos
        varObjList[i].ref=refAllele
        varObjList[i].alt=altAllele
        #cerate a unique ID since we have transcripts
        varId='_'.join([chrom,pos,refAllele,altAllele,transcriptId])
        varObjList[i].varId=varId
        varObjList[i].varId_dash='-'.join([chrom,pos,refAllele,altAllele])
        varObjList[i].zyg=rows.ZYG
        varObjList[i].geneSymbol=rows.SYMBOL
        varObjList[i].geneEnsId=rows.Gene
        varObjList[i].rsId=rows.Existing_variation
        varObjList[i].GERPpp_RS=rows.GERPpp_RS
        varObjList[i].featureType=rows.Feature_type
        varObjList[i].gnomadAF=rows.gnomAD_AF
        varObjList[i].CLIN_SIG=rows.CLIN_SIG
        varObjList[i].HGVSc=rows.HGVSc
        varObjList[i].HGVSp=rows.HGVSp

        #increment
        i=i+1

    #print('in def len of objects:', len(varObjList))
    return varObjList

#function to get OMIM info using MArrvel API
def getOMIMUsingMarrvelAPI(varObj):
    print('getOMIM:',varObj.varId_dash, varObj.geneSymbol, varObj.rsId)
    url='http://marrvel.org/data/omim/gene/symbol/{}/variant/{}'.format(varObj.geneSymbol,varObj.varId_dash)
    #Example: 6:99365567 T>C / FBXL4
    #url='http://marrvel.org/data/omim/gene/symbol/{}/variant/{}'.format('FBXL4','6-99365567-T-C')
    #print('url:', url)

    req=requests.get(url)
    #print('reqText:', req)
    #print('resqStatusCode:', req.status_code)
    varFound=0
    geneFound=0
    omimDict={}
    if req.status_code!=200:
        #print('OMIM query not found')
        pass
    else:
        #the entry for the gene exists but need to check if the variants exists also using the rs ID
        geneFound=1
        #print('reqText:', req.text)
        dat=json.loads(req.text)
        #search var in the list of dict
        retVal=next((item for item in dat['allelicVariants'] if item["dbSnps"] == varObj.rsId), None)
        if retVal is not None:
            varFound=1
        omimDict=data
        #print('varFoundFlag:',varFoundFlag)


    retList=[varFound, geneFound, omimDict]
    return retList

def getClinVarUsingMarrvelAPI(varObj):
    print('getClinvar:',varObj.varId_dash, varObj.geneSymbol, varObj.rsId)
    #the APIs
    #/clinvar/gene/symbol/:geneSymbol
    #/clinvar/gene/variant/:variant

    #check if varaint annotated
    url='http://marrvel.org/data/clinvar/variant/{}'.format(varObj.varId_dash)
    req=requests.get(url)
    #print('reqText:', req)
    #print('reqText:', req.text)
    varFound=0
    varDict={}
    if req.text:
        varFound=1
        varDict=json.loads(req.text)
    else:
        varFound=0


    #check if gene annotated
    url='http://marrvel.org/data/clinvar/symbol/{}'.format(varObj.geneSymbol)
    req=requests.get(url)
    #print('reqText:', req)
    #print('reqText:', req.text)
    geneFound=0
    geneDict={}
    if req.text:
        geneFound=1
        geneDict=json.loads(req.text)
    else:
        geneFound=0

    #return
    retList=[varFound,varDict, geneFound, geneDict]
    return retList


def getCurationScore(varObj):
    #print('in curation score')
    #get the curation score per each DB
    #omim score
    omimScore=''
    if varObj.omimVarFound==1:
        if varObj.omimSymMatchFlag==1:
            omimScore='High'
        else:
            omimScore='Medium'
    elif varObj.omimGeneFound==1:
        if varObj.omimSymMatchFlag==1:
            omimScore='Medium'
        else:
            omimScore='Low'
    else:
        omimScore='Low'

    #hgmd score
    hgmdScore=''
    if varObj.hgmdVarFound==1:
        if varObj.hgmdSymMatchFlag==1:
            hgmdScore='High'
        else:
            hgmdScore='Medium'
    elif varObj.hgmdGeneFound==1:
        if varObj.hgmdSymMatchFlag==1:
            hgmdScore='Medium'
        else:
            hgmdScore='Low'
    else:
        hgmdScore='Low'


    #clinvar score
    #if anything is not in the path/likely path list then it should have a low score
    clinVarScore=''
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
                                             
    clinVarScore=''
    if varObj.clinVarVarFound==1:
        if varObj.clinvarSignDesc in pathList:
            if varObj.clinVarSymMatchFlag==1:
                clinVarScore='High'
            else:
                clinVarScore='Medium'
        elif varObj.clinvarSignDesc in benignList:
            clinVarScore='Low'
        else:
            if varObj.clinVarSymMatchFlag==1:
                clinVarScore='Medium'
            else:
                clinVarScore='Low'
    elif varObj.clinVarGeneFound==1:
        if varObj.clinVarSymMatchFlag==1:
            clinVarScore='Medium'
        else:
            clinVarScore='Low'
    else:
       clinVarScore='Low' 
    
    
    
    #OLD
    #pathList=['Pathogenic','Likely pathogenic']
    #if varObj.clinvarSignDesc not in pathList:
    #    clinVarScore='Low'
    #elif varObj.clinVarVarFound==1:
    #    if varObj.clinVarSymMatchFlag==1:
    #        clinVarScore='High'
    #    else:
    #        clinVarScore='Medium'
    #elif varObj.clinVarGeneFound==1:
    #    if varObj.clinVarSymMatchFlag==1:
    #        clinVarScore='Medium'
    #    else:
    #        clinVarScore='Low'
    #else:
    #    clinVarScore='Low'

    #print('\tomimScore:', omimScore,'hgmdScore:',hgmdScore,'clinVarScore:',clinVarScore)

    #total score
    totalScore='-'
    #need more work on this total score
    if 0:
        #if omim, clinvar, HGMD avaiable then use omim > hgmd > clinvar for the final score
        varFound=0
        geneFound=0
        if varObj.omimVarFound==1:
            varFound=1
            geneFound=varObj.omimGeneFound
        elif varObj.hgmdVarFound==1:
            varFound=1
            geneFound=varObj.hgmdGeneFound
        elif varObj.clinVarVarFound==1:
            varFound=1
            geneFound=varObj.clinVarGeneFound
        else:
            totalScore='Low'
            retList=[omimScore, hgmdScore, clinVarScore, totalScore]
            return retList

        #print('\tgetCuration func varFound:', varFound,'gene found:', geneFound)
        #continue total score and follow the flowchart
        if varFound==1:
            if symMatchFlag==1:
                totalScore='High'
            elif symMatchFlag==0:
                totalScore='Low'
        elif geneFound==1:
            if symMatchFlag==1:
                totalScore='High'
            elif symMatchFlag==0:
                totalScore='Low'
        else:
            totalScore='Low'

    #print('\t curation total score:', totalScore)
    retList=[omimScore, hgmdScore, clinVarScore, totalScore]
    return retList

def gethg19LocFromHg38(chrom, pos):
    converter = get_lifter('hg38', 'hg19')
    pos = int(pos)
    converter[chrom][pos]

    retVal=converter.query(chrom, pos)
    #print('retVal:', retVal)
    newChrom=retVal[0][0]
    newChrom=newChrom.replace('chr', '')
    newPos=retVal[0][1]
    #print('newChr:',newChrom,'newPos:', newPos)
    retList=[newChrom, newPos]
    return retList


def getDGVUsingMarrvelAPI(varObj):
    #varObj=varObjList[0]
    print('v:',varObj.geneSymbol)
    #DGV returns a list not a dict
    url='http://marrvel.org/data/dgv/variant/{}'.format(varObj.varId_dash)
    req=requests.get(url)
    #print('req:', req.text)
    dgvDictList=[]
    #geneList=[]
    typeList=[]
    subtypeList=[]
    if req.status_code!=200:
        #print('OMIM query not found')
        pass
    else:
        dat=json.loads(req.text)
        dgvDictList=dat
        for row in dat:
            dgvType=row['type']
            typeList.append(dgvType)
            dgvSubtype=row['subType']
            subtypeList.append(dgvSubtype)
            #dgvGeneList=[d['symbol'] for d in row['genes']]
            #geneList.append(dgvGeneList)

    retList=[dgvDictList,typeList,subtypeList]
    return retList


def getDecipherUsingMarrvelAPI(varObj):
    print('v:',varObj.geneSymbol)
    url='http://marrvel.org/data/decipher/variant/{}'.format(varObj.varId_dash)#varId_dash
    #decipher/variant/:variant
    #decipher/genomloc/:hg19Chr/:hg19Start/:hg19Stop
    req=requests.get(url)
    #print('req:', req.text)
    decipherDictList=[]
    decipherDeletionObsList=[]
    decipherStudyList=[]
    if req.status_code!=200:
        #print('Decipher query not found')
        pass
    else:
        dat=json.loads(req.text)
        decipherDictList=dat
        for row in dat:
            #print('\n row:', row)
            deleteObs=row['deletion']['obs']
            decipherDeletionObsList.append(deleteObs)
            deleteStudy=row['study']
            decipherStudyList.append(deleteObs)
            #print('deleteObs:', deleteObs, 'study:',deleteStudy)

    retList=[decipherDictList,decipherDeletionObsList,decipherStudyList]
    return retList

def checkVarGeneFile(inFile, fileType):
    #check if input files are in right format
    pass

def checkPatientFile(inFile, fileType):
    #check if input files are in right format
    pass


def changeVariantFormat(inVar):
    #change the varaint format from 6:99365567 T>C to 6-99365567-T-C
    #replace : > and space with -
    inVar=inVar.replace(' ', '-')
    inVar=inVar.replace(':', '-')
    inVar=inVar.replace('>', '-')
    return inVar

def readVarGeneFile(inFile, fileType):
    #read the var/gene file and return a data frame; if it is varList only then need to call
    if fileType=='varList':
        df=pd.read_csv(inFile, sep='\t')
        print('df head:', df.head(), '\n')
        #change variant format
        df=df.apply(np.vectorize(changeVariantFormat))
        print('df head:', df.head(), '\n')
        return df

def checkVarGnomad(varId):
    url='http://marrvel.org/data/gnomad/variant/{}'.format(varId)
    req=requests.get(url)
    pass


def getVarGroups(varList, groupSize):
    #get a list of variant groups for parallel processing
    for i in range(0, len(varList), groupSize):
        yield varList[i:i + groupSize]
