## UPDATE
#     Changes made by Chaozhong noted as #CL
# from liftover import get_lifter
import pandas as pd
import numpy as np
import pprint
import json
import requests


# class for variant
class Variant:
    """Create a class to create an object for each variant
    Parameters:
        A class name for initiation call
    Output:
        An object of class Variant
    """

    def __init__(self):
        self.chrom = "-"
        self.pos = "-"
        self.start = "-"
        self.stop = "-"
        self.ref = "-"
        self.alt = "-"
        self.hg19Chrom = "-"
        self.hg19Pos = "-"
        self.hg38Chrom = "-"
        self.hg38Pos = "-"
        self.varId = "-"
        self.varId_dash = "-"  # varinat ID separated by dash like '6-99365567-T-C'

        self.patientID = "-"
        self.zyg = "-"
        self.geneSymbol = "-"
        self.geneEnsId = "-"
        self.rsId = "-"
        self.featureType = "-"
        self.gnomadAF = "-"
        self.gnomadAFg = "-"  # genome gnomad score
        # dbnsfp attributes
        self.CLIN_SIG = "-"
        self.CADD_phred = "-"  # phred score from dbnsfp
        self.CADD_PHRED = "-"
        self.GERPpp_RS = "-"  # GERP plus plus
        self.GERPpp_NR = "-"
        self.DANN_score = "-"
        self.FATHMM_pred = "-"
        self.FATHMM_score = "-"
        self.GTEx_V8_gene = "-"
        self.GTEx_V8_tissue = "-"
        self.Polyphen2_HDIV_score = "-"
        self.Polyphen2_HVAR_score = "-"
        self.REVEL_score = "-"
        self.SIFT_score = "-"

        self.clinvar_AlleleID = "-"  # CL added
        self.clinvar_clnsig = "-"
        self.clinvar_CLNREVSTAT = "-"  # CL added
        self.clinvar_CLNSIGCONF = "-"  # CL added

        self.fathmm_MKL_coding_score = "-"
        self.HGVSc = "-"
        self.HGVSp = "-"
        # LRT
        self.LRT_score = "-"
        self.LRT_Omega = "-"
        # Phylo
        self.phyloP100way_vertebrate = "-"
        self.M_CAP_score = "-"
        self.MutationAssessor_score = "-"
        self.MutationTaster_score = "-"
        self.ESP6500_AA_AC = "-"
        self.ESP6500_AA_AF = "-"
        self.ESP6500_EA_AC = "-"
        self.ESP6500_EA_AF = "-"

        # symtom info
        self.SymptomMatched = {
            "omim": 0,
            "clinvar": 0,
            "hgmd": 0,
        }  # list containing omimMatched, clinvarMatched, hgmdMatched
        self.symptomScore = {}  # for OMIM
        self.symptomName = {}  # for OMIM
        self.omimSymptomSimScore = "-"
        self.omimSymMatchFlag = "-"
        self.hgmdSymptomScore = "-"
        self.hgmdSymptomSimScore = "-"
        self.hgmdSymMatchFlag = "-"
        self.clinVarSymMatchFlag = "-"

        # gnomad gene metrics from flat file
        self.gnomadGeneZscore = "-"
        self.gnomadGenePLI = "-"
        self.gnomadGeneOELof = "-"  # O/E lof
        self.gnomadGeneOELofUpper = "-"  # O/E lof upper
        # conseqeunce and impact
        self.IMPACT = "-"
        self.Consequence = "-"
        # omim
        self.omimDict = {}
        self.omimGeneFound = "-"
        self.omimVarFound = "-"
        self.omimList = []
        self.omimGeneDict = {}
        self.omimAlleleDict = {}
        self.phenoList = []
        self.phenoInhList = []
        self.phenoMimList = []
        # HGMD
        self.HGMDDict = {}
        self.hgmdGeneFound = "-"
        self.hgmdVarFound = "-"
        self.hgmdVarPhenIdList = []
        self.hgmdVarHPOIdList = []
        self.hgmdVarHPOStrList = []
        # clinvar
        self.clinVarVarDict = {}
        self.clinVarGeneDict = {}
        self.clinVarVarFound = "-"
        self.clinVarGeneFound = "-"
        self.clinVarList = []
        self.clinvarTotalNumVars = "-"
        self.clinvarNumP = "-"  # number of clinvar variants pathogenic
        self.clinvarNumLP = "-"  # number of clinvar variants likely pathogenic
        self.clinvarNumLB = "-"
        self.clinvarNumB = "-"
        self.clinvarTitle = "-"  # title from the flat file
        self.clinvarSignDesc = "-"  # significance description from the flat file
        self.clinvarCondition = "-"
        # dgv
        self.DGVDict = {}
        self.dgvDictList = []
        self.dgvTypeList = []
        self.dgvSubtypeList = []
        self.dgvVarFound = "-"
        # Decipher
        self.DecipherDict = {}
        self.decipherDictList = []
        self.decipherDeletionObsList = []
        self.decipherStudyList = []
        self.decipherVarFound = "-"
        # the module scores. The scores will be broken down so they can be added as features
        # curation score
        self.curationScoreTotal = "-"
        self.curationScoreHGMD = "-"
        self.curationScoreOMIM = "-"
        self.curationScoreClinVar = "-"
        # conservation scores
        self.conservationScoreTotal = "-"
        self.conservationScorePhylop = "-"
        self.conservationScoreLRT = "-"
        self.conservationScoreGerpPP = "-"
        self.conservationScoreCNV = "-"
        self.conservationScoreDGV = "-"
        self.conservationScoreDecipher = "-"
        self.conservationScoreGnomad = "-"
        self.conservationScoreGeneConstZ = "-"
        self.conservationScoreGeneConstPLi = "-"
        self.conservationScoreDomino = "-"
        self.conservationScoreOELof = "-"
        # module 3 scores
        self.effectOnGeneScore = "-"
        self.geneDiseaseAssocScore = "-"
        self.modelOrganismScore = "-"


# function to get var objects from dataframe read from a VEP annoataed file


# function to get OMIM info using MArrvel API




def getCurationScore(varObj):
    # print('in curation score')
    # get the curation score per each DB
    # omim score
    omimScore = ""
    if varObj.omimVarFound == 1:
        if varObj.omimSymMatchFlag == 1:
            omimScore = "High"
        else:
            omimScore = "Medium"
    elif varObj.omimGeneFound == 1:
        if varObj.omimSymMatchFlag == 1:
            omimScore = "Medium"
        else:
            omimScore = "Low"
    else:
        omimScore = "Low"

    # hgmd score
    hgmdScore = ""
    if varObj.hgmdVarFound == 1:
        if varObj.hgmdSymMatchFlag == 1:
            hgmdScore = "High"
        else:
            hgmdScore = "Medium"
    elif varObj.hgmdGeneFound == 1:
        if varObj.hgmdSymMatchFlag == 1:
            hgmdScore = "Medium"
        else:
            hgmdScore = "Low"
    else:
        hgmdScore = "Low"

    # clinvar score
    # if anything is not in the path/likely path list then it should have a low score
    clinVarScore = ""
    pathList = [
        "Pathogenic",
        "Likely pathogenic",
        "Pathogenic, Affects",
        "Pathogenic/Likely pathogenic, other",
        "Pathogenic/Likely pathogenic",
        "Pathogenic/Likely pathogenic, drug response",
        "Pathogenic/Likely pathogenic, risk factor",
        "Likely pathogenic, drug response",
        "Likely pathogenic, risk factor",
        "Likely pathogenic, association",
        "Likely pathogenic, other",
        "Pathogenic, association, protective",
        "Pathogenic, Affects",
        "Pathogenic, association",
        "Pathogenic, other",
        "Pathogenic, protective",
        "Pathogenic, protective, risk factor",
        "Pathogenic, risk factor",
        "Pathogenic/Likely pathogenic, other",
        "Pathogenic/Likely pathogenic, risk factor",
    ]

    # pathList=['Pathogenic','Likely pathogenic', 'Pathogenic/Likely pathogenic']
    benignList = [
        "Benign",
        "Likely benign",
        "Benign/Likely benign",
        "Benign, association",
        "Benign, drug response",
        "Benign, other",
        "Benign, protective",
        " Benign/Likely benign, Affects",
        "Benign/Likely benign, association",
        "Benign/Likely benign, drug response",
        "Benign/Likely benign, drug response, risk factor",
        "Benign/Likely benign, other",
        "Benign/Likely benign, protective",
        "Benign/Likely benign, protective, risk factor",
        "Benign/Likely benign, risk factor",
        "Likely benign, drug response, other",
        "Likely benign, other",
        "Likely benign, other, risk factor",
        "Likely benign, risk factor",
    ]

    clinVarScore = ""
    if varObj.clinVarVarFound == 1:
        if varObj.clinvarSignDesc in pathList:
            if varObj.clinVarSymMatchFlag == 1:
                clinVarScore = "High"
            else:
                clinVarScore = "Medium"
        elif varObj.clinvarSignDesc in benignList:
            clinVarScore = "Low"
        else:
            if varObj.clinVarSymMatchFlag == 1:
                clinVarScore = "Medium"
            else:
                clinVarScore = "Low"
    elif varObj.clinVarGeneFound == 1:
        if varObj.clinVarSymMatchFlag == 1:
            clinVarScore = "Medium"
        else:
            clinVarScore = "Low"
    else:
        clinVarScore = "Low"

    # OLD
    # pathList=['Pathogenic','Likely pathogenic']
    # if varObj.clinvarSignDesc not in pathList:
    #    clinVarScore='Low'
    # elif varObj.clinVarVarFound==1:
    #    if varObj.clinVarSymMatchFlag==1:
    #        clinVarScore='High'
    #    else:
    #        clinVarScore='Medium'
    # elif varObj.clinVarGeneFound==1:
    #    if varObj.clinVarSymMatchFlag==1:
    #        clinVarScore='Medium'
    #    else:
    #        clinVarScore='Low'
    # else:
    #    clinVarScore='Low'

    # print('\tomimScore:', omimScore,'hgmdScore:',hgmdScore,'clinVarScore:',clinVarScore)

    # total score
    totalScore = "-"
    # need more work on this total score
    if 0:
        # if omim, clinvar, HGMD avaiable then use omim > hgmd > clinvar for the final score
        varFound = 0
        geneFound = 0
        if varObj.omimVarFound == 1:
            varFound = 1
            geneFound = varObj.omimGeneFound
        elif varObj.hgmdVarFound == 1:
            varFound = 1
            geneFound = varObj.hgmdGeneFound
        elif varObj.clinVarVarFound == 1:
            varFound = 1
            geneFound = varObj.clinVarGeneFound
        else:
            totalScore = "Low"
            retList = [omimScore, hgmdScore, clinVarScore, totalScore]
            return retList

        # print('\tgetCuration func varFound:', varFound,'gene found:', geneFound)
        # continue total score and follow the flowchart
        if varFound == 1:
            if symMatchFlag == 1:
                totalScore = "High"
            elif symMatchFlag == 0:
                totalScore = "Low"
        elif geneFound == 1:
            if symMatchFlag == 1:
                totalScore = "High"
            elif symMatchFlag == 0:
                totalScore = "Low"
        else:
            totalScore = "Low"

    # print('\t curation total score:', totalScore)
    retList = [omimScore, hgmdScore, clinVarScore, totalScore]
    return retList


















