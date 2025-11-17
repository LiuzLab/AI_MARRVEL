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
