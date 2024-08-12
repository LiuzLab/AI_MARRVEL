import re

from .utils_1 import Variant
from .utils_for_marrvel_flatfile import (
    getClinVarUsingMarrvelFlatFile,
    getHGMDUsingFlatFile,
    getAnnotateInfoRow_2,
)

# It is as it is, and supposed to be refactored.
def getAnnotateInfoRows_2(
        varDf,
        genomeRef,
        clinvarGeneDf,
        clinvarAlleleDf,
        omimGeneSortedDf,
        omimAlleleList,
        hgmdHPOScoreDf,
        moduleList,
        decipherSortedDf,
        gnomadMetricsGeneSortedDf,
):
    def f(row):
        return getAnnotateInfoRow_2(
            row,
            genomeRef,
            clinvarGeneDf,
            clinvarAlleleDf,
            omimGeneSortedDf,
            omimAlleleList,
            hgmdHPOScoreDf,
            moduleList,
            decipherSortedDf,
            gnomadMetricsGeneSortedDf,
        )
    annotateInfoDf = varDf.apply(f, axis=1, result_type='expand')
    return annotateInfoDf


def getAnnotateInfoRow_3_1(row, genomeRef):
    # CL 03-14-2023: commented all printing lines
    # print('type of row:', type(row))
    varObj = Variant()
    transcriptId = row.Feature
    # s=row.Uploaded_variation.split('_') '1_10204_-/T' 1_1588250_T_A
    ####row[0]: 21_11039079_C/A
    ####s: ['21', '11039079', 'C/A']
    # print('row[0]:', row[0])

    # two ways of input of first column either 1_1588250_T_A OR 21_11039079_C/A, so use the option flag
    optFlag = 0
    if row[0].find("/") != -1:
        optFlag = 1

    if optFlag == 0:
        s = row[0].split("_")
        # print('s:', s)
        chrom = s[0]
        pos = int(s[1])
        ref = s[2]
        alt = s[3]
    elif optFlag == 1:
        s = row[0].split("_")
        # print('s:', s)
        chrom = s[0]
        pos = int(s[1])
        s = s[2].split("/")
        ref = s[0]
        alt = s[1]

    # get the start and stop from second column like '1:10203-10204'
    if "-" in row[1]:
        s = row[1].split(":")
        tmp = s[1]
        s = tmp.split("-")
        # print('s:',s)
        start = int(s[0])
        stop = int(s[1])
    else:
        # start and stop the same
        s = row[1].split(":")
        start = int(s[1])
        stop = int(s[1])

    # print('chrom:', chrom,'pos:',pos,'ref:',ref,'alt:',alt,'start:',start,'stop:',stop)
    # change chrom X and Y and MT to numbers
    if chrom == "X":
        chrom = 23
    elif chrom == "Y":
        chrom = 24
    elif chrom == "MT":
        chrom = 25
    elif re.search(r"GL", chrom):
        chrom = 26
    chrom = int(chrom)

    # if it is hg38 get its hg19 coordinates
    # CL 03-14-2023: we have separate database for hg19 and hg38,
    #                we don't need to use LiftOver which is inaccurate
    #                related codes commented and modified
    if genomeRef == "hg38":
        varObj.hg38Chrom = chrom
        varObj.hg38Pos = pos
        varObj.chrom = chrom
        varObj.pos = pos
        varObj.start = start
        varObj.stop = stop

        """
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
        """
    else:
        varObj.hg19Chrom = chrom
        varObj.hg19Pos = pos
        varObj.chrom = chrom
        varObj.pos = pos
        varObj.start = start
        varObj.stop = stop

    geneSymbol = row.SYMBOL
    # print('gene:', geneSymbol)
    varObj.geneSymbol = geneSymbol
    varObj.CADD_phred = row.CADD_phred
    varObj.CADD_PHRED = row.CADD_PHRED

    # assign
    varObj.ref = ref
    varObj.alt = alt
    varObj.varId_dash = "-".join([str(chrom), str(start), ref, alt])
    # print('varId dash:', varObj.varId_dash)
    varId = "_".join([str(chrom), str(pos), ref, alt, transcriptId])
    varObj.varId = varId
    if "ZYG" in row:
        varObj.zyg = row.ZYG
    varObj.geneEnsId = row.Gene
    varObj.rsId = row.Existing_variation
    varObj.GERPpp_RS = row.GERPpp_RS
    varObj.featureType = row.Feature_type
    varObj.gnomadAF = row.gnomAD_AF
    varObj.gnomadAFg = row.gnomADg_AF
    varObj.CLIN_SIG = row.CLIN_SIG  # CL: useless but kept for now
    varObj.LRT_Omega = row.LRT_Omega
    varObj.LRT_score = row.LRT_score
    varObj.phyloP100way_vertebrate = row.phyloP100way_vertebrate
    varObj.IMPACT = row.IMPACT
    varObj.Consequence = row.Consequence
    varObj.HGVSc = row.HGVSc
    varObj.HGVSp = row.HGVSp
    # dbnsfp attributes
    varObj.GERPpp_NR = row.GERPpp_NR
    varObj.DANN_score = row.DANN_score
    varObj.FATHMM_pred = row.FATHMM_pred
    varObj.FATHMM_score = row.FATHMM_score
    varObj.GTEx_V8_gene = row.GTEx_V8_gene
    varObj.GTEx_V8_tissue = row.GTEx_V8_tissue
    varObj.Polyphen2_HDIV_score = row.Polyphen2_HDIV_score
    varObj.Polyphen2_HVAR_score = row.Polyphen2_HVAR_score
    varObj.REVEL_score = row.REVEL_score
    varObj.SIFT_score = row.SIFT_score

    varObj.clinvar_AlleleID = row.clinvar  # Clinvar allele ID from clinvar.vcf.gz
    varObj.clinvar_clnsig = (
        row.clinvar_CLNSIG
    )  # CL: Clinvar SIG from clinvar.vcf.gz
    # varObj.clinvar_clnsig = row.clinvar_clnsig #CL: Clinvar SIG from VEP, deleted
    varObj.clinvar_CLNREVSTAT = (
        row.clinvar_CLNREVSTAT
    )  # CL: Clinvar STAT from clinvar.vcf.gz, for interface only
    varObj.clinvar_CLNSIGCONF = (
        row.clinvar_CLNSIGCONF
    )  # CL: Clinvar SIGCONF from clinvar.vcf.gz
    varObj.clin_code = row.clinvar_CLNSIG  # CL: feature name for ai

    varObj.fathmm_MKL_coding_score = row.fathmm_MKL_coding_score
    varObj.LRT_score = row.LRT_score
    varObj.LRT_Omega = row.LRT_Omega
    varObj.phyloP100way_vertebrate = row.phyloP100way_vertebrate
    varObj.M_CAP_score = row.M_CAP_score
    varObj.MutationAssessor_score = row.MutationAssessor_score
    varObj.MutationTaster_score = row.MutationTaster_score
    varObj.ESP6500_AA_AC = row.ESP6500_AA_AC
    varObj.ESP6500_AA_AF = row.ESP6500_AA_AF
    varObj.ESP6500_EA_AC = row.ESP6500_EA_AC
    varObj.ESP6500_EA_AF = row.ESP6500_EA_AF

    varObj.VARIANT_CLASS = row.VARIANT_CLASS
    varObj.Feature = row.Feature
    varObj.hom = row.gnomADg_controls_nhomalt
    varObj.hgmd_id = row.hgmd  # CL added
    varObj.hgmd_symbol = row.hgmd_GENE  # CL added
    varObj.hgmd_rs = row.hgmd_RANKSCORE
    varObj.hgmd_PHEN = row.hgmd_PHEN  # CL added
    varObj.hgmd_CLASS = row.hgmd_CLASS  # CL added

    if row.clinvar_CLNSIGCONF != "-":
        clin_dict = dict()
        for ro in row.clinvar_CLNSIGCONF.split("|_"):
            temp = ro.split("(")
            clin_dict[temp[0]] = int(temp[1][0])
        PLP_sum = clin_dict.get("Pathogenic", 0) + clin_dict.get(
            "Likely_pathogenic", 0
        )
        varObj.clin_dict = clin_dict
        varObj.clin_PLP = PLP_sum
        varObj.clin_PLP_perc = PLP_sum / sum(clin_dict.values())
    else:
        if "benign" in row.clinvar_clnsig.lower():
            varObj.clin_PLP_perc = 0
        elif "pathogenic" in row.clinvar_clnsig.lower():
            varObj.clin_PLP_perc = 1
        else:
            varObj.clin_PLP_perc = "-"
        varObj.clin_PLP = "-"
        varObj.clin_dict = "-"

    if row.SpliceAI_pred != "-":
        varObj.spliceAI = row.SpliceAI_pred
        temp = row.SpliceAI_pred.split("|")
        varObj.spliceAImax = max(
            float(temp[1]), float(temp[2]), float(temp[3]), float(temp[4])
        )
    else:
        varObj.spliceAI = "-"
        varObj.spliceAImax = "-"

    return vars(varObj)


def getAnnotateInfoRow_3_x(
        varObj,
        genomeRef,
        clinvarGeneDf,
        clinvarAlleleDf,
        omimGeneSortedDf,
        omimAlleleList,
        hgmdHPOScoreDf,
        moduleList,
        decipherSortedDf,
        gnomadMetricsGeneSortedDf,
):
    if "conserve" in moduleList:
        # get decipher: 0.6s
        decipherDictList = []
        decipherDeletionObsList = []
        decipherStudyList = []
        decipherVarFound = 0
        deletionObs = "-"
        # get the varaint object info from varObj
        chromVal = int(varObj.chrom)
        posVal = int(varObj.pos)
        startVal = int(varObj.start)
        stopVal = int(varObj.stop)

        # CL 03-14-2023: changed column names to be compatible with hg38
        # vals=decipherDf[ ( decipherDf['hg19Chr'] == chromVal ) & ( decipherDf['hg19Start']==startVal ) & (decipherDf['hg19Stop']==stopVal) ]
        if (chromVal, startVal, stopVal) in decipherSortedDf:
            vals = decipherSortedDf.loc[(chromVal, startVal, stopVal)]

            decipherVarFound = 1
            deletionObs = vals.iloc[0]["deletion.obs"]
            decipherDeletionObsList.append(deletionObs)

        # print('\tchrom:', chromVal,'posVal:', posVal,'start:', startVal,'stopVal:', stopVal)
        # print('\tdecipherVarFound:',decipherVarFound,'decipherDeletionObs:', deletionObs)
        retList = [
            decipherDictList,
            decipherDeletionObsList,
            decipherStudyList,
            decipherVarFound,
        ]

        # [decipherDictList,decipherDeletionObsList,decipherStudyList, decipherVarFound]
        varObj.decipherDictList = retList[0]
        varObj.decipherDeletionObsList = retList[1]
        varObj.decipherStudyList = retList[2]
        varObj.decipherVarFound = retList[3]

        # get gnomad gene metrics from gnomad file: 3.1s
        if varObj.geneSymbol in gnomadMetricsGeneSortedDf.index:  # pLI, oe_lof, oe_lof_upper,mis_z
            val = gnomadMetricsGeneSortedDf.loc[varObj.geneSymbol]
            gnomadGeneZscore = val["mis_z"]
            gnomadGenePLI = val["pLI"]
            gnomadGeneOELof = val["oe_lof"]
            gnomadGeneOELofUpper = val["oe_lof_upper"]
        else:
            # get the values
            gnomadGeneZscore = "-"
            gnomadGenePLI = "-"
            gnomadGeneOELof = "-"
            gnomadGeneOELofUpper = "-"

        retList = [gnomadGeneZscore, gnomadGenePLI, gnomadGeneOELof, gnomadGeneOELofUpper]

        # [decipherDictList,decipherDeletionObsList,decipherStudyList, decipherVarFound]
        varObj.gnomadGeneZscore = retList[0]
        varObj.gnomadGenePLI = retList[1]
        varObj.gnomadGeneOELof = retList[2]  # O/E lof
        varObj.gnomadGeneOELofUpper = retList[3]  # O/E lof upper

    if "curate" in moduleList:
        # get OMIM: 2s
        # print('\nGetting OMIM')
        # varObj.omimList=jsonDict['omim']
        # retList=[varFound, geneFound, omimDict, omimGeneDict, omimAlleleDict]
        inputSnpList = []
        if "," in varObj.rsId:
            inputSnpList = varObj.rsId.split(",")
        else:
            inputSnpList = varObj.rsId
        # print('\tinputSnpList:', inputSnpList)
        varFound = 0
        geneFound = 0
        omimDict = {}
        omimGeneDict = {}
        omimAlleleDict = {}
        phenoList = []
        phenoInhList = []
        phenoMimList = []
        # check gene
        # keys: dict_keys(['phenotypes', 'allelicVariants', 'mimNumber', 'status', 'title', 'description', 'geneEntrezId', 'geneSymbol'])
        if varObj.geneSymbol in omimGeneSortedDf.index:
            # print('\tgene:', varObj.geneSymbol, 'found')
            geneFound = 1
            omimGeneDict = omimGeneSortedDf.loc[varObj.geneSymbol]
            snpList = []
            for a in omimGeneDict["allelicVariants"]:
                # print('a:', a)
                # print('type:', type(a))
                if "dbSnps" in a:
                    snpList.append(a["dbSnps"])
            # print('\tsnpList:', snpList)
            # print('\tlen snpList:', len(snpList))
            # check if input snpID matches the OMIM one
            set1 = set(inputSnpList)
            set2 = set(snpList)
            if set1.intersection(set2):
                varFound = 1
            else:
                varFound = 0

            # get disease info from OMIM
            # print('\tphenotypes:', type(omimGeneDict['phenotypes']), ' len:', len(omimGeneDict['phenotypes']) )
            for a in omimGeneDict["phenotypes"]:
                # print('type:', type(a))
                pheno = a["phenotype"]
                if "phenotypeMimNumber" in a:
                    phenoMim = a["phenotypeMimNumber"]
                else:
                    phenoMim = "-"
                if "phenotypeInheritance" in a:
                    phenoInh = a["phenotypeInheritance"]
                else:
                    phenoInh = "-"
                phenoList.append(pheno)
                phenoInhList.append(phenoInh)
                phenoMimList.append(str(phenoMim))
                # print('phenotype:', pheno,phenoMim,phenoInh)

        # print('\tvarFound:', varFound)
        # print('\tphenoList:', phenoList)
        # print('\tphenoInhList:', phenoInhList)
        # print('\tphenoMimList:', phenoMimList)

        omimRet = [
            varFound,
            geneFound,
            omimDict,
            omimGeneDict,
            omimAlleleDict,
            phenoList,
            phenoInhList,
            phenoMimList,
        ]

        varObj.omimVarFound = omimRet[0]
        varObj.omimGeneFound = omimRet[1]
        varObj.omimDict = omimRet[2]
        varObj.omimGeneDict = omimRet[3]
        varObj.omimAlleleDict = omimRet[4]
        varObj.phenoList = omimRet[5]
        varObj.phenoInhList = omimRet[6]
        varObj.phenoMimList = omimRet[7]
        # print('OMIM res:')
        # print('\tgeneFound:',varObj.omimGeneFound,'varFound:',varObj.omimVarFound )

        # get clinvar: 0.1s
        # print('\nReading clinVar')
        clinVarRet = getClinVarUsingMarrvelFlatFile(
            varObj, clinvarAlleleDf, clinvarGeneDf
        )
        varObj.clinVarVarFound = clinVarRet[0]
        varObj.clinVarVarDict = clinVarRet[1]
        varObj.clinVarGeneFound = clinVarRet[2]
        varObj.clinVarGeneDict = clinVarRet[3]
        varObj.clinvarTotalNumVars = clinVarRet[4]
        varObj.clinvarNumP = clinVarRet[5]
        varObj.clinvarNumLP = clinVarRet[6]
        varObj.clinvarNumLB = clinVarRet[7]
        varObj.clinvarNumB = clinVarRet[8]
        varObj.clinvarTitle = clinVarRet[9]
        varObj.clinvarSignDesc = (
            varObj.clinvar_clnsig
        )  # clinVarRet[10] #CL: changed to clinvar.vcf.gz annotation
        varObj.clinvarCondition = clinVarRet[11]
        # print('clinVar res:')
        """
        if debugFlag==1:
            print('\tgeneFound::',varObj.clinVarGeneFound,'varFound:',varObj.clinVarVarFound)
            print('\tnumVars:',varObj.clinvarTotalNumVars,'numPathologic:',varObj.clinvarNumP,'numBenign:',varObj.clinvarNumB)
            print('\tsignDesc:', varObj.clinvarSignDesc)
        """

        # get HGMD: 0.3s
        if "curate" in moduleList:
            # print('\nReading HGMD')
            hgmdRet = getHGMDUsingFlatFile(varObj, hgmdHPOScoreDf)
            # hgmdVarFound,hgmdGeneFound,hgmdVarPhenIdList,hgmdVarHPOIdList,hgmdVarHPOStrList
            varObj.hgmdVarFound = hgmdRet[0]
            varObj.hgmdGeneFound = hgmdRet[1]
            varObj.hgmdVarPhenIdList = hgmdRet[2]
            varObj.hgmdVarHPOIdList = hgmdRet[3]
            varObj.hgmdVarHPOStrList = hgmdRet[4]
            # print('HGMD results:')
            # print('\thgmdVarFound:',varObj.hgmdVarFound,'hgmdGeneFound:',varObj.hgmdGeneFound,
            #      'hgmdVarPhenIdList:',varObj.hgmdVarPhenIdList,'hgmdVarHPOIdList:',
            #      varObj.hgmdVarHPOIdList,
            #      'hgmdVarHPOStrList:',varObj.hgmdVarHPOStrList)

        return {
            "hg19Chrom": varObj.hg19Chrom,
            "hg19Pos": varObj.hg19Pos,
            "chrom": varObj.chrom,
            "pos": varObj.pos,
            "start": varObj.start,
            "stop": varObj.stop,
            "geneSymbol": varObj.geneSymbol,
            "CADD_phred": varObj.CADD_phred,
            "CADD_PHRED": varObj.CADD_PHRED,
            "ref": varObj.ref,
            "alt": varObj.alt,
            "varId": varObj.varId,
            "ZYG": varObj.zyg,
            "HGVSc": varObj.HGVSc,
            "HGVSp": varObj.HGVSp,
            "Gene": varObj.geneEnsId,
            "Existing_variation": varObj.rsId,
            "GERPpp_RS": varObj.GERPpp_RS,
            "Feature_type": varObj.featureType,
            "gnomadAF": varObj.gnomadAF,
            "gnomadAFg": varObj.gnomadAFg,
            "CLIN_SIG": varObj.CLIN_SIG,
            "LRT_Omega": varObj.LRT_Omega,
            "LRT_score": varObj.LRT_score,
            "phyloP100way_vertebrate": varObj.phyloP100way_vertebrate,
            # dbnsfp attributes
            "GERPpp_NR": varObj.GERPpp_NR,
            "DANN_score": varObj.DANN_score,
            "FATHMM_pred": varObj.FATHMM_pred,
            "FATHMM_score": varObj.FATHMM_score,
            "GTEx_V8_gene": varObj.GTEx_V8_gene,
            "GTEx_V8_tissue": varObj.GTEx_V8_tissue,
            "Polyphen2_HDIV_score": varObj.Polyphen2_HDIV_score,
            "Polyphen2_HVAR_score": varObj.Polyphen2_HVAR_score,
            "REVEL_score": varObj.REVEL_score,
            "SIFT_score": varObj.SIFT_score,
            "clinvar_AlleleID": varObj.clinvar_AlleleID,  # Clinvar allele ID from clinvar.vcf.gz
            "clinvar_clnsig": varObj.clinvar_clnsig,  # CL: Clinvar SIG from clinvar.vcf.gz
            "clinvar_CLNREVSTAT": varObj.clinvar_CLNREVSTAT,  # CL: Clinvar STAT from clinvar.vcf.gz, for interface only
            "clinvar_CLNSIGCONF": varObj.clinvar_CLNSIGCONF,  # CL: Clinvar SIGCONF from clinvar.vcf.gz
            "clin_code": varObj.clin_code,  # CL: feature for ai
            "fathmm_MKL_coding_score": varObj.fathmm_MKL_coding_score,
            "LRT_score": varObj.LRT_score,
            "LRT_Omega": varObj.LRT_Omega,
            "phyloP100way_vertebrate": varObj.phyloP100way_vertebrate,
            "M_CAP_score": varObj.M_CAP_score,
            "MutationAssessor_score": varObj.MutationAssessor_score,
            "MutationTaster_score": varObj.MutationTaster_score,
            "ESP6500_AA_AC": varObj.ESP6500_AA_AC,
            "ESP6500_AA_AF": varObj.ESP6500_AA_AF,
            "ESP6500_EA_AC": varObj.ESP6500_EA_AC,
            "ESP6500_EA_AF": varObj.ESP6500_EA_AF,
            # dbnsfp
            "gnomadGeneZscore": varObj.gnomadGeneZscore,
            "gnomadGenePLI": varObj.gnomadGenePLI,
            "gnomadGeneOELof": varObj.gnomadGeneOELof,  # O/E lof
            "gnomadGeneOELofUpper": varObj.gnomadGeneOELofUpper,  # O/E lof upper,
            "IMPACT": varObj.IMPACT,
            "Consequence": varObj.Consequence,
            "omimVarFound": varObj.omimVarFound,
            "omimGeneFound": varObj.omimGeneFound,
            "omimDict": varObj.omimDict,
            "omimGeneDict": varObj.omimGeneDict,
            "omimAlleleDict": varObj.omimAlleleDict,
            "phenoList": varObj.phenoList,
            "phenoInhList": varObj.phenoInhList,
            "phenoMimList": varObj.phenoMimList,
            "clinVarVarFound": varObj.clinVarVarFound,
            "clinVarVarDict": varObj.clinVarVarDict,
            "clinVarGeneFound": varObj.clinVarGeneFound,
            "clinVarGeneDict": varObj.clinVarGeneDict,
            "clinvarTotalNumVars": varObj.clinvarTotalNumVars,
            "clinvarNumP": varObj.clinvarNumP,
            "clinvarNumLP": varObj.clinvarNumLP,
            "clinvarNumLB": varObj.clinvarNumLB,
            "clinvarNumB": varObj.clinvarNumB,
            "clinvarTitle": varObj.clinvarTitle,
            "clinvarSignDesc": varObj.clinvarSignDesc,
            "clinvarCondition": varObj.clinvarCondition,
            "hgmdVarFound": varObj.hgmdVarFound,
            "hgmdGeneFound": varObj.hgmdGeneFound,
            "hgmdVarPhenIdList": varObj.hgmdVarPhenIdList,
            "hgmdVarHPOIdList": varObj.hgmdVarHPOIdList,
            "hgmdVarHPOStrList": varObj.hgmdVarHPOStrList,
            "varId_dash": varObj.varId_dash,
            "dgvDictList": varObj.dgvDictList,
            "dgvTypeList": varObj.dgvTypeList,
            "dgvSubtypeList": varObj.dgvSubtypeList,
            "dgvVarFound": varObj.dgvVarFound,
            "decipherDictList": varObj.decipherDictList,
            "decipherDeletionObsList": varObj.decipherDeletionObsList,
            "decipherStudyList": varObj.decipherStudyList,
            "decipherVarFound": varObj.decipherVarFound,
            "gnomadGeneZscore": varObj.gnomadGeneZscore,
            "gnomadGenePLI": varObj.gnomadGenePLI,
            "gnomadGeneOELof": varObj.gnomadGeneOELof,
            "gnomadGeneOELofUpper": varObj.gnomadGeneOELofUpper,
            # symptom
            "SymptomMatched": varObj.SymptomMatched,
            "symptomScore": varObj.symptomScore,
            "symptomName": varObj.symptomName,
            "omimSymptomSimScore": varObj.omimSymptomSimScore,
            "omimSymMatchFlag": varObj.omimSymMatchFlag,
            "hgmdSymptomScore": varObj.hgmdSymptomScore,
            "hgmdSymptomSimScore": varObj.hgmdSymptomSimScore,
            "hgmdSymMatchFlag": varObj.hgmdSymMatchFlag,
            "clinVarSymMatchFlag": varObj.clinVarSymMatchFlag,
            "VARIANT_CLASS": varObj.VARIANT_CLASS,
            "Feature": varObj.Feature,
            "hom": varObj.hom,
            "hgmd_rs": varObj.hgmd_rs,
            "hgmd_id": varObj.hgmd_id,  # CL added
            "hgmd_symbol": varObj.hgmd_symbol,  # CL added
            "hgmd_PHEN": varObj.hgmd_PHEN,  # CL added
            "hgmd_CLASS": varObj.hgmd_CLASS,  # CL added
            "clin_dict": varObj.clin_dict,
            "clin_PLP": varObj.clin_PLP,
            "clin_PLP_perc": varObj.clin_PLP_perc,
            "spliceAI": varObj.spliceAI,
            "spliceAImax": varObj.spliceAImax,

            "zyg": varObj.zyg,
            'geneEnsId': varObj.geneEnsId,
            'rsId': varObj.rsId
        }




def getAnnotateInfoRows_3(
        varDf,
        genomeRef,
        clinvarGeneDf,
        clinvarAlleleDf,
        omimGeneSortedDf,
        omimAlleleList,
        hgmdHPOScoreDf,
        moduleList,
        decipherSortedDf,
        gnomadMetricsGeneSortedDf,
):
    def f1(row):
        return getAnnotateInfoRow_3_1(row, genomeRef)

    def fx(row):
        return getAnnotateInfoRow_3_x(
            row,
            genomeRef,
            clinvarGeneDf,
            clinvarAlleleDf,
            omimGeneSortedDf,
            omimAlleleList,
            hgmdHPOScoreDf,
            moduleList,
            decipherSortedDf,
            gnomadMetricsGeneSortedDf,
        )

    df = varDf.apply(f1, axis=1, result_type='expand')
    annotateInfoDf = df.apply(fx, axis=1, result_type='expand')
    return annotateInfoDf
