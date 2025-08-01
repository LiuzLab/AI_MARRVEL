import re

from .utils_1 import Variant
from .utils_for_marrvel_flatfile import (
    getClinVarUsingMarrvelFlatFile,
    getHGMDUsingFlatFile,
    getAnnotateInfoRow_2,
)


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
    # NOTE(JL): It is old implementation and not used.
    # But left to for tracing purpose. Feel free to remove
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
    varObj = Variant()
    transcriptId = row.Feature
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
    # NOTE(JL): varId_dash should have been renamed. the start value is needed for simple_repeat
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


def getAnnotateInfoRow_3_2(
        varObj,
        decipherSortedDf,
):
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
    return {
        "decipherDictList": retList[0],
        "decipherDeletionObsList": retList[1],
        "decipherStudyList": retList[2],
        "decipherVarFound": retList[3],
    }


def getAnnotateInfoRow_3_3(
        varObj,
        gnomadMetricsGeneSortedDf,
):
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

    return {
        "gnomadGeneZscore": retList[0],
        "gnomadGenePLI": retList[1],
        "gnomadGeneOELof": retList[2],  # O/E lof
        "gnomadGeneOELofUpper": retList[3],  # O/E lof upper
    }


def getAnnotateInfoRow_3_4(
        varObj,
        omimGeneSortedDf,
):
    # get OMIM: 2s
    inputSnpList = []
    if "," in varObj.rsId:
        inputSnpList = varObj.rsId.split(",")
    else:
        inputSnpList = varObj.rsId
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
            if "dbSnps" in a:
                snpList.append(a["dbSnps"])
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

    return {
        "omimVarFound": omimRet[0],
        "omimGeneFound": omimRet[1],
        "omimDict": omimRet[2],
        "omimGeneDict": omimRet[3],
        "omimAlleleDict": omimRet[4],
        "phenoList": omimRet[5],
        "phenoInhList": omimRet[6],
        "phenoMimList": omimRet[7],
    }


def getAnnotateInfoRow_3_5(
        varObj,
        clinvarGeneDf,
        clinvarAlleleDf,
):
    clinVarRet = getClinVarUsingMarrvelFlatFile(
        varObj, clinvarAlleleDf, clinvarGeneDf
    )

    clinVarRet[10] = varObj.clinvar_clnsig  # clinVarRet[10] #CL: changed to clinvar.vcf.gz annotation

    return {
        "clinVarVarFound": clinVarRet[0],
        "clinVarVarDict": clinVarRet[1],
        "clinVarGeneFound": clinVarRet[2],
        "clinVarGeneDict": clinVarRet[3],
        "clinvarTotalNumVars": clinVarRet[4],
        "clinvarNumP": clinVarRet[5],
        "clinvarNumLP": clinVarRet[6],
        "clinvarNumLB": clinVarRet[7],
        "clinvarNumB": clinVarRet[8],
        "clinvarTitle": clinVarRet[9],
        "clinvarSignDesc": clinVarRet[10],
        "clinvarCondition": clinVarRet[11],
    }


def getAnnotateInfoRow_3_6(
        varObj,
        hgmdHPOScoreGeneSortedDf,
):
    hgmdRet = getHGMDUsingFlatFile(varObj, hgmdHPOScoreGeneSortedDf)

    return {
        "hgmdVarFound": hgmdRet[0],
        "hgmdGeneFound": hgmdRet[1],
        "hgmdVarPhenIdList": hgmdRet[2],
        "hgmdVarHPOIdList": hgmdRet[3],
        "hgmdVarHPOStrList": hgmdRet[4],
    }


def getAnnotateInfoRows_3(
        vepDf,
        genomeRef,
        clinvarGeneDf,
        clinvarAlleleDf,
        omimGeneSortedDf,
        omimAlleleList,
        hgmdHPOScoreGeneSortedDf,
        moduleList,
        decipherSortedDf,
        gnomadMetricsGeneSortedDf,
):
    def f1(row):
        return getAnnotateInfoRow_3_1(row, genomeRef)

    def f2(row):
        if "curate" not in moduleList:
            return row
        return getAnnotateInfoRow_3_2(row, decipherSortedDf)

    def f3(row):
        if "conserve" not in moduleList:
            return row
        return getAnnotateInfoRow_3_3(row, gnomadMetricsGeneSortedDf)

    def f4(row):
        if "curate" not in moduleList:
            return row
        return getAnnotateInfoRow_3_4(row, omimGeneSortedDf)

    def f5(row):
        if "curate" not in moduleList:
            return row
        return getAnnotateInfoRow_3_5(row, clinvarGeneDf, clinvarAlleleDf)

    def f6(row):
        if "curate" not in moduleList:
            return row
        return getAnnotateInfoRow_3_6(
            row, hgmdHPOScoreGeneSortedDf
        )

    annotateInfoDf = vepDf.apply(f1, axis=1, result_type='expand')
    df = annotateInfoDf.apply(f2, axis=1, result_type='expand')
    annotateInfoDf[df.columns] = df
    df = annotateInfoDf.apply(f3, axis=1, result_type='expand')
    annotateInfoDf[df.columns] = df
    df = annotateInfoDf.apply(f4, axis=1, result_type='expand')
    annotateInfoDf[df.columns] = df
    df = annotateInfoDf.apply(f5, axis=1, result_type='expand')
    annotateInfoDf[df.columns] = df
    df = annotateInfoDf.apply(f6, axis=1, result_type='expand')
    annotateInfoDf[df.columns] = df
    return annotateInfoDf
