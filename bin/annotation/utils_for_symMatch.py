## UPDATE
#     Changes made by Chaozhong noted as #CL

from annotation.utils_1 import *


def omimSymMatch(varObj, omimHPOScoreDf, inFileType):
    """
    Find OMIM symptom match score
    Param:
    omimHPOScoreDf: this is read from user input OMIM symptom mathc file. Foe example:/houston_10t/dongxue/MARRVEL_AI/BG/HPO2Gene/out_May12021/UDN/HPOsimi_UDN630665_simi_0.tsv
    inFileType:the type of input file. Will be used later.
    Return:
    Calculate the variant symptom score
    """
    # print('\nin omimSymMatch')
    # print('type varObj:', type(varObj))
    # print('\tvar:', varObj.varId_dash)
    simScore = 0
    if "phenotypes" in varObj.omimGeneDict.keys():
        # print('\t phenotypes key found')
        for pheno in varObj.omimGeneDict["phenotypes"]:
            if "phenotypeMimNumber" in pheno.keys():
                PminNum = pheno["phenotypeMimNumber"]
                # print('\tPminNum:', PminNum)
                # print('\ttmp1:', omimHPOScoreDf[omimHPOScoreDf['Pheno_ID'] == PminNum])
                valDf = omimHPOScoreDf[omimHPOScoreDf["Pheno_ID"] == PminNum]
                numRows = len(valDf.index)
                if numRows != 0:
                    simScore = omimHPOScoreDf[
                        omimHPOScoreDf["Pheno_ID"] == PminNum
                    ].iloc[0]["Similarity_Score"]
                else:
                    simScore = 0
                # print('\tsimScore:', simScore)
            else:
                simScore = 0
            if simScore >= 0.2:  # change the threshold to a user arg later
                varObj.symptomScore[PminNum] = simScore
                symptomList = (
                    omimHPOScoreDf[omimHPOScoreDf["Pheno_ID"] == PminNum]
                    .iloc[0]["Disease_Name"]
                    .split(";")
                )
                symptomList = [
                    i.strip().upper() for i in symptomList
                ]  # remove space and tranfer all words to capital
                varObj.symptomName[PminNum] = symptomList
                varObj.omimSymMatchFlag = 1
            else:
                continue
    # store the sim score we get
    varObj.omimSymptomSimScore = simScore
    # print('omimSymMatch results:')
    # print('\tsymtomSimScore:', varObj.omimSymptomSimScore)
    # print('\tsymptomScore:', varObj.symptomScore)
    # print('\tsymptomName:', varObj.symptomName)
    # print('\tomimSymMatchFlag:', varObj.omimSymMatchFlag)


def hgmdSymMatch(varObj, hgmdHPOScoreAccSortedDf, hgmdHPOScoreGeneSortedDf):
    """
    Find HGMD symptom match score
    Param:
    omimHPOScoreDf: this is read from user input OMIM symptom mathc file. For example:/six_tera/chaozhong/module_1/out/UDN/HPOsimi_UDNUDN630665_simi_0.tsv
    inFileType:the type of input file. Will be used later.
    Return:
    Calculate the variant symptom score
    """

    # print('\nin HGMDSymMatch')
    # print('\tvar:', varObj.varId_dash)
    hgmdSymptomSimScore = "-"
    if varObj.hgmd_id in hgmdHPOScoreAccSortedDf.index:
        varScore = hgmdHPOScoreAccSortedDf.loc[varObj.hgmd_id].Similarity_Score

        varObj.hgmdSymptomScore = varScore
        if varScore >= 0.2:
            varObj.hgmdSymMatchFlag = 1
        hgmdSymptomSimScore = varScore
    elif varObj.hgmdGeneFound:
        geneScore = hgmdHPOScoreGeneSortedDf.loc[varObj.geneSymbol].Similarity_Score

        if geneScore >= 0.2:
            varObj.hgmdSymMatchFlag = 1
        hgmdSymptomSimScore = geneScore

    varObj.hgmdSymptomSimScore = hgmdSymptomSimScore
    # print('hgmdSymMatch results:')
    # print('\thgmdSymMatchFlag:', varObj.hgmdSymMatchFlag)
    # print('\thgmdSymptomSimScore:', varObj.hgmdSymptomSimScore)


def clinVarSymMatch(varObj, inFileType):
    """
    Find clinvar symptom match score
    Param:
    Return:
    Calculate the variant symptom score
    """
    # print('\nin clinvarSymMatch')
    if varObj.clinVarVarFound:
        # pheno = varObj.clinVarVarDict['condition'].strip().upper()
        # print('clinvar condition:', varObj.clinvarCondition)
        if type(varObj.clinvarCondition) is str:
            pheno = varObj.clinvarCondition.strip().upper()
            for PminNum in varObj.symptomName.keys():
                if (
                    pheno in varObj.symptomName[PminNum]
                    or "#%s %s" % (PminNum, pheno) in varObj.symptomName[PminNum]
                ):
                    varObj.clinVarSymMatchFlag = 1
        else:
            varObj.clinVarSymMatchFlag = 0

    # Need more checking on this one per gene condition
    if varObj.clinVarGeneFound:
        for var in varObj.clinVarGeneDict:
            pheno = var["condition"].strip().upper()
            for PminNum in varObj.symptomName.keys():
                if (
                    pheno in varObj.symptomName[PminNum]
                    or "#%s %s" % (PminNum, pheno) in varObj.symptomName[PminNum]
                ):
                    varObj.clinVarGeneSymMatchFlag = 1

    # print('\tvarObj.clinVarVarSymMatchFlag:', varObj.clinVarSymMatchFlag)
