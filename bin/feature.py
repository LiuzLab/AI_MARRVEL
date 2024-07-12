#!/usr/bin/env python3.8
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

# local methods
from annotation.utils_1 import *

# from utils_for_marrvel_annotate_script import *
from annotation.utils_for_symMatch import *
from annotation.utils_for_marrvel_flatfile import *

# from utils_for_marrvel_flatfile_module_3 import *
from functools import partial
from annotation.marrvel_score_recalc import *


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
    # check the user arguments and if anything necessary is missing
    pass


def main():
    # example of running the script:
    # python3 main.py \
    # -patientID UDN630665 -patientHPOsimiOMIM  /houston_10t/dongxue/MARRVEL_AI/BG/HPO2Gene/out_May12021/UDN/HPOsimi_UDN630665_simi_0.tsv \
    # -patientHPOsimiHGMD  /six_tera/chaozhong/module_1/out/UDN/HPOsimi_UDNUDN630665_simi_0.tsv \
    # -varFile /houston_20t/rami/vep_PassVarQual_filt_UDN_Apr2521/540200-UDN630665-P_HGYTFBCXY-2-ID01-vep_gnomADg_gnomADe_splice_clin_hgmd_filt.txt \
    # -inFileType vepAnnotTab -patientFile testPatientFile.txt -patientFileType one -genomeRef hg19 -diseaseInh AD -modules curate,conserve

    parser = argparse.ArgumentParser()
    parser.add_argument("-varFile", "--varFile", help="Proivde input variant file")
    parser.add_argument(
        "-inFileType",
        "--inFileType",
        help="Proivde type of input file:vcf, vepAnnotTab",
    )
    parser.add_argument(
        "-patientFile", "--patientFile", help="Proivde HPO IDs"
    )  # currently not used
    parser.add_argument(
        "-patientFileType", "--patientFileType", help="Proivde type of file:one, two"
    )  # currently not used
    parser.add_argument(
        "-patientHPOsimiOMIM",
        "--patientHPOsimiOMIM",
        help="Proivde patient HPO similarity file-OMIM",
    )
    parser.add_argument(
        "-patientHPOsimiHGMD",
        "--patientHPOsimiHGMD",
        help="Proivde patient HPO similarity file-HGMD",
    )
    parser.add_argument(
        "-diseaseInh", "--diseaseInh", help="Proivde disease Inheritance:AD, AR, XD, XR"
    )
    parser.add_argument(
        "-modules",
        "--modules",
        help="Select modules to run:all,curate,conserve,effectOnGene",
    )
    parser.add_argument(
        "-genomeRef", "--genomeRef", help="Proivde genome ref: hg19, hg38"
    )
    args = parser.parse_args()
    # check the user args
    checkUserArgs(args)
    print("input file:", args.varFile)
    print("type of input file:", args.inFileType)
    print("modules:", args.modules)
    moduleList = args.modules.split(",")
    print("modules list:", moduleList)

    # start time
    start_time = time.time()

    # need to write user options to a log file

    # read the gnomad gene metrics file
    fileName = "annotate/anno_hg19/gnomad.v2.1.1.lof_metrics.by_gene.txt"
    # fileName='/database_test/gnomad.v2.1.1.lof_metrics.by_gene.txt'
    gnomadMetricsGeneDf = pd.read_csv(fileName, sep="\t")
    # print('gnomadMetricsGeneDf shape:', gnomadMetricsGeneDf.shape)

    # set up flags
    # debug flag; if 1 then can check/test the code
    debugFlag = 0
    # print('debugFlag:', debugFlag)
    # use flat files to read the marrvel databases instead of using API
    flatFileFlag = 1
    # use marrvel annotation script and API to get marrvel info. this is slower
    annotateFlag = 0

    # the directory for marrvel database flat file
    # flatFilesDir='/home/dongxue/hdd_10t/MARRVEL_AI/resource/MARRVEL/'

    # initialization
    dgvDf = []
    decipherDf = []
    hgmdDf = []
    omimHPOScoreDf = []
    hgmdHPOScoreDf = []
    clinvarGeneDf = []
    clinvarAlleleDf = []
    omimAlleleList = []

    # read HPO files
    if "curate" in moduleList:
        fileName = args.patientHPOsimiOMIM
        omimHPOScoreDf = pd.read_csv(fileName, sep="\t")
        print("patientHPOsimi-OMIM dimension:", omimHPOScoreDf.shape)

        fileName = args.patientHPOsimiHGMD
        hgmdHPOScoreDf = pd.read_csv(fileName, sep="\t")
        print("patientHPOsimi-HGMD dimension:", hgmdHPOScoreDf.shape)

    # read the clinvar genes file
    if "curate" in moduleList:

        # Q: Why we use the same for hg19 and hg38?
        # A:
        if args.genomeRef == "hg38":
            fileName = "annotate/anno_hg19/gene_clinvar.csv"
        else:
           fileName = "annotate/anno_hg19/gene_clinvar.csv"

        clinvarGeneDf = pd.read_csv(fileName, sep=",")
        # sort by gene name
        clinvarGeneDf.sort_values("symbol", inplace=True)
        clinvarGeneDf.set_index(["symbol"], inplace=True, drop=False)

    # read OMIM
    if "curate" in moduleList:
        # read the OMIM gene file
        fileName = "annotate/anno_hg19/gene_omim.json"

        with open(fileName) as f:
            omimGeneList = json.load(f)

        if debugFlag == 1:
            for omimGeneDict in omimGeneList:
                print("type of omimGeneDict:", type(omimGeneDict))
                print("keys:", omimGeneDict.keys())
                for keyVal in omimGeneDict.keys():
                    print("keyVal:", keyVal)
                    print("\tsubkeys type:", type(omimGeneDict[keyVal]))
                    if isinstance(omimGeneDict[keyVal], list):
                        print("\n\t\tfound list")
                        print("\t\t type:", type(omimGeneDict[keyVal]))

                break

        # read the OMIM allele file
        fileName = "annotate/anno_hg19/omim_alleric_variants.json"

        with open(fileName) as f:
            omimAlleleList = json.load(f)
        if debugFlag == 1:
            print(
                "type of omim allelic :",
                type(omimAlleleList),
                "len:",
                len(omimAlleleList),
            )
            for omimAlleleDict in omimAlleleList:
                print("type of omimAlleleDict:", type(omimAlleleDict))
                # check the keys
                print("keys:", omimAlleleDict.keys())
                for keyVal in omimAlleleDict.keys():
                    print("keyVal:", keyVal)
                    print("\tsubkeys type:", type(omimAlleleDict[keyVal]))
                break

    # read DGV database
    if "conserve" in moduleList:
        print("reading DGV flat file")
        if args.genomeRef == "hg38":
            fileName = "annotate/anno_hg38/dgv.csv"
        else:
            fileName = "annotate/anno_hg19/dgv.csv"
        dgvDf = pd.read_csv(fileName, sep=",")
        # hg19Chr	hg19Start	hg19Stop
        dgvDf = dgvDf.fillna(0)

        # if args.genomeRef == 'hg19':
        dgvDf.columns = ["Chr", "Start", "Stop"] + dgvDf.columns.tolist()[3:]

        dgvDf["Start"] = dgvDf["Start"].astype(int)
        dgvDf["Stop"] = dgvDf["Stop"].astype(int)
        dgvDf["Chr"] = dgvDf["Chr"].replace("X", 23)
        dgvDf["Chr"] = dgvDf["Chr"].replace("Y", 24)
        dgvDf["Chr"] = dgvDf["Chr"].replace("MT", 25)
        dgvDf["Chr"] = dgvDf["Chr"].replace("GL.*", 26, regex=True)

        # make chr as int
        dgvDf["Chr"] = dgvDf["Chr"].astype(int)

        print("finsihed reading DGV")

    # read DECIPHER database
    if "conserve" in moduleList:
        print("reading Decipher flat file")
        if args.genomeRef == "hg38":
            fileName = "annotate/anno_hg38/decipher.csv"
        else:
            fileName = "annotate/anno_hg19/decipher.csv"
        decipherDf = pd.read_csv(fileName, sep=",")
        # hg19Chr,hg19Start,hg19Stop
        decipherDf = decipherDf.fillna(0)

        decipherDf.columns = ["Chr", "Start", "Stop"] + decipherDf.columns.tolist()[3:]

        decipherDf["Start"] = decipherDf["Start"].astype(int)
        decipherDf["Stop"] = decipherDf["Stop"].astype(int)
        decipherDf["Chr"] = decipherDf["Chr"].replace("X", 23)
        decipherDf["Chr"] = decipherDf["Chr"].replace("Y", 24)
        decipherDf["Chr"] = decipherDf["Chr"].replace("MT", 25)
        decipherDf["Chr"] = decipherDf["Chr"].replace("GL.*", "26", regex=True)
        # make chr as int
        decipherDf["Chr"] = decipherDf["Chr"].astype(int)

        """ CL: change column name to be compatible with hg38
        decipherDf['hg19Start']=decipherDf['hg19Start'].astype(int)
        decipherDf['hg19Stop']=decipherDf['hg19Stop'].astype(int)
        decipherDf['hg19Chr']=decipherDf['hg19Chr'].replace('X',23  )
        decipherDf['hg19Chr']=decipherDf['hg19Chr'].replace('Y',24  )
        decipherDf['hg19Chr']=decipherDf['hg19Chr'].replace('MT',25  )
        decipherDf['hg19Chr']=decipherDf['hg19Chr'].replace('GL.*','26',regex=True )
        #make chr as int
        decipherDf['hg19Chr']=decipherDf['hg19Chr'].astype(int)
        """
        print("finsihed reading DECIPHER")

    # read the input variant file
    if args.inFileType == "vepAnnotTab" and flatFileFlag == 1:
        # numHeaderSkip=354
        # this depends on the type of annotation used
        numHeaderSkip = 0
        with open(args.varFile, "r") as F:
            for line in F:
                if line.startswith("##"):
                    numHeaderSkip += 1
                else:
                    break
        print("input annoatated varFile:", args.varFile)
        t1 = time.time()
        varDf = pd.read_csv(
            args.varFile, sep="\t", skiprows=numHeaderSkip, error_bad_lines=False
        )
        # varDf=varDf[0:10]
        # #do this if need to have a small test
        print("shape:", varDf.shape)
        t2 = time.time()
        inputReadTime = t2 - t1
        inputNumRows = len(varDf.index)

        # TODO: change with `.rename`
        if "GERP++_RS" in varDf.columns:
            print("found GERP++RS")
            varDf["GERPpp_RS"] = varDf["GERP++_RS"]
        if "GERP++_NR" in varDf.columns:
            print("found GERP++NR")
            varDf["GERPpp_NR"] = varDf["GERP++_NR"]
        # update column names which have - in it
        if "fathmm-MKL_coding_score" in varDf.columns:
            varDf["fathmm_MKL_coding_score"] = varDf["fathmm-MKL_coding_score"]
        if "M-CAP_score" in varDf.columns:
            varDf["M_CAP_score"] = varDf["M-CAP_score"]

        if "conserve" in moduleList:
            decipherSortedDf = decipherDf.set_index(['Chr', 'Start', 'Stop']).sort_index()
            gnomadMetricsGeneSortedDf = gnomadMetricsGeneDf.groupby('gene').first().sort_index()

        if "curate" in moduleList:
            omimGeneDf = pd.DataFrame(omimGeneList)
            omimGeneSortedDf = omimGeneDf.set_index('geneSymbol').sort_index()

        def f(row):
            return getAnnotateInfoRow_2(
                row,
                args.genomeRef,
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
        if "conserve" in moduleList:
            dgvDfV = dgvDf[['Chr', 'Start', 'Stop']].values.astype('int32')
            def f(varObj):
                # print('\nGetting DGV')
                dgvDictList = []
                dgvTypeList = []
                dgvSubtypeList = []
                dgvVarFound = 0
                dgvType = "-"
                dgvSubtype = "-"
            
                chromVal = int(varObj.chrom)
                posVal = int(varObj.pos)
                startVal = int(varObj.start)
                stopVal = int(varObj.stop)
            
                # CL 03-14-2023: changed column names to be compatible with hg38
                vals = dgvDf[(dgvDfV[:, 0] == chromVal) & (dgvDfV[:, 1] <= startVal) & (dgvDfV[:, 2] >= stopVal)]
                numRows = len(vals.index)
            
                if numRows > 0:
                    dgvVarFound = 1
                    dgvType = vals.iloc[0]["type"]
                    dgvSubtype = vals.iloc[0]["subType"]
            
                dgvTypeList.append(dgvType)
                dgvSubtypeList.append(dgvSubtype)
                return {
                    "dgvDictList": dgvDictList,
                    "dgvTypeList": dgvTypeList,
                    "dgvSubtypeList": dgvSubtypeList,
                    "dgvVarFound": dgvVarFound,
                }
                
            
            resDf = annotateInfoDf.apply(f, axis=1, result_type='expand')
            annotateInfoDf[resDf.columns] = resDf

    # we now have the variant object list "varObjList" with annotations; loop thru it and get the scores
    if "curate" in moduleList:
        df = []
        for i, varObj in annotateInfoDf.iterrows():
            # the curate score is under the utils_1.py file
            omimSymMatch(varObj, omimHPOScoreDf, args.inFileType)
            hgmdSymMatch(varObj, hgmdHPOScoreDf)
            clinVarSymMatch(varObj, args.inFileType)
            # OMIM and clinvar info
            retList = getCurationScore(
                varObj
            )
            df.append({
                "curationScoreOMIM": retList[0],
                "curationScoreHGMD": retList[1],
                "curationScoreClinVar": retList[2],
                "curationScoreTotal": retList[3],

                # Due to direct update inside the functions, here we additionally copy the columns
                "symptomScore": varObj.symptomScore,
                "symptomName": varObj.symptomName,
                "omimSymMatchFlag": varObj.omimSymMatchFlag,
                "omimSymptomSimScore": varObj.omimSymptomSimScore,
                "hgmdSymptomScore": varObj.hgmdSymptomScore,
                "hgmdSymMatchFlag": varObj.hgmdSymMatchFlag,
                "hgmdSymptomSimScore": varObj.hgmdSymptomSimScore,
            })

        df = pd.DataFrame(df)
        annotateInfoDf[df.columns] = df

    # conserve score is under
    if "conserve" in moduleList:
        df = []
        for i, varObj in annotateInfoDf.iterrows():
            retList = getConservationScore(varObj, args.diseaseInh)

            df.append({
                "conservationScoreGnomad": retList[0],
                "conservationScoreDGV": retList[1],
                "conservationScoreOELof": retList[2],
            })
            
        df = pd.DataFrame(df)
        annotateInfoDf[df.columns] = df

    end_time = time.time()
    process_time = end_time - start_time

    print("pipeline time:", process_time)
    # log file for times
    fileName = "log.txt"
    f = open(fileName, "w")
    f.write(
        "Process time:"
        + str(process_time)
        + " seconds and in mins:"
        + str(process_time / 60)
        + "\n"
    )
    f.close()
    print("log file name:", fileName)
    print("input read time:", inputReadTime)
    print("input num rows:", inputNumRows)
    print("m:", moduleList)

    print("Score re-calculation:")
    score = load_raw_matrix(annotateInfoDf)
    score = omimCurate(score)
    score = hgmdCurate(score)
    score = clinvarCurate(score)
    score = conservationCurate(score)

    score.to_csv("scores.csv", index=False)

    exit()


if __name__ == '__main__':
    main()
