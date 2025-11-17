## UPDATE
#     Changes made by Chaozhong noted as #CL

from annotation.variant_models import *
import numpy as np
import time
import re


def getConservationScore(varObj, diseaseInh):
    """
    Get conservation score
    Params:
    varObj:variant object
    diseaseInh: the disease inheritance
    Return:
    return three conservation scores: conservationScoreGnomad, conservationScoreDGV, conservationScoreOELof
    """
    # print('in conservation score call')
    # create conservation scores
    # set up thresholds
    gnomadGFreqCut = 0.01  # genome g
    gnomadEFreqCut = 0.01  # exome e
    # GERP ranges
    gerpPPCut = 4
    # LRT cut
    # phylo cut

    # gnomad
    conservationScoreGnomad = "-"
    gnomadAFVal = _parse_numeric(str(varObj.gnomadAF), "min")
    gnomadAFgVal = _parse_numeric(str(varObj.gnomadAFg), "min")
    if gnomadAFVal != "-" and gnomadAFgVal != "-":
        # gnomadAFVal=float(varObj.gnomadAF)
        # gnomadAFgVal=float(varObj.gnomadAFg)
        if gnomadAFVal < 0.01 and gnomadAFgVal < 0.01:
            conservationScoreGnomad = "High"
        else:
            conservationScoreGnomad = "Low"
    else:
        conservationScoreGnomad = "-"

    # DGV
    conservationScoreDGV = "-"
    if "deletion" in varObj.dgvSubtypeList or "loss" in varObj.dgvSubtypeList:
        conservationScoreDGV = "Low"
    else:
        conservationScoreDGV = "High"

    # gene O/E score
    conservationScoreOELof = "-"
    if varObj.gnomadGeneOELofUpper != "-":
        gnomadGeneOELofUpperVal = float(varObj.gnomadGeneOELofUpper)
        if gnomadGeneOELofUpperVal < 0.35:
            conservationScoreOELof = "High"
        else:
            conservationScoreOELof = "Low"

    # return
    retList = [conservationScoreGnomad, conservationScoreDGV, conservationScoreOELof]
    return retList


def _parse_numeric(valStr: str, select: str = "min"):
    """
    Function to convert string to float,
    and takes care of situation when multiple values exist
    """
    select_method = {"min": min, "max": max}
    vals = valStr.split(",")
    if "-" in vals:
        return "-"
    else:
        vals = [float(i) for i in vals]
        return select_method[select](vals)
