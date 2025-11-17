def getCurationScore(varObj):
    """
    Determine per-database curation scores for a variant.
    Returns list: [omimScore, hgmdScore, clinVarScore, totalScore]
    """
    # OMIM score
    if varObj.omimVarFound == 1:
        omimScore = "High" if varObj.omimSymMatchFlag == 1 else "Medium"
    elif varObj.omimGeneFound == 1:
        omimScore = "Medium" if varObj.omimSymMatchFlag == 1 else "Low"
    else:
        omimScore = "Low"

    # HGMD score
    if varObj.hgmdVarFound == 1:
        hgmdScore = "High" if varObj.hgmdSymMatchFlag == 1 else "Medium"
    elif varObj.hgmdGeneFound == 1:
        hgmdScore = "Medium" if varObj.hgmdSymMatchFlag == 1 else "Low"
    else:
        hgmdScore = "Low"

    # ClinVar score
    pathogenic_terms = {
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
        "Pathogenic, association",
        "Pathogenic, other",
        "Pathogenic, protective",
        "Pathogenic, protective, risk factor",
        "Pathogenic, risk factor",
        "Pathogenic/Likely pathogenic, risk factor",
    }
    benign_terms = {
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
    }

    clinVarScore = "Low"
    if varObj.clinVarVarFound == 1:
        if varObj.clinvarSignDesc in pathogenic_terms:
            clinVarScore = "High" if varObj.clinVarSymMatchFlag == 1 else "Medium"
        elif varObj.clinvarSignDesc in benign_terms:
            clinVarScore = "Low"
        else:
            clinVarScore = "Medium" if varObj.clinVarSymMatchFlag == 1 else "Low"
    elif varObj.clinVarGeneFound == 1:
        clinVarScore = "Medium" if varObj.clinVarSymMatchFlag == 1 else "Low"

    # Total score placeholder (requires future logic)
    totalScore = "-"
    return [omimScore, hgmdScore, clinVarScore, totalScore]
