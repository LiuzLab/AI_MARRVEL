def getClinVarUsingMarrvelFlatFile(varObj, clinvarAlleleDf, clinvarGeneDf):
    """
    function to get clinvar info using marrrvel flat file
    Params:a varaint object read from VEP annotation

    Returns:
    List of clinvar annotations
    """
    # print('in clinvar using flatfile')
    # print('\tvar:',varObj.varId_dash)
    varFound = 0
    varDict = {}
    geneFound = 0
    geneDict = {}
    clinvarTotalNumVars = 0
    clinvarNumP = 0
    clinvarNumLP = 0
    clinvarNumLB = 0
    clinvarNumB = 0
    ## NOTE(CL): no use but kept
    clinvarTitle = ""
    clinvarSignDesc = ""
    clinvarCondition = ""

    # print('in function type clinavrDf:', type(clinvarAlleleDf) )
    chromVal = int(varObj.chrom)
    posVal = int(varObj.pos)

    """ # NOTE(CL): Old implementation
    #using inetger column not index
    if 1:
        vals=clinvarAlleleDf.loc[(clinvarAlleleDf['chr']==chromVal ) &
                             (clinvarAlleleDf['start'] == posVal ) &
                             (clinvarAlleleDf['stop'] == posVal )]
        numRows=len(vals.index)
        print('columns numRows:', numRows)
    #using index
    if 0:
        idVal=str(chromVal)+'_'+str(posVal)+'_'+str(posVal)
        #vals=clinvarAlleleDf.loc[(clinvarAlleleDf['id1']==idVal)]
        try:
            vals=clinvarAlleleDf.loc[idVal]
            numRows=len(vals.index)
            print('index numRows:', numRows)
        except:
            numRows=0

    print('\tcheck var numRows:', numRows)
    if numRows > 0:
        varFound=1
        print('\tclinvar var found')
        print('\t clinvar vals:', vals)
        clinvarTitle=vals.iloc[0]['title']
        clinvarSignDesc=vals.iloc[0]['significance.description']
        clinvarCondition=vals.iloc[0]['condition']
        print('\ttitle:', clinvarTitle, 'signDes:', clinvarSignDesc, 'condition:', clinvarCondition)
    else:
        varFound=0
    """

    # NOTE(CL): check if var annotated in clinvar using VEP clinvar.vcf.gz custom annotation
    if varObj.clinvar_AlleleID != "-":
        varFound = 1
        # print('\tclinvar var found')
    else:
        varFound = 0

    # check if gene annotated in clinvar using clinvarGeneDf
    geneSymbol = varObj.geneSymbol
    # print('\tgene symbol:', geneSymbol)
    try:
        vals = clinvarGeneDf.loc[geneSymbol]
        numRows = len(vals.index)
    except:
        numRows = 0
    # print('\tcheck gene numRows:', numRows)
    if numRows > 0:
        # print('\tclinvar gene found')
        geneFound = 1
        clinvarTotalNumVars = vals["totalClinvarVars"]
        clinvarNumP = vals["P"]
        clinvarNumLP = vals["LP"]
        clinvarNumLB = vals["LB"]
        clinvarNumB = vals["B"]
    # return
    retList = [
        varFound,
        varDict,
        geneFound,
        geneDict,
        clinvarTotalNumVars,
        clinvarNumP,
        clinvarNumLP,
        clinvarNumLB,
        clinvarNumB,
        clinvarTitle,
        clinvarSignDesc,
        clinvarCondition,
    ]
    return retList


def getHGMDUsingFlatFile(varObj, hgmdHPOScoreGeneSortedDf):
    """
    function to get HGMD from local flat file
    Params:
    varObj:a varaint object read from VEP annotation
    hgmdHPOScoreGeneSortedDf: HGMD data frame read from local file

    Returns:
    List of HGMD annotations
    """
    hgmdGeneFound = 0
    hgmdVarFound = 0
    hgmdVarPhenIdList = []
    hgmdVarHPOIdList = []
    hgmdVarHPOStrList = []

    # NOTE(CL): check VarFound
    if varObj.hgmd_id != "-":
        hgmdVarFound = 1
    else:
        hgmdVarFound = 0

    # NOTE(CL): check geneFound
    if varObj.geneSymbol in hgmdHPOScoreGeneSortedDf.index:
        hgmdGeneFound = 1
    else:
        hgmdGeneFound = 0

    retList = [
        hgmdVarFound,
        hgmdGeneFound,
        hgmdVarPhenIdList,
        hgmdVarHPOIdList,
        hgmdVarHPOStrList,
    ]
    return retList

