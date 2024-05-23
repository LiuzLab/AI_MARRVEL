from os.path import exists
from scipy import stats 
from os import mkdir
import numpy as np
import pandas as pd

__author__ = "Dongxue Mao; Linhua Wang"
__email__ = "dmao@bcm.edu; linhuaw@bcm.edu"


#runs diffusion module 

def diffusion(nn, y, alpha = 0.5, max_iter = 100):
    """Helper function that runs the network diffusion.

    Parameters:
    ----------
    nn: 2D-numpy array, normalized version of the adjacency matrix of the network of interest
    y: list of float, initial heat of each node
    alpha: float, parameter for diffusion, representing the percentage of heat diffused out at each iteration
    max_iter: integer, max number of iterations of the diffusion

    Return:
    -------
    F: 2D-numpy array, diffused network
    """
    alpha = 0.5
    max_iter = 100
    y = y.to_numpy()
    F = y
    ## initialize an empty array for the heat
    Fs = np.empty((nn.shape[0], 0), float)
    ## restart
    fY = (1-alpha) * y
    for i in range(0, max_iter):
        Fs = np.append(Fs, F, axis=1)
        F = alpha * nn @ Fs[:,[i]] + fY
    return(F)

def diffuseSample(ID, Anno_df, Phrank_folder):
    """Method that runs diffusion network for a specific sample.

    Parameters:
    ----------
    ID: str, sample ID.
    Anno_folder: str, path to the folder that contains the annotation file for module 1 & 2.
    Phrank_folder: str, path to the folder that contains Phrank gene ranking score for diffusion.

    Return:
    -------
    Pandas Data Frame, index indicates variant ID, 
        columnns contatins diffused phrank scores normalized as percentile.
    """

    ## Get the PPI network for diffusion
    cor_path = "/run/data_dependencies/mod5_diffusion/combined_score.hdf5"
    cor_df = pd.read_hdf(cor_path, mode='r')
    cor_df = cor_df[(cor_df.T != 0).sum() > 1]
    cor_df = cor_df.loc[cor_df.index, cor_df.index]
    cor = cor_df.values
    cor_GeneID = pd.DataFrame({"ID":cor_df.columns.tolist()})

    Phrank_path = Phrank_folder + ID.split('.')[0] + ".txt"
    ## Validate Phrank file's existence
    if not exists(Phrank_path):
        print("File %s not exists, skipping..." %Phrank_path)
        return

    Phrank = pd.read_csv(Phrank_path, sep = "\t", header = None)
    Phrank = Phrank.rename({0:"Ensembl_Gene_ID",1:"Score"}, axis = "columns")
    ## remove rows without Ensembl_Gene_ID
    Phrank["Ensembl_Gene_ID"].replace('', np.nan, inplace=True)
    Phrank.dropna(subset=["Ensembl_Gene_ID"], inplace= True)
    ## remove duplicated genes and keep the highest similarity score
    Phrank = Phrank.sort_values('Score', ascending=False).drop_duplicates('Ensembl_Gene_ID').sort_index()
    simi = Phrank.rename({0:"Ensembl_Gene_ID","Score":"Similarity_Score"}, axis = "columns")


    ## load module 1 & 2 file as reference
    #geneFiltFn = "%s/%s.txt.gz" %(Anno_folder, ID)

    ## Validate annotation file's existence
    #if not exists(geneFiltFn):
        #print("Module 1-2 file not exists for %s, skipping..." %ID)
        #return

    #m12_df = pd.read_csv(geneFiltFn, sep="\t", compression="gzip")
    m12_df=Anno_df
    m12_genes = [g for g in m12_df.loc[:, "geneEnsId"].tolist() if "ENSG" in g]
    m12_gc = len(set(m12_genes))

    ## Filter Genes
    ## Keep all genes in filtered GeneID in module 1 and 2 and Phrank, but not in cor_GeneID
    m12_wSimi = list(set(Phrank["Ensembl_Gene_ID"]) & set(m12_genes))
    m12_wSimi_woCor = list(set(m12_wSimi) - (set(cor_GeneID["ID"])))
    m12_wSimi_woCor = simi[simi["Ensembl_Gene_ID"].isin(m12_wSimi_woCor)]
    ## normalized the cor matrix
    net = abs(cor)
    D = 1/np.sqrt(net.sum(axis = 1))
    D2 = np.diag(D)
    net_norm = D2 @ net @ D2

    ## Set Y as similarity score
    Y = cor_GeneID.merge(simi, left_on = "ID", right_on = "Ensembl_Gene_ID", how = "left")
    Y = Y[["ID","Similarity_Score"]].fillna(0)
    Y.set_index("ID", inplace = True)

    ## Start diffusion
    diff_res = diffusion(net_norm,Y,0.5,100)

    ## Combine diffusion result with GeneID
    diff_wGeneID = pd.DataFrame(diff_res, columns=["Final_Heat"])
    diff_wGeneID["GeneID"] = cor_GeneID

    ## remove the genes with zero heat
    diff_wGeneID = diff_wGeneID[diff_wGeneID["Final_Heat"] != 0]
    diff_wGeneID = diff_wGeneID[["GeneID","Final_Heat"]]

    ## combine diffusion result with similarity scores
    FinalHeat_wSimi = pd.concat([m12_wSimi_woCor,diff_wGeneID],ignore_index=True, sort=True)

    ## remove diffusion result with non-vep genes
    diff_wGeneID_in_m12 = diff_wGeneID[diff_wGeneID['GeneID'].isin(m12_genes)]
    m12_wSimi_woCor = m12_wSimi_woCor.rename(columns = {"Similarity_Score":"Final_Heat"})
    diff_wGeneID_in_m12 = diff_wGeneID_in_m12.rename(columns = {"GeneID":"Ensembl_Gene_ID"})
    FinalHeat_wSimi = pd.concat([m12_wSimi_woCor,diff_wGeneID_in_m12],ignore_index=True, sort=True)
    FinalHeat_wSimi = FinalHeat_wSimi.sort_values(by = 'Final_Heat', ascending = False)

    ## Get the diffusion scores for genes in the annotation data, 0 if not existed
    score_ordered = []
    gene_ordered = m12_df.loc[:, "geneEnsId"].tolist()
    for gene in gene_ordered:
        if gene in FinalHeat_wSimi.Ensembl_Gene_ID.tolist():
            score_ordered.append(FinalHeat_wSimi.loc[FinalHeat_wSimi.Ensembl_Gene_ID == gene, 'Final_Heat'].to_numpy()[0])
        else:
            score_ordered.append(0)

    ## Normalize to have a z-score (max 1, min 0)
    scores = stats.rankdata(score_ordered, 'max')/len(score_ordered)

    m12_df['diffuse_Phrank_STRING'] = scores
    m12_df = m12_df.drop_duplicates(subset = ["varId_dash"])
    m12_df = m12_df[["varId_dash", "diffuse_Phrank_STRING"]]
    m12_df = m12_df.set_index("varId_dash", drop=True)
    return m12_df


