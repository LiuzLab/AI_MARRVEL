#!/usr/bin/env python3.8
import numpy as np
import pandas as pd

def main():
    ## Get the PPI network for diffusion
    cor_path = "mod5_diffusion/combined_score.hdf5"
    cor_df = pd.read_hdf(cor_path, mode="r")
    cor_df = cor_df[(cor_df.T != 0).sum() > 1]
    cor_df = cor_df.loc[cor_df.index, cor_df.index]

    cor = cor_df.values
    cor_GeneID_arr = pd.DataFrame({"ID": cor_df.columns.tolist()}).values

    net = abs(cor).astype('float32')
    D = 1 / np.sqrt(net.sum(axis=1))
    D2 = np.diag(D)
    net_norm = D2 @ net @ D2
    np.savez("net_norm_cor_GeneID.npz", net_norm=net_norm, cor_GeneID_arr=cor_GeneID_arr)

if __name__ == "__main__":
    main()
