#!/usr/bin/env python3.8

import pandas as pd
import sys

###merge score files with final rank matrix

#read rank matrix and score file
feature=pd.read_csv("./"+sys.argv[1]+".trio.csv")
ndg = pd.read_csv("./"+sys.argv[1]+".trio.NDG.csv")
ndg = ndg.loc[:,["Unnamed: 0","predict","ranking"]]
ndg.columns = ["Unnamed: 0","predict (nd)","ranking (nd)"]

feature = feature.merge(ndg, on="Unnamed: 0", how="left")
feature.to_csv("./"+sys.argv[1]+".trio.prediction.csv",index=False)