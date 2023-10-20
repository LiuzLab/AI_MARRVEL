from utilities import *
import warnings
import json
import os
from glob import glob
from datetime import datetime
import joblib
from confidence import *
import sys

#read in precalulated model
model=joblib.load("/run/data_dependencies/predict_new/hg19/final_model_wo_bg_val.job")

#read in features to use
features=list(pd.read_csv("/run/data_dependencies/predict_new/hg19/features.csv"))

#predict with Linhua's code
out=rank_patient(model, "/out/"+sys.argv[1]+".matrix.txt", features, train_causals=None)

#write out
out.to_csv("/out/final_matrix/"+sys.argv[1]+".csv")