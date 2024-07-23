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
model=joblib.load("/run/predict_trio/rf_trio_base.job")

#read in features to use
features=list(pd.read_csv("/run/predict_trio/features.csv"))

#predict with Linhua's code
out=rank_patient(model, "/out/inheritance/"+sys.argv[1]+".trio.mtx.csv", features, train_causals=None)

#write out
out.to_csv("/out/final_matrix/"+sys.argv[1]+".trio.csv")
