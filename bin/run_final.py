#!/usr/bin/env python3.8
from predict_new.utilities import *
import warnings
import json
import os
from glob import glob
from datetime import datetime
import joblib
from predict_new.confidence import *
import sys

# read in precalulated model
model = joblib.load("predict_new/hg19/final_model_wo_bg_val.job")

# read in features to use
features = list(pd.read_csv("predict_new/hg19/features.csv"))

# predict with Linhua's code
out = rank_patient(model, sys.argv[1] + ".matrix.txt", features, train_causals=None)

# write out
out.to_csv(sys.argv[1] + ".csv")
