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

jobfile_path = sys.argv[1]
featurenames_path = sys.argv[2]
matrix_path = sys.argv[3]
output_path = sys.argv[4]

# read in precalulated model
model = joblib.load(jobfile_path)

# read in features to use
features = list(pd.read_csv(featurenames_path))

# predict with Linhua's code
out = rank_patient(model, matrix_path, features, train_causals=None)

# write out
out.to_csv(output_path)
