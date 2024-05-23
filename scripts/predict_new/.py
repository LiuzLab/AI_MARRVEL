from utilities import *
import warnings
import json
import os
from glob import glob
from datetime import datetime
import joblib
from confidence import *

model=joblib.load("/home/lul0614/pipeline/predict_new/final_model_wo_bg_val.job")
features=pd.read_csv("features.csv").columns
out=rank_patient(model, "/local/LucianLi/final-all/989068_fillna_tierFalse.txt", features, train_causals=None)

out.to_csv("test.csv")