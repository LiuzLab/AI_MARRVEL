from utilities import *
import warnings
import json
import os
from glob import glob
from datetime import datetime
import joblib
from confidence import *


causal_distr = get_causal_distribution('/houston_30t/alexw/MARRVEL_ML/results/rf_fillNA_tierNDB_is_causal_version1_tuneFalse_05-31-2022/BG_pred_raw_rank', 
                                       label_col='is_strong')

pd.DataFrame(causal_distr).to_csv('causal_distr.txt',index=False,header=False)