import joblib
import pandas as pd
import numpy as np
import shap
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
import xgboost as xgb
import json

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from utils.numpy_to_python import convert_numpy_to_python

def create_shap_json(variant_ids, feature_names, shap_values, output_file = "./results/shap_values.json"):
    """ 
    Create JSON File from shap values generated from shap.TreeExplainer
    Example Usage where processed_df is feature matrix:
        variant_ids = list(processed_df.index)
        feature_names = list(processed_df.columns)
        json_output = create_shap_json(variant_ids, feature_names, shap_values)
    """
    json_data = []
    
    for idx in range(len(shap_values.values)):
        entry = {
            "variant_id": variant_ids[idx],
            "base_value": float(shap_values.base_values[idx]),
            "model_output_score": {},
            "feature_values": {}
        }
        
        # Add feature contributions and values
        for feat_idx, feat_name in enumerate(feature_names):
            entry["model_output_score"][feat_name] = float(shap_values.values[idx][feat_idx])
            entry["feature_values"][feat_name] = float(shap_values.data[idx][feat_idx])
            
        json_data.append(entry)
    
    # Convert to JSON string with proper formatting
    json_string = json.dumps(json_data, indent=2, default=convert_numpy_to_python)
    
    # Optionally save to file
    with open(output_file, 'w') as f:
        f.write(json_string)
    
    return json_string