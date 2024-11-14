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

def create_shap_json(variant_ids, feature_names, shap_values, output_file="./results/shap_values.json"):
    """ 
    Create JSON file from SHAP values.
    """
    if not variant_ids or not feature_names:
        raise ValueError("variant_ids or feature_names cannot be None")

    json_data = []

    # Determine if base_values is an array and extract the relevant value for binary classification
    is_multiclass = isinstance(shap_values.base_values[0], (list, np.ndarray))

    for idx in range(len(shap_values.values)):
        # Handle the case where base_values is an array (e.g., binary classification)
        if is_multiclass:
            base_value = float(shap_values.base_values[idx][1])  # Use the base value for class 1
        else:
            base_value = float(shap_values.base_values[idx])

        entry = {
            "variant_id": variant_ids[idx],
            "base_value": base_value,
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