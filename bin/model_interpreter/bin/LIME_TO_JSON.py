# Import necessary libraries
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from utils.numpy_to_python import convert_numpy_to_python

from lime import lime_tabular
import numpy as np
import json

def create_lime_explainer(X_train, feature_names, categorical_features=None, class_names=['negative', 'positive'], mode="classification"):
    """
    Create a LIME explainer object
    """
    if categorical_features is None:
        categorical_features = []
    
    explainer = lime_tabular.LimeTabularExplainer(
        training_data=np.array(X_train),
        feature_names=feature_names,
        class_names=class_names,
        categorical_features=categorical_features,
        mode=mode,
        verbose=False,
        random_state=42
    )
    
    return explainer

def create_lime_json(variant_ids, feature_names, processed_df, model, explainer, output_file="./results/lime_values.json", num_features=20):
    """
    Create JSON File from LIME explanations with preserved exact values
    """
    json_data = []
    
    # Get explanations for all instances
    for idx, instance in enumerate(processed_df.values):
        # Get explanation for this instance
        explanation = explainer.explain_instance(
            instance, 
            model.predict_proba,
            num_features=num_features
        )
        
        # Get the feature explanations as list of (feature, importance) tuples
        exp_list = explanation.as_list()
        
        # Handle intercept value (it's a dict with class label as key)
        intercept_value = explanation.intercept[1] if isinstance(explanation.intercept, dict) else explanation.intercept
        
        # Handle local prediction (it's a numpy array)
        local_pred = explanation.local_pred[0] if isinstance(explanation.local_pred, (list, np.ndarray)) else explanation.local_pred
        
        # Create entry for this instance
        entry = {
            "variant_id": variant_ids[idx],
            "base_value": float(intercept_value),
            "confidence": float(explanation.score),
            "prediction_local": float(local_pred),
            "model_output_score": {},
            "feature_values": {},
            "feature_explanations": []  # Add this to store the original explanation format
        }
        
        # Store the original explanation format
        for feature_text, importance in exp_list:
            entry["feature_explanations"].append({
                "feature": feature_text,  # This preserves the full text with conditions
                "importance": float(importance)
            })
            
            # Extract just the feature name for the model_output_score
            # Assuming format is "feature_name > value" or similar
            feature_name = feature_text.split(' >')[0].strip()
            entry["model_output_score"][feature_name] = float(importance)
            
        # Add actual feature values
        for feature_name in feature_names:
            entry["feature_values"][feature_name] = float(instance[feature_names.index(feature_name)])
        
        json_data.append(entry)
    
    # Convert to JSON string with proper formatting
    json_string = json.dumps(json_data, indent=2, default=convert_numpy_to_python)
    
    # Save to file
    with open(output_file, 'w') as f:
        f.write(json_string)
    
    return json_string