from .bin import create_lime_explainer, create_lime_json, create_shap_json

import joblib
import pandas as pd
import numpy as np
import shap
from lime import lime_tabular
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
import xgboost as xgb
from typing import Dict
import logging
import json
from typing import List, Optional, Union, Dict, Any


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ModelInterpreter:
    def __init__(self):                
        self.model = None
        self.model_type = None
        self.shap_explainer = None
        self.shap_values = None
        self.lime_explainer = None
        self.model_features = None
        self.is_model_loaded = False
        self.variant_ids = None
        self.feature_names = None

    def set_model(self, model) -> None:
        """
        Accept a pre-loaded model directly and extract its feature names.
        """
        try:
            # Check the model type and extract feature names
            if isinstance(model, (RandomForestClassifier, RandomForestRegressor)):
                self.model_type = "random_forest"
                self.task_type = "classification" if isinstance(model, RandomForestClassifier) else "regression"
                self.model_features = list(model.feature_names_in_)
            elif isinstance(model, xgb.XGBModel):
                self.model_type = "xgboost"
                self.task_type = "classification" if model.objective.startswith('binary:') else "regression"
                self.model_features = list(model.get_booster().feature_names)
            else:
                raise ValueError(f"Unsupported model type: {type(model)}")
            
            self.model = model
            self.is_model_loaded = True
            logger.info(f"Model set successfully with {len(self.model_features)} features")
            logger.info(f"Model features: {self.model_features}")
            
        except Exception as e:
            logger.error(f"Error setting model: {str(e)}")
            raise
        
    def load_model(self, model_path: str) -> None:
        """
        Load a model from a file using joblib.
        """
        try:
            # Load the model from a file
            loaded_model = joblib.load(model_path)
            
            # Delegate to set_model to extract features and set properties
            self.set_model(loaded_model)
        
        except Exception as e:
            logger.error(f"Error loading model: {str(e)}")
            raise
    
    def get_feature_differences(self, df: pd.DataFrame) -> Dict[str, set]:
        """
        Utility method to show feature differences between input data and model
        """
        if not self.is_model_loaded:
            raise ValueError("Model not loaded. Please load model first.")
            
        df_features = set(df.columns)
        model_features = set(self.model_features)
        
        return {
            'missing_features': model_features - df_features,
            'extra_features': df_features - model_features
        }         
      
    def prepare_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Match feature names between model and feature dataframe
        """
        try:
            if not self.is_model_loaded:
                raise ValueError("Model not loaded. Please load model first.")

            # Store variant IDs
            self.variant_ids = list(df.index)
            
            # Get current features
            df_features = df.columns
            
            # Find missing and extra features
            missing_features = set(self.model_features) - set(df_features)
            extra_features = set(df_features) - set(self.model_features)
            
            # Log the differences
            logger.info(f"Missing features: {missing_features}")
            logger.info(f"Extra features: {extra_features}")
            
            # Subset to model features
            subset_df = df[self.model_features]
            
            # Store feature names
            self.feature_names = self.model_features
            
            return subset_df
            
        except Exception as e:
            logger.error(f"Error preparing data: {str(e)}")
            raise

    def calculate_shap_values(self, df: pd.DataFrame) -> None:
        """
        Calculate SHAP values with optimization for speed
        """
        try:
            if not self.is_model_loaded:
                raise ValueError("Model not loaded. Please load model first.")

            if self.shap_explainer is None:
                # Use approximate=True for faster computation
                self.shap_explainer = shap.TreeExplainer(
                    self.model,
                    feature_perturbation='interventional',
                    approximate=True
                )
                logger.info(f"Created SHAP explainer for {self.model_type} model")
                
            self.shap_values = self.shap_explainer(df)
            logger.info("SHAP values calculated successfully")
            
            # Conditionally adjust SHAP values for RandomForestClassifier
            if self.model_type == "random_forest" and self.task_type == "classification":
                # Extract SHAP values for the positive class
                self.shap_values.values = self.shap_values.values[:, :, 1]
                logger.info("Extracted SHAP values for positive class (class 1)")

            # Set variant_ids and feature_names if not already set
            self.variant_ids = list(df.index)
            self.feature_names = list(df.columns)

        except Exception as e:
            logger.error(f"Error calculating SHAP values: {str(e)}")
            raise
    
    def run_shap_analysis(self, output_file: str = './results/shap_values.json') -> str:
        """
        Generate SHAP analysis JSON
        
        Parameters:
        -----------
        output_file : str, optional
            Path where to save the JSON file (default: './shap_values.json')
        """
        if self.shap_values is None:
            raise ValueError("SHAP values not calculated. Run calculate_shap_values first.")
    
        return create_shap_json(
            variant_ids=self.variant_ids,
            feature_names=self.feature_names,
            shap_values=self.shap_values,
            output_file=output_file  # Pass the output path
        )
    
    def run_lime_analysis(self, processed_df: pd.DataFrame, num_features = 20,
                         categorical_features: Optional[List[int]] = None, 
                         output_file: str = './results/lime_values.json') -> str:
        """
        Run LIME analysis on processed data
        
        Parameters:
        -----------
        processed_df : pd.DataFrame
            Processed feature matrix
        categorical_features : List[int], optional
            List of categorical feature indices
        output_file : str, optional
            Path where to save the JSON file
            
        Returns:
        --------
        str
            JSON string containing the LIME explanations
        """
        if not self.is_model_loaded:
            raise ValueError("Model not loaded. Please load model first.")
    
        if self.lime_explainer is None:
            self.lime_explainer = create_lime_explainer(
                X_train=processed_df.values,
                feature_names=self.feature_names,
                categorical_features=categorical_features
            )
        
        return create_lime_json(
            variant_ids=self.variant_ids,
            feature_names=self.feature_names,
            processed_df=processed_df,
            num_features = num_features,
            model=self.model,
            explainer=self.lime_explainer,
            output_file=output_file
        )