#!/usr/bin/env python3
# bin/model_interpreter/main.py

import argparse
from variant_model_interpreter import ModelInterpreter
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Run model interpretation (SHAP/LIME) analysis')
    parser.add_argument('--analysis-type', type=str, default='both',  # Set default to 'both'
                      choices=['shap', 'lime', 'both'],
                      help='Type of analysis to run (shap, lime, or both)')
    parser.add_argument('--model-path', type=str, required=True,
                      help='Path to the trained model .job file')
    parser.add_argument('--input-path', type=str, required=True,
                      help='Path to the input CSV file')
    parser.add_argument('--output-path', type=str, required=True,
                      help='Path where to save the analysis JSON')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Initialize interpreter
    interpreter = ModelInterpreter()
    
    # Load model
    interpreter.load_model(args.model_path)
    
    # Load and prepare data
    input_df = pd.read_csv(args.input_path, index_col=0)
    processed_df = interpreter.prepare_data(input_df)
    
    # Run specified analysis
    if args.analysis_type in ['shap', 'both']:
        interpreter.calculate_shap_values(processed_df)
        shap_output = args.output_path if args.analysis_type == 'shap' else args.output_path.replace('.json', '_shap.json')
        interpreter.run_shap_analysis(output_file=shap_output)
        
    if args.analysis_type in ['lime', 'both']:
        lime_output = args.output_path if args.analysis_type == 'lime' else args.output_path.replace('.json', '_lime.json')
        interpreter.run_lime_analysis(processed_df, output_file=lime_output)

if __name__ == "__main__":
    main()