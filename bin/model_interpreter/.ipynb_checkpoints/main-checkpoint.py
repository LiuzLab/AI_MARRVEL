# modules/model_interpreter/run_shap.py
import argparse
import pandas as pd
from variant_model_interpreter import ModelInterpreter

def parse_args():
    parser = argparse.ArgumentParser(description='Run SHAP analysis on model predictions')
    parser.add_argument('--model-path', type=str, required=True,
                        help='Path to the trained model .joblib file (XGB preferred, RF is slow)')
    parser.add_argument('--input-path', type=str, required=True,
                        help='Path to the input feature matrix file')
    parser.add_argument('--output-path', type=str, required=True,
                        help='Path where to save the SHAP analysis JSON')
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Initialize interpreter
    interpreter = ModelInterpreter()
    
    # Load model
    interpreter.load_model(args.model_path)
    
    # Load and prepare input data
    input_df = pd.read_csv(args.input_path, index_col=0)
    processed_df = interpreter.prepare_data(input_df)
    
    # Calculate SHAP values
    interpreter.calculate_shap_values(processed_df)
    
    # Run SHAP analysis and save results
    shap_json = interpreter.run_shap_analysis(output_file=args.output_path)

if __name__ == "__main__":
    main()