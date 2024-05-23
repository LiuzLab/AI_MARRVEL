import argparse
import os
import pandas as pd
from confidence import *
import joblib
from time import time
import generate_bivar_data
from datetime import datetime
from multiprocessing import Pool
from tqdm import tqdm
from scipy.stats import rankdata
from integrate_output import *

parser = argparse.ArgumentParser()

parser.add_argument('-id', metavar='I', type=str,
                    help = 'sample ID')

parser.add_argument('-n_cpu',  type=int, default=10,
                    help = 'folders containing all extended final matrices')

args = parser.parse_args()

#st_time = time()
#prj_name = args.project
sample_id = args.id
n_cpu = args.n_cpu

out_folder = '/out/conf_4Model'

if not os.path.exists(out_folder):
    os.mkdir(out_folder)

print(f"### Results will be saved at: {out_folder}.")
print(f"### Process started at: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}.")

model_names = ['default', 'recessive', 'nd', 'nd_recessive']
model_dict  = {}

for mn in model_names:
    model = joblib.load(f"/run/data_dependencies/model_inputs/{mn}/final_model.job")
    ref_panel = joblib.load(f"/run/data_dependencies/model_inputs/{mn}/reference_panel.job")
    features = open(f"/run/data_dependencies/model_inputs/{mn}/features.csv", 'r').readline().split(",")
    model_dict[mn] = {'model': model, 'ref': ref_panel, 'features': features}
    #if not os.path.exists(f"./{out_folder}/{mn}"):
    #    os.mkdir(f"./{out_folder}/{mn}")

print(f"### All models, features, reference panels loaded.")

print(f"### Start generating recessive matrix and running predictions.")

def assign_ranking(df):
    pred_df = df.copy()
    pred_df = pred_df.sort_values("predict", ascending=False)
    min_rankings = rankdata(1 - pred_df["predict"].to_numpy(), method='min').astype(int)
    max_rankings = rankdata(1 - pred_df["predict"].to_numpy(), method='max').astype(int)
    pred_df['min_ranking'] = min_rankings
    pred_df['max_ranking'] = max_rankings
    pred_df['ranking'] = max_rankings
    return pred_df
    

def AIM(data_folder, sample_id, n_thread):
    feature_fn = f'/out/final_matrix/{sample_id}.csv'

    if not os.path.exists(feature_fn):
        print(f"{feature_fn} does not exist.")
        return
    
    df = pd.read_csv(feature_fn,  index_col=0)    

    for mn in ['default', 'nd']:
        df_pred = df.copy()
        if 'predict' in df_pred.columns:
            df_pred = df_pred.drop(columns=['predict'])
        predict = model_dict[mn]['model'].predict_proba(df_pred.loc[:, model_dict[mn]['features']])[:, 1]
        df_pred.insert(loc=df_pred.shape[1]-1,column='predict',value=predict)
        df_pred = assign_confidence_score(model_dict[mn]['ref'], df_pred)
        df_pred = df_pred.sort_values('confidence', ascending=False)
        df_pred = assign_ranking(df_pred)
        df_pred.to_csv(f"{out_folder}/{sample_id}_{mn}_predictions.csv")
    
    print(f"### Start processing recessive data and make predictions.")

    
    default_pred = pd.read_csv(f"{out_folder}/{sample_id}_default_predictions.csv", index_col=0)

    generate_bivar_data.process_sample( data_folder = out_folder,
                                        sample_id = sample_id,
                                        default_pred = default_pred,                                           
                                        labeling=False, n_thread = n_thread)

    recessive_feature_file = f"{out_folder}/recessive_matrix/{sample_id}.csv"
    if os.path.exists(recessive_feature_file):
        # AIM found recessive variant pairs and generated the matrix
        df = pd.read_csv(recessive_feature_file,  index_col=0)
        for mn in ['recessive', 'nd_recessive']:
            df_pred = df.copy()
            if 'predict' in df_pred.columns:
                df_pred = df_pred.drop(columns=['predict'])
            predict = model_dict[mn]['model'].predict_proba(df_pred.loc[:, model_dict[mn]['features']])[:, 1]
            df_pred.insert(loc=df_pred.shape[1]-1,column='predict',value=predict)
            df_pred = assign_confidence_score(model_dict[mn]['ref'], df_pred)
            df_pred = df_pred.sort_values('confidence', ascending=False)
            df_pred = assign_ranking(df_pred)
            df_pred.to_csv(f"{out_folder}/{sample_id}_{mn}_predictions.csv")
    else:
        # AIM found no recessive variant pair
        pass

    print(f"Integrating all information for sample {sample_id}...")
    ######### Construction Paused here #########
    integrated_df = integrate_output(out_folder, data_folder, sample_id)
    if not os.path.exists(f"{out_folder}/integrated"):
        os.mkdir(f"{out_folder}/integrated")
    integrated_df.to_csv(f'{out_folder}/integrated/{sample_id}_integrated.csv')
    return

#for sample_id in tqdm(sample_folders):
AIM(out_folder, sample_id, n_thread = n_cpu)
