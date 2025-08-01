import pandas as pd
from scipy.stats import rankdata
from predict_new.confidence import *

def rank_patient(model, fn, features, train_causals=None):
    df = pd.read_csv(
        fn, sep="\t", index_col=0
    )  ## File name format can be a parameter in the future
    

    identifier = fn.split("/")[-1].split(".")[0]
    patient_features = df.loc[:, features]
    patient_preds = model.predict_proba(patient_features)[:, 1]
    pred_df = df.copy()
    pred_df["predict"] = patient_preds

    if train_causals is not None:
        pred_df = pred_df.loc[~pred_df.index.isin(train_causals), :]
    pred_df = pred_df.sort_values("predict", ascending=False)
    min_rankings = rankdata(1 - pred_df["predict"].to_numpy(), method="min").astype(int)
    max_rankings = rankdata(1 - pred_df["predict"].to_numpy(), method="max").astype(int)
    pred_df["min_ranking"] = min_rankings
    pred_df["max_ranking"] = max_rankings
    pred_df["ranking"] = max_rankings
    pred_df["identifier"] = identifier
    return pred_df
