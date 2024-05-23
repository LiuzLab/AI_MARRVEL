#!/usr/bin/env python
"""Script to estimate the confidence of whether a case is solved or not and assign 
a confidence score to each prediction.
"""
__author__ = "Linhua Wang"
__email__ = "linhuaw@bcm.edu"
import pandas as pd
import os
import numpy as np
import scipy

def get_causal_distribution(pred_fd, label_col='is_strong'):
    # pred_fd: folder for prediction files only
    fns = os.listdir(pred_fd)
    # get top20 predictions
    top20_preds = []
    for fn in fns:
        top20_preds.append(pd.read_csv(f'{pred_fd}/{fn}', index_col=0).iloc[:20,:])
    top20_preds = pd.concat(top20_preds)
    # get scores of causal variants
    top20_causals = np.ravel(top20_preds.loc[top20_preds[label_col]==1, 'predict'].to_numpy())
    return top20_causals

def assign_confidence_score(causal_distr, pred_df):
    # pred_fd: folder for prediction files only
    # pred_df: prediction data frame of the sample of interest
    
    prediction_scores = pred_df.predict.to_numpy()

    # assign confidence score using percentile
    conf_scores = [scipy.stats.percentileofscore(causal_distr, i) for i in prediction_scores]
    # assign confidence category by discretizing
    conf_categories = []
    for conf_score in conf_scores:
        if conf_score < 25:
            conf_categories.append("Unsolved")
        elif conf_score < 50:
            conf_categories.append("Solved (Low)")
        elif conf_score < 75:
            conf_categories.append("Solved (Medium)")
        else:
            conf_categories.append("Solved (High)")
    pred_df['confidence'] = conf_scores
    pred_df['confidence level'] = conf_categories
    pred_df.sort_values("confidence", ascending=False)
    return pred_df

# def determine_case_confidence(causal_distr, pred_df, label_col='is_strong'):
#     pred_df = assign_confidence_score(causal_distr, pred_df, label_col)
#     return pred_df['confidence level'].tolist()[0]