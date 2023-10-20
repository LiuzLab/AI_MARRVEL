from math import ceil
import os
import pandas as pd
import numpy as np
from scipy.stats import rankdata
import sys
import xgboost as xgb
from tqdm import trange
from time import time
import json
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier as RF
import itertools
from confidence import *

def get_training_data(fns, label_col):
    """Read in training data, drop features with specified mode in the meta['feature_mode'] json.
        Allowed modes are ['all', 'raw', 'var_disase_db_off', 'gene_disase_db_off']
    """
    np.random.seed(0)
    ### READ PARAMETERS
    f = open("meta.json")
    meta = json.load(f)
    feature_mode = meta["feature_mode"]
    feature_meta = pd.read_csv(meta["feature_meta_fn"], sep=',')
    ### READ training samples
    t0 = time()
    print("[TrainVal-Train] Start reading training data.")
    dfs = []
    for i in trange(len(fns)):
        fn = fns[i]
        df = pd.read_csv(fn, index_col=0, sep="\t")
        dfs.append(df)
    dfs = pd.concat(dfs)

    t1 = time()
    print(f"[TrainVal-Train] Completed reading training data in {(t1-t0):.2f} seconds.")
    train_Y = dfs[label_col].to_numpy().astype(int)
    neg_pos_ratios = [(train_Y == 0).sum()/(train_Y == 1).sum()]

    ## Drop features based on feature mode
    if feature_mode == "all":
        features = feature_meta.Feature.tolist()
    elif feature_mode == "raw":
        features = feature_meta.loc[feature_meta.Engineered != 1,:].Feature.tolist()
    elif feature_mode == "var_disase_db_off":
        features = feature_meta.loc[feature_meta.DiseaseDB_VariantLevel != 1,:].Feature.tolist()
    elif feature_mode == "gene_disase_db_off":
        features = feature_meta.loc[feature_meta.DiseaseDB_GeneLevel != 1,:].Feature.tolist()
    else:
        print(f"Wrong feature mode entered, should be one of\
        ['all', 'raw', 'var_disease_db_off', 'gene_disease_db_off'],\
        your input is {feature_mode}."
        )
        sys.exit(0)
    train_X = dfs[features]
    print("[TrainVal-Train] Training data obtained.")
    return train_X, train_Y, np.mean(neg_pos_ratios)

def tune_params(fns, classifier='xgb'):
    ### READ PARAMETERS
    f = open("meta.json")
    meta = json.load(f)
    train_label = meta["train_label"]
    test_label = meta["test_label"]
    train_val_frac = meta["train_val_frac"]
    sep2 =  "\t"

    train_fns_sampled = list(np.random.choice(fns, int(train_val_frac * len(fns))))
    val_fns_sampled = list(set(fns) - set(train_fns_sampled))

    print(f"[TrainVal-Train] train {len(train_fns_sampled)} samples with label {train_label},\
         validate {len(val_fns_sampled)} samples with label {test_label}.")

    para_train_dfs = [pd.read_csv(train_fn, usecols=[train_label], sep=sep2) for train_fn in train_fns_sampled]
    para_train_dfs = pd.concat(para_train_dfs)
    if classifier == 'xgb':
        ## tune parameters
        n_estimators = [200, 300]
        learning_rates = list(np.logspace(0.001, 0.2, num=10, endpoint=True))
        max_depths = [None]
        params = list(itertools.product(n_estimators, learning_rates, max_depths))
    elif classifier == 'rf':
        n_estimators = [100, 200, 300, 400]
        max_depths = [None]
        params = list(itertools.product(n_estimators, max_depths))
    
    train_X, train_Y, neg_pos_ratio = get_training_data(train_fns_sampled, train_label)
    perf_dfs = []
    for param in params:
        t0 = time()
        if classifier == 'xgb':
            n_est, lr, md = param
            print(f"[TrainVal-Train] n_trees: {n_est}, lr: {lr}, max_depth: {md}, n/p ratio: {neg_pos_ratio}.")
            clf = train_XGBoost(train_X, train_Y, n_est, lr, md, neg_pos_ratio)
            
        elif classifier == 'rf':
            n_est, md = param
            print(f"[TrainVal-Train] n_trees: {n_est},  max_depth: {md}, n/p ratio: {neg_pos_ratio}.")
            clf, _, _ = train_minibatch_RF(train_X,train_Y, n_est, md, neg_pos_ratio)
        t1 = time()
        print(f"[TrainVal-Train] done in {(t1-t0):.2f} seconds.")
        rank_dfs = rank_test_patients(clf, val_fns_sampled, train_X.columns.tolist(), test_label)
        obj_df = objective(rank_dfs)
        if classifier == 'xgb':
            perf_df = pd.DataFrame(data=[[n_est, md,lr, obj_df['objective'],
                                        obj_df['var_top1_acc'], obj_df['var_top5_acc'], obj_df['var_top10_acc'] ]],
                                columns=['n_estimator', 'max_depth', 'learning_rate',
                                'objective', 'var_top1', 'var_top5', 'var_top10'])
        else:
            perf_df = pd.DataFrame(data=[[n_est, md, obj_df['objective'],
                                        obj_df['var_top1_acc'], obj_df['var_top5_acc'], obj_df['var_top10_acc'] ]],
                                columns=['n_estimator', 'max_depth',
                                'objective', 'var_top1', 'var_top5', 'var_top10'])
        perf_dfs.append(perf_df)
        t2 = time()
        print(f"[TrainVal-Validate] done in {(t2-t1):.2f} seconds, objective = {obj_df['objective']}.")           
    perf_dfs = pd.concat(perf_dfs).sort_values(["objective", "var_top1", "var_top5", "var_top10"], ascending=False)
    return perf_dfs

def quantify_perf(test_perfs, topK=5, mode='variant_level'):
    if mode == 'sample_level':
        sample_perfs = test_perfs.groupby('identifier')[['ranking']].min()
    else:
        sample_perfs = test_perfs.copy()
        
    topK_acc = (sample_perfs['ranking'] <= topK).sum() / sample_perfs.shape[0]
    return topK_acc

def plot_topK_performance(rank_dfs):
    f = plt.figure(figsize=(10,4))
    topks = range(1, 101)
    var_accs = []
    sample_accs = []

    for i in topks:
        var_accs.append(quantify_perf(rank_dfs, topK=i))
        sample_accs.append(quantify_perf(rank_dfs, topK=i, mode='sample_level'))
    df = pd.DataFrame({"top_k": topks, "var_accuracy": var_accs, "sample_accuracy": sample_accs})
    ax1 = plt.subplot(1,2,1)
    ax1.plot(topks, var_accs, linestyle='--', marker='o', color='blue')
    ax1.set_xlabel("Top K", fontsize=14)
    ax1.set_ylabel("Variant Accuracy", fontsize=14)
    ax2 = plt.subplot(1,2,2)
    ax2.plot(topks, sample_accs, linestyle='--', marker='o', color='green')
    ax2.set_xlabel("Top K", fontsize=14)
    ax2.set_ylabel("Sample Accuracy", fontsize=14)
    plt.subplots_adjust(wspace=0.3)
    return f, df

def objective(rank_dfs):
    ### READ PARAMETERS
    f = open("meta.json")
    meta = json.load(f)
    alpha = meta["objective_params"]["alpha"]
    beta = meta["objective_params"]["beta"]
    gamma = meta["objective_params"]["gamma"]
    f.close()
    ### Get accuracy and calculate objective
    vl1 = quantify_perf(rank_dfs, topK=1)
    vl5 = quantify_perf(rank_dfs, topK=5)
    vl10 = quantify_perf(rank_dfs, topK=10)
    sl1 = quantify_perf(rank_dfs, topK=1, mode='sample_level')
    sl5 = quantify_perf(rank_dfs, topK=5, mode='sample_level')
    sl10 = quantify_perf(rank_dfs, topK=10, mode='sample_level')
    obj_score = alpha * vl1 + beta * (vl5 - vl1) + gamma * (vl10 - vl5)
    result = {"var_top1_acc": vl1,
    "var_top5_acc": vl5,
    "var_top10_acc": vl10,
    "sample_top1_acc": sl1,
    "sample_top5_acc": sl5,
    "sample_top10_acc": sl10,
    "objective": obj_score}
    return result

def train_XGBoost(X, Y, n_estimator, lr, max_depth, ratio):
    xgb.config_context(verbosity=0)
    bst = xgb.XGBClassifier(max_depth = max_depth,
                            n_estimators=n_estimator,
                            learning_rate= lr, 
                            scale_pos_weight = ratio,
                            verbosity=0, n_jobs=-1)
    bst.fit(X, Y)
    return bst

def train_RF(X, Y, n_estimator=None, max_depth=None, ratio=None, model=None):
    if model is None:
        model = RF(class_weight = {0:1, 1:ratio}, random_state=2021, n_jobs=-1, 
                        n_estimators=n_estimator, max_depth=max_depth, warm_start=True)
    model.fit(X, Y)
    return model

def train_minibatch_RF(fns, n_estimator, max_depth, label_col, batch=500):
    k = ceil(len(fns) / batch)
    model = RF(random_state=2021, n_jobs=-1, 
                        n_estimators=n_estimator,
                         max_depth=max_depth, warm_start=True)
    train_causals = []
    for i in range(k):
        if i < k-1:
            batch_fns = fns[i*batch:(i+1)*batch]
        else:
            batch_fns = fns[i*batch:]
        trainX, trainY, neg_pos_ratio = get_training_data(batch_fns, label_col=label_col)
        features = trainX.columns.tolist()
        model.set_params(**{'class_weight': {0:1, 1:neg_pos_ratio},
                            'n_estimators': n_estimator})
        model.fit(trainX, trainY)
        n_estimator += 10
        train_causals += trainX.index[np.where(trainY==1)[0]].tolist()
    return model, features, train_causals


def rank_patient(model, fn, features, train_causals=None):
    df = pd.read_csv(fn, sep="\t", index_col=0) ## File name format can be a parameter in the future
    identifier = fn.split("/")[-1].split(".")[0]
    patient_features = df.loc[:, features]
    patient_preds = model.predict_proba(patient_features)[:, 1]
    pred_df = df.copy()
    pred_df["predict"] = patient_preds

    if train_causals is not None:
        pred_df = pred_df.loc[~pred_df.index.isin(train_causals),:]
    pred_df = pred_df.sort_values("predict", ascending=False)
    min_rankings = rankdata(1 - pred_df["predict"].to_numpy(), method='min').astype(int)
    max_rankings = rankdata(1 - pred_df["predict"].to_numpy(), method='max').astype(int)
    pred_df['min_ranking'] = min_rankings
    pred_df['max_ranking'] = max_rankings
    pred_df['ranking'] = max_rankings
    pred_df['identifier'] = identifier
    return pred_df

def rank_test_patients(model, fns, features, label_col, causal_distr = None,
                       out_folder=None, train_causals=None):
    print(out_folder)
    causal_rank_dfs = []
    for i in trange(len(fns)):
        fn = fns[i]

        if not os.path.exists(fn):
            continue

        identifier = fn.split("/")[-1].split("_")[0]
        rank_df = rank_patient(model, fn, features, train_causals)
        
        if causal_distr is not None:
            rank_df =  assign_confidence_score(causal_distr, rank_df, label_col='is_strong')
        
        if out_folder is not None:
            if not os.path.exists(out_folder):
                os.mkdir(out_folder)

            rank_df.to_csv(f"{out_folder}/{identifier}_variant_rank.csv")
        
        causal_rank = get_causal_rank(rank_df, label_col)
        
        if causal_rank is not None:
            causal_rank_dfs.append(causal_rank)
            
    causal_rank_dfs = pd.concat(causal_rank_dfs)
    return causal_rank_dfs

def get_causal_rank(pred_df, label_col):
    if (pred_df[label_col] == 1).sum() == 0:
        return None
    return pred_df.loc[pred_df[label_col] == 1,:]