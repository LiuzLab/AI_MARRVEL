from utilities import *
import warnings
import json
import os
from glob import glob
from datetime import datetime
import joblib
from confidence import *

def test(final_model, val_fns, features, test_label, out_folder, mode, train_causals, train_pred_folder):
    if train_pred_folder is not None:
        causal_distr = get_causal_distribution(train_pred_folder, label_col=test_label)
    else:
        causal_distr = None
        
    causal_rank_raw = rank_test_patients(final_model, val_fns, features, test_label, causal_distr = causal_distr,
                            out_folder=f"{out_folder}/{mode}_pred_raw_rank", train_causals=None)

    f_raw, eval_df_raw = plot_topK_performance(causal_rank_raw)
    f_raw.savefig(f"{out_folder}/{mode}_test_perf_raw.png", dpi=100)
    eval_df_raw.to_csv(f"{out_folder}/{mode}_test_perf_raw.csv")

    ## evaluate only on non-seen variants
#     causal_rank_exclude = rank_test_patients(final_model, val_fns, features, test_label,train_pred_folder,
#                             out_folder=f"{out_folder}/{mode}_pred_exclude_rank",
#                              train_causals=train_causals)

#     f_exclude, eval_df_exclude = plot_topK_performance(causal_rank_exclude)
#     f_exclude.savefig(f"{out_folder}/{mode}_test_perf_exclude.png", dpi=100)
#     eval_df_exclude.to_csv(f"{out_folder}/{mode}_test_perf_exclude.csv")

def main(feature_fillna=None, feature_tier_with_clinvar=None,
        train_label=None, test_label=None, out_folder=None,
        tune=True, train_test=None):
    """
    feature_fillna: bool, fill na feature or not
    feature_tier_with_clinvar: bool, tier feature use clinvar/hgmd info or not
    train_label: str, label to be used in training (is_causal, is_strong)
    test_label: str, label to be used in training (is_causal, is_strong)
    out_folder: if not specified, will be created automatically
    tune: bool, whether or not tuning is needed
    train_test: train/test split strategy, 
                    if str, must be one of ['BG_val', 'version1', 'leave_one_out']
                    if int/float, represent number/fraction of test samples to be randomly used.
    """
    ### -> load parameters
    f = open("meta.json")
    meta = json.load(f)
    data_folder = meta["data_folder"]
    classifier = meta['classifier']
    if train_test is None:
        train_test = meta['train_test']
    if feature_fillna is None:
        feature_fillna = meta["feature_fillna"]
    if feature_tier_with_clinvar is None:
        feature_tier_with_clinvar = meta["feature_tier_with_clinvar"]
    if train_label is None:
        train_label = meta["train_label"]
    if test_label is None:
        test_label = meta["test_label"]
    if out_folder is None:
        fillna = "fillNA" if feature_fillna else "withNA"
        tier = 'tierDB' if feature_tier_with_clinvar else 'tierNDB'
        dateid = datetime.now().strftime('%m-%d-%Y')
        fd = '_'.join([classifier, fillna, tier, train_label, train_test, f'tune{tune}', dateid])
        out_folder = f"./results/{fd}"

    if not os.path.exists("./results/"):
        os.mkdir("./results/")
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    
    print(f"### Classifier: {classifier}.")
    print(f"### DataHub: {data_folder}.")
    print(f"### TrainTestMode: {train_test}.")
    print(f"### Features' NA filled: {feature_fillna}.")
    print(f"### Tier used disease DB: {feature_tier_with_clinvar}.")
    print(f"### Results will be saved at: {out_folder}.")
    print(f"### Process started at: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}.")

    ### <- load parameters
    if feature_fillna and feature_tier_with_clinvar:
        ext = "fillna_tierTrue.txt"
        fns = glob(f"{data_folder}/*_{ext}")
    elif feature_fillna and (not feature_tier_with_clinvar):
        ext = "fillna_tierFalse.txt"
        fns = glob(f"{data_folder}/*_{ext}")
    elif (not feature_fillna) and feature_tier_with_clinvar:
        ext = "keepna_tierTrue.txt"
        fns = glob(f"{data_folder}/*_{ext}")
    else:
        ext = "keepna_tierFalse.txt"
        fns = glob(f"{data_folder}/*_{ext}")

    ### [Start] ---> Get training and testing files
    ## Files BG provided for additional validation are taken out of training
    fns = [os.path.abspath(fn) for fn in fns] # take absolute paths
    sample_meta = pd.read_csv(meta["sample_meta_fn"], sep=" ")
    bg_samples = sample_meta.loc[sample_meta.val==1,'BG_ID'].tolist()
    print(len(bg_samples))
    ## Testing files
    bg_val_fns = [os.path.abspath(f"{data_folder}/{vid}_{ext}") for vid in bg_samples]
    udn_val_fns = [fn for fn in fns if 'UDN' in fn]
    ddd_val_fns = [fn for fn in fns if 'DDD' in fn]
    print(f'#fns={len(fns)}, #BG validation={len(bg_val_fns)}, #UDN validation={len(udn_val_fns)}, #DDD validation={len(ddd_val_fns)}')
    ## Training files
    train_fns = list(set(fns) - set(udn_val_fns) - set(bg_val_fns) - set(ddd_val_fns))
    ### <- Get training and testing files [End]

    # if isinstance(train_sample_num, int) and (len(train_fns) > train_sample_num):
    #     train_fns = np.random.choice(train_fns, train_sample_num, replace=False)
    if tune:
        validation_performances = tune_params(train_fns,  classifier=classifier)
        validation_performances.to_csv(f"{out_folder}/param_tuning_perfs.csv")
        val_best_params = validation_performances.iloc[0, :].tolist()
        if classifier == 'xgb':
            n_est, md, classifier = int(val_best_params[0]), int(val_best_params[1]), float(val_best_params[2])
        else:
            n_est, md = int(val_best_params[0]), int(val_best_params[1])
    else:
        ### Fix parameter tuned, TO BE REMOVED --> 
        if classifier == 'xgb':
            n_est, lr, md = 300, 0.1, None
            print(f"n_est={n_est}, lr={lr}, md={md}")
        elif classifier == 'rf':
            n_est, md = 300, None
        ### <-- Fix parameter tuned, TO BE REMOVED

    if classifier == 'xgb':
        trainX, trainY, neg_pos_ratio = get_training_data(train_fns, train_label)
        features = trainX.columns.tolist()
        print(f"[Final XGB] Training with {len(train_fns)} samples.")
        final_model = train_XGBoost(trainX, trainY, n_est, lr, md, neg_pos_ratio)
        causal_idx = np.where(trainY == 1)[0]
        train_causals = trainX.iloc[causal_idx,:]
        train_causals.to_csv(f"{out_folder}/train_causal_df.csv")
        train_causals = train_causals.index.tolist()
    elif classifier == 'rf':
        final_model, features, train_causals = train_minibatch_RF(train_fns, n_estimator=n_est, max_depth=md,
         label_col=train_label, batch=2000)
        print(f"[Final RF] Trained, #trees = {final_model.get_params('n_estimators')}")

    joblib.dump(final_model, f"{out_folder}/final_model_wo_bg_val.job")
    
    with open(f"{out_folder}/features.csv", "w") as feature_f:
        feature_f.write(",".join(features))
    
    # ### predict on training samples
    # test(final_model, train_fns, features, test_label, out_folder, 'Training', train_causals, train_pred_folder=None)
    ### BG validation
    test(final_model, bg_val_fns, features, test_label, out_folder, 'BG', train_causals, train_pred_folder=None)
    ### UDN validation
    test(final_model, udn_val_fns, features, test_label, out_folder, 'UDN', train_causals, train_pred_folder=f"{out_folder}/BG_pred_raw_rank")
    ### DDD validation
    test(final_model, ddd_val_fns, features, test_label, out_folder, 'DDD', train_causals, train_pred_folder=f"{out_folder}/BG_pred_raw_rank")
    


if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main(tune=False)
    