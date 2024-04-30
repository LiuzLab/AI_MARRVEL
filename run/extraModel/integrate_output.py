import pandas as pd
from glob import glob
import os

def process_recessive_matrix(df):
    df['varPair'] = df.index.copy()
    df['var1'] = df.varPair.apply(lambda x: x.split("_")[0])
    df['var2'] = df.varPair.apply(lambda x: x.split("_")[1])
    var_grps = df.groupby("var1")

    result = []
    for _, var_grp in var_grps:
        var_grp = var_grp.sort_values("predict", ascending=False)
        result.append(var_grp.iloc[[0],:])
    result = pd.concat(result)
    result.set_index("var1", inplace=True)
    return result

def integrate_output(prj_folder, data_folder, sample_id):
    ## Group variants by Ens IDs
    expanded_fn = f"/out/final_matrix_expanded/{sample_id}.expanded.csv.gz"
    #expanded_df = []
   # for expanded_fn in expanded_fns:
    if 'csv' in expanded_fn:
        expanded_df = pd.read_csv(expanded_fn, sep=",", index_col=0, compression='infer')
    else:
        expanded_df = pd.read_csv(expanded_fn, sep="\t", index_col=0, compression='infer')
    #expanded_df.append(df)
    #expanded_df = pd.concat(expanded_df)
    # expanded_df = expanded_df.loc[~expanded_df.index.duplicated(keep='first')]
    ## read predictions for default model
    default_fn = f"{data_folder}/{sample_id}_default_predictions.csv"
    default_pred = pd.read_csv(default_fn, index_col=0)  
    default_raw_pred_dict = dict(zip(default_pred.index.tolist(), default_pred.predict.tolist()))
    default_conf_dict = dict(zip(default_pred.index.tolist(), default_pred['confidence'].tolist()))
    default_conf_lvl_dict = dict(zip(default_pred.index.tolist(), default_pred['confidence level'].tolist()))
    default_rank_dict = dict(zip(default_pred.index.tolist(), default_pred.ranking.tolist()))

    expanded_df['predict (default)'] = expanded_df.index.map(default_raw_pred_dict)
    expanded_df['confidence (default)'] = expanded_df.index.map(default_conf_dict)
    expanded_df['confidence level (default)'] = expanded_df.index.map(default_conf_lvl_dict)
    expanded_df['ranking (default)'] = expanded_df.index.map(default_rank_dict)

    ## read predictions for novel disease model
    nd_fn = f"{data_folder}/{sample_id}_nd_predictions.csv"
    nd_pred = pd.read_csv(nd_fn, index_col=0)  
    nd_raw_pred_dict = dict(zip(nd_pred.index.tolist(), nd_pred.predict.tolist()))
    nd_conf_dict = dict(zip(nd_pred.index.tolist(), nd_pred['confidence'].tolist()))
    nd_rank_dict = dict(zip(nd_pred.index.tolist(), nd_pred.ranking.tolist()))
    nd_conf_lvl_dict = dict(zip(nd_pred.index.tolist(), nd_pred['confidence level'].tolist()))

    expanded_df['predict (nd)'] = expanded_df.index.map(nd_raw_pred_dict)
    expanded_df['confidence (nd)'] = expanded_df.index.map(nd_conf_dict)
    expanded_df['ranking (nd)'] = expanded_df.index.map(nd_rank_dict)
    expanded_df['confidence level (nd)'] = expanded_df.index.map(nd_conf_lvl_dict)

    ## recessive predictions for default model
    recessive_fn = f"{data_folder}/{sample_id}_recessive_predictions.csv"
    if os.path.exists(recessive_fn):
        recessive_pred = pd.read_csv(recessive_fn, index_col=0)
        recessive_pred = process_recessive_matrix(recessive_pred)

        recessive_pred_dict = dict(zip(recessive_pred.index.tolist(), recessive_pred.predict.tolist()))
        recessive_conf_dict = dict(zip(recessive_pred.index.tolist(), recessive_pred['confidence'].tolist()))
        recessive_rank_dict = dict(zip(recessive_pred.index.tolist(), recessive_pred.ranking.tolist()))
        recessive_other_var = dict(zip(recessive_pred.index.tolist(), recessive_pred.var2.tolist()))
        recessive_conf_lvl_dict = dict(zip(recessive_pred.index.tolist(), recessive_pred['confidence level'].tolist()))

        rec_preds, rec_confs, rec_ranks, rec_var2s = [], [], [], []
        rec_conf_lvls = []
        for var in expanded_df.index.tolist():
            if var in recessive_pred_dict:
                rec_preds.append(recessive_pred_dict[var])
                rec_confs.append(recessive_conf_dict[var])
                rec_ranks.append(recessive_rank_dict[var])
                rec_var2s.append(recessive_other_var[var])
                rec_conf_lvls.append(recessive_conf_lvl_dict[var])
            else:
                rec_preds.append(-1)
                rec_confs.append(-1)
                rec_ranks.append(99999)
                rec_var2s.append('NA')
                rec_conf_lvls.append('Unsolved')

        expanded_df['predict (recessive)'] = rec_preds
        expanded_df['confidence (recessive)'] = rec_confs
        expanded_df['ranking (recessive)'] = rec_ranks
        expanded_df['recessive var2'] = rec_var2s
        expanded_df['confidence level (recessive)'] = rec_conf_lvls
    else:
        expanded_df['predict (recessive)'] = -1
        expanded_df['confidence (recessive)'] = -1
        expanded_df['ranking (recessive)'] = 99999
        expanded_df['recessive var2'] = 'NA'
        expanded_df['confidence level (recessive)'] = 'Unsolved'

    ## novel disease recessive 
    nd_recessive_fn = f"{data_folder}/{sample_id}_nd_recessive_predictions.csv"
    if os.path.exists(nd_recessive_fn):
        nd_recessive_pred = pd.read_csv(nd_recessive_fn, index_col=0)
        nd_recessive_pred = process_recessive_matrix(nd_recessive_pred)

        nd_recessive_pred_dict = dict(zip(nd_recessive_pred.index.tolist(), nd_recessive_pred.predict.tolist()))
        nd_recessive_conf_dict = dict(zip(nd_recessive_pred.index.tolist(), nd_recessive_pred['confidence'].tolist()))
        nd_recessive_rank_dict = dict(zip(nd_recessive_pred.index.tolist(), nd_recessive_pred.ranking.tolist()))
        nd_recessive_conf_lvl_dict = dict(zip(nd_recessive_pred.index.tolist(), nd_recessive_pred['confidence level'].tolist()))
        nd_recessive_other_var = dict(zip(nd_recessive_pred.index.tolist(), nd_recessive_pred.var2.tolist()))

        nd_rec_preds, nd_rec_confs, nd_rec_ranks, nd_rec_var2s = [], [], [], []
        nd_rec_lvls = []
        for var in expanded_df.index.tolist():
            if var in nd_recessive_pred_dict:
                nd_rec_preds.append(nd_recessive_pred_dict[var])
                nd_rec_confs.append(nd_recessive_conf_dict[var])
                nd_rec_ranks.append(nd_recessive_rank_dict[var])
                nd_rec_var2s.append(nd_recessive_other_var[var])
                nd_rec_lvls.append(nd_recessive_conf_lvl_dict[var])
            else:
                nd_rec_preds.append(-1)
                nd_rec_confs.append(-1)
                nd_rec_ranks.append(99999)
                nd_rec_var2s.append('NA')
                nd_rec_lvls.append('Unsolved')

        expanded_df['predict (nd recessive)'] = nd_rec_preds
        expanded_df['confidence (nd recessive)'] = nd_rec_confs
        expanded_df['confidence level (nd recessive)'] = nd_rec_lvls
        expanded_df['ranking (nd recessive)'] = nd_rec_ranks
        expanded_df['nd recessive var2'] = nd_rec_var2s
    else:
        expanded_df['predict (nd recessive)'] = -1
        expanded_df['confidence (nd recessive)'] = -1
        expanded_df['confidence level (nd recessive)'] = 99999
        expanded_df['ranking (nd recessive)'] = 'NA'
        expanded_df['nd recessive var2'] = 'Unsolved'

    return expanded_df

if __name__ == "__main__":
    from tqdm import tqdm
    import os
    prj_folder = 'projects/Texome/'
    data_folder = 'data/Texome/'
    sample_ids = os.listdir(data_folder)
    for sample_id in tqdm(sample_ids):
        integrated_df = integrate_output(prj_folder, data_folder, sample_id)
        if not os.path.exists(f"{prj_folder}/AIM_integrated"):
            os.mkdir(f"{prj_folder}/AIM_integrated")
        integrated_df.to_csv(f'{prj_folder}/AIM_integrated/{sample_id}_integrated.csv')

