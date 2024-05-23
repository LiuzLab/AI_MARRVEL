import numpy as np
import pandas as pd

__author__ = "Chaozhong Liu"
__email__ = "chaozhol@bcm.edu"

###adds known strong/weak labels to matrix (not used in production pipeline)

def add_label(feature_mtx, label_file, sample_id):
	"""
	feature_mtx: file path of matrix after feature engineering

	label_file: path of file containing all diagnosed causal variants with columns "BG_ID", "var" and "Confidence"

	sample_id: id of sample to be labeled
	"""

	patient = feature_mtx
	# is_causal - including strong and weak
	label_df = pd.read_csv(label_file)
	label_df = label_df[label_df['BG_ID'] == sample_id]
	label = label_df['var'].values.tolist()
	label = list(set(label))

	if len(label) == 0:
	    print('No causal variant found in label file: %s'%(sample_id))
	    return

	if sum(patient.index.isin(label)) == 0:
	    print('True label not found in matrix: %s.'%(sample_id))
	    return


	patient['is_causal'] = 0
	patient.loc[patient.index.isin(label),'is_causal'] = 1


	# is_strong - including strong only
	label_df = pd.read_csv(label_file)
	label_df = label_df.loc[label_df['Confidence']=='Strong',:]
	label_df = label_df[label_df['BG_ID'] == sample_id]
	label = label_df['var'].values.tolist()
	label = list(set(label))

	patient['is_strong'] = 0

	if len(label) == 0:
	    print('No strong causal variant found in label file: %s'%(sample_id))
	    return patient


	if sum(patient.index.isin(label)) == 0:
	    print('Strong label not found in matrix: %s.'%(sample_id))
	    return patient


	patient.loc[patient.index.isin(label),'is_strong'] = 1


	return patient


