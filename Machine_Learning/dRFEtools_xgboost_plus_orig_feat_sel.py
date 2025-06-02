#!/usr/bin/env python

"""
This script modifies the functions from the dRFEtools python package available in github (https://github.com/LieberInstitute/dRFEtools) for its application using XGBoost and Dask. In addition, it incorporates the ranking system described in https://link.springer.com/chapter/10.1007/11893295_1.

Modifications made by Gemma I. Martínez-Redondo

"""

import numpy as np
import pandas as pd
import xgboost
import time
from os.path import join, exists
from itertools import chain
from sklearn.metrics import (
		balanced_accuracy_score
	)

def features_rank_fnc(features, feature_importances_stats, n_features_to_keep, out_dir, epsilon=0.00001):
	"""
	Ranks features and writes the results to a file
	Args:
		features: A vector of feature names
		feature_importances_stats: list containing feature importances for each of the folds as well as accuracy and other stats
		n_features_to_keep (int): Number of features to keep.
		out_dir (str): Output directory for text file. Default is current
					   directory.
		epsilon (int): Small number used to avoid division by 0 when calculating new features ranking. Default is 10e-5

	Returns:
		List of selected features

	Writes:
	   Text file: Ranked features by fold tab-delimitated text file
	"""
	if not isinstance(n_features_to_keep, int) or n_features_to_keep < 0:
		raise ValueError("n_features_to_keep must be a non-negative integer")
	features = np.array(features)
	#Get feature importance (gain)
	feature_importances_stats_np=np.array([feature_importances_stats[f]['feat_importance'] for f in feature_importances_stats])
	#Get train and test accuracy
	Acc_t=np.array([feature_importances_stats[f]['accuracy_score_train'] for f in feature_importances_stats])
	Acc_v=np.array([feature_importances_stats[f]['accuracy_score'] for f in feature_importances_stats])
	#Compute the per-fold common term based on accuracies
	G = (Acc_t + Acc_v)/(abs((Acc_t - Acc_v))+epsilon)
	#Calculate ranking of features
	features_importance_new_rank = feature_importances_stats_np.T@G
	rank = np.argsort(features_importance_new_rank)[::-1] # reverse sort. Take into account that if number of features to keep is higher than the number of features that were used in all folds, it will still select features that have a value of 0.
	eliminated = rank[n_features_to_keep:]
	selected = rank[:n_features_to_keep]
	if len(eliminated) == 0:
		rank_df = pd.DataFrame({
			'Geneid': features[rank],
			'Rank': 1
		})
	else:
		rank_df = pd.DataFrame({
			'Geneid': features[eliminated],
			'Rank': np.arange(n_features_to_keep + 1,
							  n_features_to_keep + 1 + len(eliminated))
		})
	output_file = join(out_dir, "rank_features.txt")
	rank_df.sort_values('Rank', ascending=False).to_csv(
		output_file,
		sep='\t',
		mode='a',
		index=False,
		header=not exists(output_file)
	)
	return selected


def get_predictions(client, model, X):
	"""
	Calculates the prediction of probabilities and labels using trained model.

	Args:
		client: Dask client
		model: XGBoost model trained on training data set
		X: dask array of features from data set

	Returns:
		tuple: predicted probabilities and labels from data set
	"""
	probs=xgboost.dask.predict(client, model, X)
	return probs, probs.round()


def score_accuracy(predictions, y):
	"""
	Calculates the balanced accuracy score.

	Args:
		client: Dask client
		predictions: dask array or numpy array predicted sample labels from data set
		y: dask array or numpy arrayof sample labels from data set

	Returns:
		float: balanced accuracy score
	"""
	if not isinstance(predictions, np.ndarray):
		pred = np.array([value for value in predictions.compute()]) #Converting it to np array is faster than using predictions directly)
	else:
		pred = predictions
	if not isinstance(y, np.ndarray):
		y_np = np.array([value for value in y.compute()]) #Converting it to np array is faster than using predictions directly)
	else:
		y_np = y
	return balanced_accuracy_score(y_np, pred) #sklearn.metrics.balanced_accuracy_score(y_true, y_pred)


def score_logloss(predictions_probs, y_test):
	"""
	Calculates the Logloss of test dataset.

	Args:
		prediction_probss: dask array or numpy array of predicted probabilities from test data set
		y_test: dask array or np array of sample labels from training data set

	Returns:
		float: logloss score
	"""
	if not isinstance(y_test, np.ndarray):
		y_test_np=np.array([value for value in y_test.compute()])
	else:
		y_test_np=y_test
	if not isinstance(predictions_probs, np.ndarray):
		predictions_probs_np=np.array([value for value in predictions_probs.compute()])
	else:
		predictions_probs_np=predictions_probs
	test_logloss=-(y_test_np*np.log(predictions_probs_np)+(1-y_test_np)*np.log(1-predictions_probs_np))
	#above formula is calculated for each value. We need the average
	return test_logloss.mean()
	

def get_strat_train_test(X, y, idx_train, idx_test):   
    # Directly index the DataFrame to create train and test subsets
    X_train = X.loc[list(idx_train)] #Error if not changed to list
    X_test = X.loc[list(idx_test)]
    y_train = y.loc[list(idx_train)]
    y_test = y.loc[list(idx_test)]
    return X_train, X_test, y_train, y_test


def xgboost_fe_step(client, X_train, y_train, X_test, y_test, n_features_to_keep, features, out_dir):
	"""
	Eliminates features step-by-step.

	Args:
		client: Dask client
		X_train: a data frame of training data
		y_train: np.ndarray of sample labels from training data set
		n_features_to_keep: number of features to keep
		features: np.ndarray of feature names
		out_dir: output directory. default '.'

	Returns:
		dict: a dictionary containing feature elimination results

	Raises:
		ValueError: If n_features_to_keep is greater than the number of features in X
	"""
	if n_features_to_keep > X_train.shape[1]:
		raise ValueError("n_features_to_keep cannot be greater than the number of features in X_train")
	d_train = xgboost.dask.DaskDMatrix(client, X_train, y_train, enable_categorical=True)
	model = xgboost.dask.train(
	    client,
	    {"tree_method": "hist",
	    'objective': 'binary:logistic',
	    'eval_metric':"logloss",},
	    d_train,
	    num_boost_round=100,
	    evals=[(d_train, "train")],
	)
	feature_importances = model['booster'].get_score(importance_type='gain')
	feature_importances = {f: feature_importances.get(f, 0) for f in features}  #Fill with 0s, not used features
	feature_importances_np=np.array(list(feature_importances.values()))
	#Get training accuracy
	train_pred_probs, train_pred = get_predictions(client, model, X_train)
	train_pred_np = np.array([value for value in train_pred.compute()])
	#Get validation accuracy
	predictions_probs, predictions = get_predictions(client, model, X_test)
	y_test_np=np.array([value for value in y_test.compute()])
	pred_np = np.array([value for value in predictions.compute()])
	result = {
		#'n_features': X_train.shape[1],
		#'selected': selected
		'logloss_train': model['history']['train']['logloss'][-1],
		'logloss_test': score_logloss(predictions_probs, y_test_np),
		'feat_importance': feature_importances_np,
		'accuracy_score': score_accuracy(pred_np, y_test_np),
		'accuracy_score_train': score_accuracy(train_pred_np, y_train)
	}
	return result


def xgboost_fe(client, X, y, idx_skf, n_features_iter, features, out_dir, starting_time):
	"""
	Iterates over features to be eliminated step-by-step.

	Args:
   		client: Dask client
		X: DataFrame or np.ndarray of full data
		y: np.ndarray of sample labels from all data set
		idx_skf: Indexes to split dataset in all folds 
		n_features_iter: Iterator for number of features to keep loop
		features: np.ndarray of feature names
		out_dir (str): Output directory.
		starting_time: Starting time of whole pipeline. Used to calculate execution time at each iteration.

	Returns:
		tuple: Feature elimination results for each iteration

	Raises:
		ValueError: If X and features have different number of columns
	"""
	if X.shape[1] != len(features):
		raise ValueError("Number of columns in X must match the length of features")
	indices_feats = np.arange(X.shape[1])
	#For each iteration of RFE
	for nf in chain(n_features_iter, [1]):
		print(nf)
		p={}
		fold=1
		#For each fold
		for idx_train, idx_test in idx_skf:
			print("\nFold {}".format(fold))
			#Get data subsets
			X_train, X_test, y_train, y_test = get_strat_train_test(X,y,idx_train,idx_test)
			p[fold] = xgboost_fe_step(client, X_train, y_train, X_test, y_test, nf, features, out_dir)
			fold = fold+1
		#Obtain features selected given ranking of features accross folds
		selected = features_rank_fnc(features, p, nf, out_dir)
		elapsed_time = time.time() - starting_time
		yield X_train.shape[1], [p[f]['logloss_train'] for f in p], [p[f]['logloss_test'] for f in p], [p[f]['accuracy_score'] for f in p], [p[f]['accuracy_score_train'] for f in p], elapsed_time, indices_feats
		indices_feats = indices_feats[selected]
		features = np.array(features)[selected].tolist()
		X = X.loc[:, features]


def n_features_iter(nf: int, keep_rate: float) -> int:
	"""
	Determines the features to keep.

	Args:
		nf (int): Current number of features
		keep_rate (float): Percentage of features to keep

	Returns:
		int: Number of features to keep
	"""
	while nf != 1:
		nf = max(1, int(nf * keep_rate))
		yield nf


def xgboost_rfe(client, X, y, idx_skf, features, out_dir='.', elimination_rate=0.2):
	"""
	Runs XGBoost feature elimination step over iterator process.

	Args:
		client: Dask client
		X (DataFrame): Full data
		y (array-like): Sample labels from data set
		idx_skf: Indexes to split dataset in all folds 
		features (array-like): Feature names
		out_dir (str): Output directory. Default '.'
		elimination_rate (float): Percent rate to reduce feature list. Default 0.2

	Returns:
		tuple: Dictionary with elimination results, and first elimination step results
	"""
	if not 0 < elimination_rate < 1:
		raise ValueError("elimination_rate must be between 0 and 1")
	d = {}
	pfirst = None
	keep_rate = 1 - elimination_rate
	starting_time = time.time()
	for p in xgboost_fe(client, X, y, idx_skf, n_features_iter(X.shape[1], keep_rate), features, out_dir, starting_time):
		if pfirst is None:
			pfirst = p
		d[p[0]] = p
	return d, pfirst
