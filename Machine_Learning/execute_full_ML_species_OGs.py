#!/mnt/netapp2/Store_csbyegim/Metazoa_OGs/ML_Species-OGs/full_dataset/ml_venv/bin/python

import dask.dataframe as dd
import dask.array as da
from os.path import join, exists
from sklearn.model_selection import StratifiedKFold

from dask.distributed import LocalCluster

import xgboost,dask_ml,sklearn,time,random,pickle
from xgboost import dask as dxgb
from xgboost.dask import DaskDMatrix
from dask_ml.metrics import mean_squared_error,accuracy_score
from sklearn.metrics import balanced_accuracy_score,confusion_matrix

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import dask
import seaborn as sns

from statistics import mean
from itertools import chain

from dRFEtools_xgboost_plus_orig_feat_sel import *

#Function to make it work in CESGA
def giant_function():
	global pd
	global plt
	global np
	global dd
	global da
	global sns
	#Start dask client
	dask.config.set({"distributed.scheduler.worker-saturation":  1.1})

	cluster = LocalCluster(n_workers=10,
	                       memory_limit='100GB')
	client = cluster.get_client()

	RANDOM_STATE = 42 #To ensure reproducibility
	random.seed(RANDOM_STATE)
	ogs_sp_file="/mnt/netapp2/Store_csbyegim/Metazoa_OGs/ML_Species-OGs/full_dataset/filtered_transposed_Species_OGs_habitat_phylum.tsv"

	print("Reading file and preparing input")
	#Read file
	df = dd.read_csv(ogs_sp_file+"/*", sep='\t', sample=5000000, blocksize=None)
	#REMEMBER: Group_ID (OGs) start in 1

	#Select columns corresponding to OGs and ID for index (all but 'Code', 'Habitat', 'Phylum' and 'File')
	X = df.drop(columns=['Code','Habitat','Phylum','File'])
	#Select column to predict ('Habitat'), the extra column for stratification ('Phylum'), and ID for index.
	y = df[['id','Habitat','Phylum']]

	#Set index for X and y
	X=X.set_index('id')
	y=y.set_index('id')

	#Get stratification indexes from ML Species-GOs. When using traditional StratifiedKFold function (which comes from sklearn) there is some problem with shuffle and the dask dataframes. See Debugging page
	with open("idx_skf_sp.pickle","rb") as f:
		idx_skf_sp=pickle.load(f)

	#Modify idx_skf to get new ids
	codes_sp=list(list(df[['Code']].compute().values))
	id_sp=list(list(df[['id']].compute().values))
	mapping_dict=pd.Series([id_sp[i][0] for i in range(964)], index=[codes_sp[i][0] for i in range(964)]).to_dict()

	# Replace the numbers in the copy (idx_skf_sp) using the mapping
	idx_skf = []
	for fold in idx_skf_sp:
	    fold_sp = []  # This will hold the replaced arrays for each fold
	    for arr in fold:
	        # Replace numbers with strings using the mapping
	        fold_sp.append(np.array([mapping_dict[num] for num in arr]))
	    idx_skf.append(fold_sp)

	#Stratification plots
	df_stratkfold=y.compute()
	fold=1
	for tr,tt in idx_skf:
	    indices = np.array([np.nan] * (len(tr)+len(tt)))
	    indices[tt] = 1 #Testing
	    indices[tr] = 0 #Training
	    df_stratkfold["F"+str(fold)]=indices
	    fold=fold+1

	df_stratkfold=df_stratkfold.sort_values('Habitat')
	df_stratkfold['Phylum'] = pd.Categorical(df_stratkfold['Phylum'], ["NO_BILATERIA","XENACOELOMORPHA","DEUTEROSTOMIA","CRANIATA","ECDYSOZOA","NEMATODA","TARDIGRADA","ONYCHOPHORA","ARTHROPODA","LOPHOTROCHOZOA","PLATYHELMINTHES","MOLLUSCA","ANNELIDA","NEMERTEA"])
	df_stratkfold=df_stratkfold.sort_values("Phylum")
	df_stratkfold=df_stratkfold.replace(["NO_BILATERIA","XENACOELOMORPHA","DEUTEROSTOMIA","CRANIATA","ECDYSOZOA","NEMATODA","TARDIGRADA","ONYCHOPHORA","ARTHROPODA","LOPHOTROCHOZOA","PLATYHELMINTHES","MOLLUSCA","ANNELIDA","NEMERTEA"], [0,1,2,3,4,5,6,7,8,9,10,11,12,13])

	# Assume df_stratkfold is your DataFrame, cmap_cv for train/test, and cmap_data for Habitat/Phylum.
	cmap_cv = plt.cm.coolwarm
	cmap_data = plt.cm.tab20
	# Manual color mapping for Habitat
	habitat_colors = {0: 'skyblue', 1: 'saddlebrown'}

	# Phylum names corresponding to their numerical codes
	phylum_names = ["NO_BILATERIA", "XENACOELOMORPHA", "DEUTEROSTOMIA", "CRANIATA", "ECDYSOZOA", 
	                "NEMATODA", "TARDIGRADA", "ONYCHOPHORA", "ARTHROPODA", "LOPHOTROCHOZOA", 
	                "PLATYHELMINTHES", "MOLLUSCA", "ANNELIDA", "NEMERTEA"]

	plt.figure()
	fig, ax = plt.subplots(figsize=(12, 8))

	# Scatter plot for train/test sets (df_stratkfold.columns[2:])
	for col in df_stratkfold.columns[2:]:
	    ax.scatter(
	        range(len(df_stratkfold[col])),
	        np.repeat(int(col[-1]) + 1, len(df_stratkfold[col])),
	        c=df_stratkfold[col],
	        marker="_",
	        lw=10,
	        cmap=cmap_cv,
	        vmin=-0.2,
	        vmax=1.2,
	    )

	# Scatter plot for Habitat
	for i, habitat_value in enumerate(df_stratkfold['Habitat']):
	    ax.scatter(
	        i, 1,  # Y position is 1 for Habitat
	        color=habitat_colors[habitat_value],
	        marker="_",
	        lw=10
	    )

	# Scatter plot for Phylum
	ax.scatter(
	    range(len(df_stratkfold['Phylum'])),
	    np.repeat(0, len(df_stratkfold['Phylum'])),
	    c=df_stratkfold['Phylum'],
	    marker="_",
	    lw=10,
	    cmap=cmap_data
	)

	# Formatting y-ticks and labels
	yticklabels = ["Phylum", "Habitat"] + list(df_stratkfold.columns[2:])
	ax.set(
	    yticks=np.arange(len(df_stratkfold.columns)),
	    yticklabels=yticklabels
	)
	ax.set_title("StratifiedKFold", fontsize=15)

	# 1. First legend for Training and Testing sets
	train_test_legend = ax.legend(
	    [Patch(color=cmap_cv(0.8)), Patch(color=cmap_cv(0.02))],
	    ["Testing set", "Training set"],
	    loc="center left",
	    bbox_to_anchor=(1.05, 0.8)
	)

	# 2. Second legend for Habitat (Assuming 0 is terrestrial and 1 is aquatic)
	habitat_legend = plt.legend(
	    [Patch(color='skyblue'), Patch(color='saddlebrown')],
	    ["Aquatic", "Terrestrial"],
	    loc="center left",
	    bbox_to_anchor=(1.05, 0.7)
	)

	# 3. Third legend for Phylum (Phylum names instead of numbers)
	phylum_patches = [Patch(color=cmap_data(i / 13.0)) for i in range(14)]
	phylum_legend = plt.legend(
	    phylum_patches,
	    phylum_names,  # Use the names of the Phylum instead of numbers
	    loc="center left",
	    bbox_to_anchor=(1.05, 0.4)
	)

	# Add the first and second legends manually back to the plot
	ax.add_artist(train_test_legend)
	ax.add_artist(habitat_legend)

	# Save the plot
	plt.savefig("StratifiedKFold_data.png", bbox_inches="tight")
	plt.close()

	#Habitat plot
	# Number of folds
	n_folds = len(idx_skf)

	# Create subplots for visualizing Habitat distributions across folds
	fig, axes = plt.subplots(n_folds, 1, sharex='col', sharey='row', figsize=(12, 5 * n_folds))
	plt.subplots_adjust(hspace=0.7)  # Increase spacing between subplots

	for fold in range(n_folds):
	    idx_train = idx_skf[fold][0]
	    idx_test = idx_skf[fold][1]    
	    # Get the training and testing data for this fold and convert to pandas DataFrame
	    y_train = y.loc[list(idx_train)].compute()  # Convert Dask DataFrame to pandas DataFrame
	    y_test = y.loc[list(idx_test)].compute()    # Convert Dask DataFrame to pandas DataFrame
	    # Map Habitat values to 'Aquatic' and 'Terrestrial'
	    y_train['Habitat'] = y_train['Habitat'].map({0: 'Marine', 1: 'NonMarine'})
	    y_test['Habitat'] = y_test['Habitat'].map({0: 'Marine', 1: 'NonMarine'})
	    # Compute proportions for Habitat in Train and Test sets
	    habitat_train_prop = y_train['Habitat'].value_counts(normalize=True).reset_index()
	    habitat_train_prop.columns = ['Habitat', 'Proportion']
	    habitat_train_prop['Set'] = 'Train'   
	    habitat_test_prop = y_test['Habitat'].value_counts(normalize=True).reset_index()
	    habitat_test_prop.columns = ['Habitat', 'Proportion']
	    habitat_test_prop['Set'] = 'Test' 
	    # Concatenate Train and Test proportions for Habitat
	    habitat_combined = pd.concat([habitat_train_prop, habitat_test_prop])
	    # Plot Habitat proportions for Train and Test in the same subplot with custom colors
	    sns.barplot(x='Habitat', y='Proportion', hue='Set', data=habitat_combined, palette='plasma', ax=axes[fold])
	    axes[fold].set_title(f'Fold {fold + 1} Habitat Proportions (Train vs Test)')
	    axes[fold].set_xlabel("Habitat")  # X-axis label
	    axes[fold].set_ylabel("Proportion")  # Y-axis label

	plt.tight_layout()

	# Save the plot
	plt.savefig("Strat_Habitat_distribution.png", bbox_inches="tight")
	plt.close()

	#Phylum plot
	# Number of folds
	n_folds = len(idx_skf)

	# Create subplots for visualizing Phylum distributions across folds
	fig, axes = plt.subplots(n_folds, 1, sharex='col', sharey='row', figsize=(10, 5 * n_folds))
	plt.subplots_adjust(hspace=0.5)

	for fold in range(n_folds):
	    idx_train = idx_skf[fold][0]
	    idx_test = idx_skf[fold][1]    
	    
	    # Get the training and testing data for this fold and convert to pandas DataFrame
	    y_train = y.loc[list(idx_train)].compute()  # Convert Dask DataFrame to pandas DataFrame
	    y_test = y.loc[list(idx_test)].compute()    # Convert Dask DataFrame to pandas DataFrame
	    
	    # Compute proportions for Phylum in Train and Test sets
	    phylum_train_prop = y_train['Phylum'].value_counts(normalize=True).reset_index()
	    phylum_train_prop.columns = ['Phylum', 'Proportion']
	    phylum_train_prop['Set'] = 'Train'
	    
	    phylum_test_prop = y_test['Phylum'].value_counts(normalize=True).reset_index()
	    phylum_test_prop.columns = ['Phylum', 'Proportion']
	    phylum_test_prop['Set'] = 'Test'
	    
	    # Concatenate Train and Test proportions for Phylum
	    phylum_combined = pd.concat([phylum_train_prop, phylum_test_prop])
	    
	    # Plot Phylum proportions for Train and Test in the same subplot
	    sns.barplot(x='Phylum', y='Proportion', hue='Set', data=phylum_combined, palette='plasma', ax=axes[fold])
	    axes[fold].set_title(f'Fold {fold + 1} Phylum Proportions (Train vs Test)')
	    
	    # Rotate x-axis labels by 45 degrees
	    axes[fold].tick_params(axis='x', rotation=45)

	plt.tight_layout()

	# Save the plot
	plt.savefig("Strat_Phylum_distribution.png", bbox_inches="tight")
	plt.close()

	#Name of features. In this dataset, we are working with numbers
	#features = ["%d" % x for x in range(1,X.shape[1]+1)] #["f%d" % x for x in range(X.shape[1])] #with f it adds an f to feature name
	features = np.array(X.columns) #Faster than previous option (above)
	out_dir="."

	starting_time = time.time()
	starting_hour = time.asctime(time.localtime())

	print("First model execution started at: "+starting_hour)

	#Let's execute manually the xgboost_rfe and _xgboost_fe functions with just one iteration to get intermediate data
	d = {}
	pfirst = None

	#From _xgboost_fe, let's train the model for all folds
	indices_feats = np.arange(X.shape[1])
	nf = len(features)
	p={}
	fold=1
	#For each fold
	for idx_train, idx_test in idx_skf:
		print("\nFold {}".format(fold))
		#Get data subsets
		X_train, X_test, y_train, y_test = get_strat_train_test(X,y['Habitat'],idx_train,idx_test)
		p[fold] = xgboost_fe_step(client, X_train, y_train, X_test, y_test, nf, features, out_dir)
		#Save results (p) at each fold (progressive, just in case something happens).
		with open("p_folds_until_"+str(fold)+"_stats_feat_importance_full_model.pickle","wb") as f:
			pickle.dump(p,f,pickle.HIGHEST_PROTOCOL)
		fold = fold+1

	elapsed_time = time.time() - starting_time
	ending_hour = time.asctime(time.localtime())
	print("End of model training.\nElapsed time:", elapsed_time)

	starting_time = time.time()
	starting_hour = time.asctime(time.localtime())

	#Now, we need to combine the results and see how many features are 0
	features = np.array(features)
	#Get feature importance (gain)
	feature_importances_stats_np=np.array([p[f]['feat_importance'] for f in p])
	#Get train and test accuracy
	Acc_t=np.array([p[f]['accuracy_score_train'] for f in p])
	Acc_v=np.array([p[f]['accuracy_score'] for f in p])
	#Compute the per-fold common term based on accuracies
	epsilon=0.00001
	G = (Acc_t + Acc_v)/(abs((Acc_t - Acc_v))+epsilon)
	#Calculate ranking of features
	features_importance_new_rank = feature_importances_stats_np.T@G
	#Save results (features_importance_new_rank)
	with open("ranked_feat_importance_full_model.pickle","wb") as f:
		pickle.dump(features_importance_new_rank,f,pickle.HIGHEST_PROTOCOL)

	elapsed_time = time.time() - starting_time
	ending_hour = time.asctime(time.localtime())
	print("End of features ranking.\nElapsed time:", elapsed_time)

	#Check features that are 0
	sum(features_importance_new_rank==0)
	#Features that are not 0
	len(features_importance_new_rank)-sum(features_importance_new_rank==0)

	#First, we'll remove all features that are 0

	#Last part of selected = features_rank_fnc(features, p, nf, out_dir)
	rank = np.argsort(features_importance_new_rank)[::-1] # reverse sort. Take into account that if number of features to keep is higher than the number of features that were used in all folds, it will still select features that have a value of 0.
	#Get number of features that are not 0
	n_features_to_keep=len(features_importance_new_rank)-sum(features_importance_new_rank==0)
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
	print("Features that are not 0 have been selected")



if __name__ == '__main__':
        giant_function()
