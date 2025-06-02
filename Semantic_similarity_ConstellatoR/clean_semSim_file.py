#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: klara_el
This script takes semantic similarity matrix as input and
removes the columns and rows with NA values
"""

import os
import pandas as pd
import numpy as np

BasePath = "./constellatoR/Per_Experiment_PerOG"

experiment = "CD" #Change accordingly
experiment_file = "semSim_perOG_" + experiment + ".txt"
experiment_file_path = os.path.join(BasePath, experiment, experiment_file)


data = []
row_labels = []

with open(experiment_file_path, 'r') as file:
    lines = file.readlines()

# First line = column headers
col_labels = lines[0].strip().split()

# Parse the rest
for line in lines[1:]:
    parts = line.strip().split()
    if not parts:
        continue
    row_labels.append(parts[0])  # first column = row label
    row_values = [float(x) if x != "NA" else np.nan for x in parts[1:]]
    data.append(row_values)

# Create DataFrame
df = pd.DataFrame(data, index=row_labels, columns=col_labels)

# Get rows where all values are NaN
rows_all_na = df[df.isna().all(axis=1)]

# Get columns where all values are NaN
cols_all_na = df.loc[:, df.isna().all(axis=0)]

# Extract the row labels for rows that are all NaN
rows_with_all_na = rows_all_na.index.tolist()

# Extract the column labels for columns that are all NaN
cols_with_all_na = cols_all_na.columns.tolist()

df_cleaned = df.dropna(axis=0, how='all')  # Remove rows with all NaN
df_cleaned = df_cleaned.dropna(axis=1, how='all')  # Remove columns with all NaN
print("Shape before cleaning:", df.shape)
print("Shape after cleaning:", df_cleaned.shape)
# Create the cleaned file path with "cleaned" appended to the filename
cleaned_file_name = experiment_file.replace(".txt", "_cleaned.txt")
cleaned_file_path = os.path.join(BasePath, experiment, cleaned_file_name)

# Save the cleaned DataFrame to the new file
df_cleaned.to_csv(cleaned_file_path, sep='\t', header=True, index=True)

print(f"Cleaned file saved at: {cleaned_file_path}")
