#!/usr/bin/env/python
'''
Script created by Gemma I. Martinez-Redondo to obtain semantic similarity of GO terms between two lists.
Usage:
        python calculate_semsim.py -i inputpath [-o outfile.txt] [-c BP] [-g go.obo]
'''

from pygosemsim import graph
from pygosemsim import similarity
from pygosemsim import term_set
import functools, argparse, re
from pathlib import Path
import pandas as pd
import numpy as np

#Define parsing of the command
def parser():
	args = argparse.ArgumentParser(description='Obtain the semantic similarity matrix of GO terms between a group of entities (OGs, species...). Wang similarity will be used with BMA algortithm for getting a semantic similarity per pair entities (and not per pair of GO terms).')
	args.add_argument('-i', '--inpath', required=True, help="Name of the path containing one file per identity you want to compare functionally (OGs, species,...). One row per gene; first column contains the gene name, second column the list of GO terms separated by ', '. No header.")
	args.add_argument('-o', '--outfile', required=False, help="Path to the output file. If not provided, 'SimMat.out' will be used as default.")
	args.add_argument('-g', '--gofile', required=False, help="Path to the Gene Ontology go.obo file. If not provided, './go' will be used.")
	args.add_argument('-c', '--gocategory', required=False, default="BP", help="GO category for which to obtain the semantic similarity matrix. If not provided, biological process (BP) is used by default. Options: {BP, MF, CC}")
	args=args.parse_args()
	return args

#Obtain arguments
inpath = parser().inpath
outfile = parser().outfile
gofile = parser().gofile
gocategory = parser().gocategory

#Use defaults if optional arguments not given
if not outfile:
	outfile='SimMat.out'

if not gofile:
	gofile='/mnt/netapp2/Store_csbyegim/Metazoa_OGs/constellatoR/go_2022'
else:
	p=Path(gofile)
	extensions="".join(p.suffixes)
	gofile=str(p).replace(extensions, "")

if not gocategory:
	gocategory = "BP"
else:
	if gocategory not in ["BP","MF","CC"]:
		raise Exception("Invalid GO category provided! Please, select one among the following: {BP, MF, CC}")

#Read GO file
G=graph.from_resource(gofile)

#Define functions
#Method for semantic similarity
sf = functools.partial(term_set.sim_func, G, similarity.wang)
#Obtain go terms per input file (output: DataFrame)
def obtain_go_terms_per_file(file):
	df = pd.read_csv(file, sep="\t", header=None)
	df.columns = ['Gene','GO_terms']
	# Convert the GO_terms column into a single list
	all_go_terms = df['GO_terms'].str.split(', ').explode().dropna().tolist()
	# Convert it into a DataFrame
	df_go = pd.DataFrame({'GO_terms': all_go_terms})
	return df_go

def parse_obo(FILEPATH):
    regex = re.compile("alt_id: (GO:.*)\n")
    with open(FILEPATH, "r") as fread:
        data = [
            [y]
            + [
                y.split(": ")[1].strip()
                for y in x.strip().split("\n")
                if y and ":" in y
            ][1:3]
            for x in fread.read().replace(";", ",").split("\n\n")
            if x and ("[Term]" in x)
            for y in (
                [x.strip().split("\n")[1].split(": ")[
                    1].strip()] + regex.findall(x)
            )
        ]
    return {x[0]: x[1::] for x in data}

def obtain_category_for_goterms(golist,G):
	new_golist={}
	for goterm in golist:
		if obo.get(goterm, None) is not None:
			category = obo[goterm][1]
			if category not in new_golist.keys():
				new_golist[category]=[goterm]
			else:
				new_golist[category].append(goterm)
	return new_golist

obo=parse_obo(gofile+".obo")

# Load GO terms per OG (file)
file_list = [f for f in Path(inpath).iterdir()]
og_dict = {}
for file in file_list:
	OG=file.stem
	og_dict[OG]=obtain_go_terms_per_file(file)

# Get OG-GOs dataframe
og_final = pd.DataFrame({
    'OG': list(og_dict.keys()),
    'GO_terms': [df['GO_terms'].tolist() for df in og_dict.values()]
})

# Get information of GO categories
og_final[['BP', 'MF', 'CC']] = og_final['GO_terms'].apply(
    lambda go_list: pd.Series(obtain_category_for_goterms(go_list, obo)))
og_final = og_final.drop(columns=['GO_terms'])

# Now, let's get for the requested category, the OG-OG matrix to get the semantic similarity
# Extract OG names and GO term lists of the specified categgory
valid_rows = og_final[gocategory].notna()  # Mask for non-NaN rows. Otherwise, if NaN, it will raise an error
og_list = og_final.loc[valid_rows, 'OG'].values  # Keep only valid OGs (not NaN)
go_terms = np.array(og_final.loc[valid_rows, gocategory].values, dtype=object)  # Keep only valid GO terms and store in NumPy
# Initialize an empty NumPy matrix (fill with NaN to indicate uncomputed values)
sim_matrix = np.full((len(og_list), len(og_list)), np.nan)
# Compute only the upper triangle (excluding diagonal)
rows, cols = np.triu_indices(len(og_list))  # k=1 excludes diagonal
# Vectorized computation for upper triangle
vec_sim = np.frompyfunc(lambda x, y: term_set.sim_bma(x, y, sf), 2, 1)
sim_matrix[rows, cols] = vec_sim(go_terms[rows], go_terms[cols])
# Fill lower triangle
sim_matrix[cols, rows] = sim_matrix[rows, cols]
# Convert to a Pandas DataFrame
sim_df = pd.DataFrame(sim_matrix, index=og_list, columns=og_list)

# Save to file
sim_df.to_csv(outfile, sep="\t", header=True, index=True, na_rep='NA')
