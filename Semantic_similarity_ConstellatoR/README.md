## 1. Description

This directory contains information on the functional convergence analyses using semantic similarity of GO terms and constellatoR.

## 2. Workflow

### Building the Semantic Similarity matrix for ConstellatoR

To test functional convergence between different datasets, we are using **constellatoR** (https://github.com/MetazoaPhylogenomicsLab/constellatoR).

There are some incongruences between FANTASIA GO version and the version that constellatoR uses (which depends on the time a dependency is installed and can’t be easily changed). So, we first created a script to estimate the OG-OG semantic similarity matrix that constellatoR will use for creating the functional clusters. This uses pygosemsim (https://github.com/mojaie/pygosemsim/tree/master).

[calculate_semsim_fixed.py](calculate_semsim_fixed.py)

</aside>

This script uses as input a directory with all the files (one per OG) containing the GO terms per gene. Then, it estimates the semantic similarity between the different OGs and builds a matrix. It’s optimized as it uses vectorized computation with numpy.

### Functional convergence between gained/lost OGs in independent terrestrialisation events

Let’s try to use only the enriched GO terms and compare directly between clades.

- calculate_semsim_enriched_gos.py
    
    ```python
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
    	args = argparse.ArgumentParser(description='Obtain the semantic similarity matrix of enriched GO terms between a group of entities (OGs, species...). Wang similarity will be used with BMA algortithm for getting a semantic similarity per pair entities (and not per pair of GO terms).')
    	args.add_argument('-i', '--inpath', required=True, help="Name of the path containing one file per identity you want to compare functionally (OGs, species,...). One row per GO term. No header.")
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
    	df = pd.read_csv(file, sep=" ", header=None)
    	df.columns = ['GO_terms','p_value']
    	# Convert the GO_terms column into a single list
    	all_go_terms = df['GO_terms'].explode().dropna().tolist()
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
    og_final[gocategory] = og_final['GO_terms']
    og_final = og_final.drop(columns=['GO_terms'])
    
    # Now, let's get for the requested category, the OG-OG matrix to get the semantic similarity
    # Extract OG names and GO term lists of the specified categgory
    og_list = og_final['OG'].values
    go_terms = np.array(og_final[gocategory].values, dtype=object)  # Store GO terms as a NumPy array
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
    ```
    

In this case, the input is the output file of the topGO enrichment (plus ordmeta). Files have GOs in first column and p-values in the second. Just for one GO category, meaning that we will need to put as argument if it’s BP, MF or CC.

1. Gains
    
- Semantic similarity between clades using enriched GO terms
    
    We have previously done a GO enrichment per terrestrial clade. So, we have a subset of terms that we can use to test for convergence, to see if the GO terms that are significantly enriched on those clades are also similar between clades.
    
    Let’s rename the files and execute a modified version of the calculate_semsim.py script above.
    
    ```bash
    # Copy files (in local computer)
    rsync -av --progress /home/metazomics/Metazoa_analyses/Phylogenetic_trees_and_analyses/GO_enrichment/Results/gene_gain/per_node/BP_gained_*_topgo_enrichment.txt csbyegim@ft3.cesga.es:/mnt/netapp2/Store_csbyegim/Metazoa_OGs/constellatoR/gains/enriched_GO_terms_per_node
    
    # Remove species-specific files (in /mnt/netapp2/Store_csbyegim/Metazoa_OGs/constellatoR/gains/enriched_GO_terms_per_node)
    rm BP_gained_ANGR1_topgo_enrichment.txt BP_gained_CYCO1_topgo_enrichment.txt BP_gained_DRAW1_PELO1_Crassi_topgo_enrichment.txt BP_gained_ECRY2_topgo_enrichment.txt BP_gained_ETES1_topgo_enrichment.txt BP_gained_HRPE1_topgo_enrichment.txt BP_gained_OVER1_topgo_enrichment.txt BP_gained_PELE1_topgo_enrichment.txt BP_gained_PELO1_Crass_topgo_enrichment.txt BP_gained_PHEI1_topgo_enrichment.txt BP_gained_PPUN2_topgo_enrichment.txt BP_gained_RVAR1_topgo_enrichment.txt BP_gained_TRLO1_topgo_enrichment.txt
    
    # Rename files (in /mnt/netapp2/Store_csbyegim/Metazoa_OGs/constellatoR/gains/enriched_GO_terms_per_node)
    rename "BP_gained_" "" *
    rename "_topgo_enrichment" "" *
    
    # Get the semantic similarity matrix
    module load cesga/2020 python/3.9.9
    
    source $STORE2/Dark_proteome/Isoforms_Semsim/pygosemsim_venv/bin/activate
    
    WD=/mnt/netapp2/Store_csbyegim/Metazoa_OGs/constellatoR/gains
    INPUT_FOLDER=$WD/enriched_GO_terms_per_node/
    
    cd $WD
    
    python $WD/../calculate_semsim_enriched_gos.py -i $INPUT_FOLDER -o $WD/gains_GO_enrich_nodes_SimMat.out
    ```
    
    This script finished in 1h 40 mins. Let’s plot these results:
    
    ```python
    import seaborn as sns
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
    import matplotlib.colorbar
    import pandas as pd
    
    df_go = pd.read_csv("gains_GO_enrich_nodes_SimMat.out", sep="\t", index_col=0)
    #Change order of columns/rows (to follow phylogeny)
    desired_order = ["Tetrapoda","Rhabditida","Tardigrada","Onychophora","Arachnida","Myriapoda","Isopoda","Brachyura","Anomura","Hexapoda","Geoplanidae","Ellobiida","Stylomatophora","Crassiclitellata","Acteonemertidae"]  # Replace with actual row/column names
    df_go = df_go.reindex(index=desired_order, columns=desired_order)
    
    #For 2 triangular heatmaps
    #sns.heatmap(df_go,annot=True,mask=matrix_go,cmap="BuPu")
    #sns.heatmap(df_clusters,annot=True,mask=matrix_clusters,cmap="YlGnBu")
    
    plt.figure(figsize=(10, 10))
    sns.heatmap(df_go,annot=True,cmap="YlGnBu")
    #plt.xticks(rotation=45)
    #plt.yticks(rotation=45)
    plt.savefig("GO_enrich_gains_convergence_nodes.svg",bbox_inches="tight")
    ```
    
    ![GO_enrich_gains_convergence_nodes.png](images/GO_enrich_gains_convergence_nodes.png)
    

2. Losses

As losses are much higher than gains, there are too many OGs. We will proceed with the same approach that we used in the end for gains: use the GO terms enriched in the losses.

```bash
# Copy files (in local computer)
rsync -av --progress /home/metazomics/Metazoa_analyses/Phylogenetic_trees_and_analyses/GO_enrichment/Results/gene_loss/per_node/BP_lost_*_topgo_enrichment.txt csbyegim@ft3.cesga.es:/mnt/netapp2/Store_csbyegim/Metazoa_OGs/constellatoR/losses/enriched_GO_terms_per_node

# Remove species-specific files (in /mnt/netapp2/Store_csbyegim/Metazoa_OGs/constellatoR/losses/enriched_GO_terms_per_node)
rm BP_lost_ANGR1_topgo_enrichment.txt BP_lost_CYCO1_topgo_enrichment.txt BP_lost_DRAW1_PELO1_Crassi_topgo_enrichment.txt BP_lost_ECRY2_topgo_enrichment.txt BP_lost_ETES1_topgo_enrichment.txt BP_lost_HRPE1_topgo_enrichment.txt BP_lost_OVER1_topgo_enrichment.txt BP_lost_PELE1_topgo_enrichment.txt BP_lost_PELO1_Crass_topgo_enrichment.txt BP_lost_PHEI1_topgo_enrichment.txt BP_lost_PPUN2_topgo_enrichment.txt BP_lost_RVAR1_topgo_enrichment.txt BP_lost_TRLO1_topgo_enrichment.txt

# Rename files (in /mnt/netapp2/Store_csbyegim/Metazoa_OGs/constellatoR/losses/enriched_GO_terms_per_node)
rename "BP_lost_" "" *
rename "_topgo_enrichment" "" *

# Get the semantic similarity matrix
module load cesga/2020 python/3.9.9

source $STORE2/Dark_proteome/Isoforms_Semsim/pygosemsim_venv/bin/activate

WD=/mnt/netapp2/Store_csbyegim/Metazoa_OGs/constellatoR/gains
INPUT_FOLDER=$WD/enriched_GO_terms_per_node/

cd $WD

python $WD/../calculate_semsim_enriched_gos.py -i $INPUT_FOLDER -o $WD/gains_GO_enrich_nodes_SimMat.out
```

```python
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.colorbar
import pandas as pd

df_go = pd.read_csv("losses_GO_enrich_nodes_SimMat.out", sep="\t", index_col=0)
#Change order of columns/rows (to follow phylogeny)
desired_order = ["Tetrapoda","Rhabditida","Tardigrada","Onychophora","Arachnida","Myriapoda","Isopoda","Brachyura","Anomura","Hexapoda","Geoplanidae","Ellobiida","Stylomatophora","Crassiclitellata","Acteonemertidae"]  # Replace with actual row/column names
df_go = df_go.reindex(index=desired_order, columns=desired_order)

#For 2 triangular heatmaps
#sns.heatmap(df_go,annot=True,mask=matrix_go,cmap="BuPu")
#sns.heatmap(df_clusters,annot=True,mask=matrix_clusters,cmap="YlGnBu")

plt.figure(figsize=(10, 10))
sns.heatmap(df_go,annot=True,cmap="YlGnBu")
#plt.xticks(rotation=45)
#plt.yticks(rotation=45)
plt.savefig("GO_enrich_losses_convergence_nodes.svg",bbox_inches="tight")
```

![GO_enrich_losses_convergence_nodes.png](images/GO_enrich_losses_convergence_nodes.png)