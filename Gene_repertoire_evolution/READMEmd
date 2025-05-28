## 1. Description

.

## 2. Files in this repository

## 3. Workflow

### Gene gain analysis (age of origin of each OG - phylostratigraphy)

We calculated the age of OGs by inferring the MRCA of the species represented in each OG.

- Script for age estimation:
    
    ```python
    #!/usr/bin/venv/python
    
    '''
    Script created by Gemma I. Martinez Redondo to obtain the age of each orthogroup and add the number of genes appearing at eacn node in the species tree
    Usage:
    	python ogs_age.py -i ortholog_groups/ortholog_counts_per_species.stats.tsv -s metazoa_outgroups_sp_tree.nwk [-c species_conversion.txt] -o metazoa_ogs.nwk
    '''
    
    import argparse
    from ete3 import PhyloTree
    from itertools import compress
    import numpy as np
    import pickle
    
    def parser():
    	args=argparse.ArgumentParser(description="This program calculates the age of the SonicParanoid orthogroups and calculates the number of genes that have been gained at each node of the given species tree.")
    	args.add_argument('-i','--infile',required=True,help="Input file ortholog_counts_per_species.stats.tsv from SonicParanoid that contains the number of sequences from each species in each orthogroup.")
    	args.add_argument('-s','--speciestree',required=True,help="Species tree in newick format.")
    	args.add_argument('-c','--conversion',required=False,help="Tab-separated species conversion file. First column must contain name of species file in the orgholog_counts_per_species.stats.tsv. Second column contains the species as appearing in the species tree. Use when name of species in species tree and in the SonicParanoid output files do not match.")
    	args.add_argument('-o','--outtree',required=False,help="Output tree with number of OGs that appeared at each node mapped. If not provided 'out_tree.nwk' will be used by default.")
    	args=args.parse_args()
    	return args
    
    #Obtain parsed data
    og_file=parser().infile
    species_tree=parser().speciestree
    conversion_file=parser().conversion
    outtree=parser().outtree
    
    '''
    og_file="ortholog_groups/ortholog_counts_per_species.stats.tsv"
    species_tree="metazoa_outgroups_sp_tree.nwk"
    conversion_file="species_conversion.txt"
    outtree="metazoa_ogs_age.nwk"
    '''
    
    #Check that names from species tree match the ones from orthogroups file if conversion file not used. If not match, raise an error
    if not conversion_file:
    	with open(og_file,"rt") as og_count:
    		tree_sps=set(PhyloTree(species_tree).get_leaves())
    		og_sps=set(og_count.readline().strip().split("\t")[1:])
    		if len(tree_sps.symmetric_difference(og_sps)) !=0:
    			raise Exception("Species tree and output of SonicParanoid do not contain the same species. Please, run again including the correct conversion file.")
    		else:
    			species_conversion={key: value for key, value in zip(og_sps,tree_sps)}
    else:
    	species_conversion={}
    	with open(conversion_file,"rt") as conv:
    		line=conv.readline().strip()
    		while line:
    			sp_sonic,sp_tree=line.split("\t")
    			species_conversion[sp_sonic]=[sp_tree]
    			line=conv.readline().strip()
    
    if not outtree:
    	outtree="out_tree.nwk"
    
    #Obtain file saying which species are present at each OG (add all species where not 0 value). Then use ete to take the MRCA and assign age
    
    OG_species={}
    with open(og_file,"rt") as og_count:
    	line=og_count.readline().strip()
    	species=[]
    	for sp_file in line.split("\t")[1:-1]:
    		species.append(species_conversion[sp_file])
    	line=og_count.readline().strip()
    	species=np.array(species)
    	while line: #Convert line into True, False if values > 0 and multiply to species list to obtain species present
    		values=line.split("\t")
    		og,sp_count=values[0],values[1:-1]
    		sp_bool=[int(x)>0 for x in sp_count]
    		species_og=species[sp_bool].tolist()
    		OG_species[int(og)]=species_og
    		line=og_count.readline().strip()
    
    sp_per_og=[]
    for og in OG_species.keys():
    	sp_per_og.append(len(OG_species[og]))
    	OG_species[og]=[item for sublist in OG_species[og] for item in sublist] # In the original OG_species species were lists inside a list
    
    OG_age={}
    t=PhyloTree(species_tree, format=8)
    for OG in OG_species.keys():
    	species=OG_species[OG]
    	if len(species) !=1:
    		node=t.get_common_ancestor(species) #Obtain node where gene appeared
    	else:
    		node=t&species[0]
    	try:
    		node.OG_num=node.OG_num+1 #Add number of genes that are in each node
    	except AttributeError:
    		node.add_feature("OG_num",1)
    	OG_age[OG]=node.name
    
    with open("OG_age.pickle","wb") as f:
    	pickle.dump(OG_age,f,pickle.HIGHEST_PROTOCOL)
    
    with open("OG_age.txt","wt") as f:
    	for OG, age in OG_age.items():
    		f.write(f"{OG}\t{age}\n")
    
    with open(outtree.split(".")[0]+".pickle","wb") as f:
            pickle.dump(t,f,pickle.HIGHEST_PROTOCOL)
    
    t.write(format=8,outfile=outtree, features=["OG_num"])
    ```
    

```bash
python ogs_age.py -i ortholog_groups/ortholog_counts_per_species.stats.tsv -s metazoa_sp_tree_node_names.nwk -c species_conversion.txt -o metazoa_ogs_origin.nwk
```

Remember that the script does not count OGs at ROOT

![OGs_age_origin_no_root.png](images/OGs_age_origin_no_root.png)

### Protostomia gains plot for paper

To better show in the paper’s figure the gains in the Protostomia node, we selected one species (*Norana najaformis -*NNAJ1-) and plotted the age of gain of all the OGs in that species.

- Code to plot phylostratigraphy of OGs in a given species
    
    ```python
    #!/usr/bin/env/python
    '''
    Script created by Gemma I. Martinez-Redondo to create a phylostratigraphy plot for a given species (4 letter [+ number] code).
    Usage:
            python get_phylostratigraphy_plot.py species
    '''
    
    import pandas as pd
    import sys
    import os.path
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Set species of interest
    species="NNAJ1" #sys.argv[1]
    
    # Define the input and output file paths
    wdpath="/mnt/netapp2/Store_csbyegim/Metazoa_OGs/"
    ogs_age=wdpath+"OG_age.txt"
    sp_ogs_table_file=wdpath+"ortholog_groups/flat.ortholog_groups.tsv"
    species_conversion_file=wdpath+"species_conversion.txt"
    
    # Get species conversion
    sp_conversion={}
    with open(species_conversion_file,"rt") as f:
    	line=f.readline().strip()
    	while line:
    		sp_file,code=line.split("\t")
    		sp_conversion[sp_file]=code
    		line=f.readline().strip()
    
    code_to_num_sp_conversion={v: k for k, v in sp_conversion.items()}
    
    # Get OG-species information
    sp_ogs_table = pd.read_csv(sp_ogs_table_file, sep="\t", usecols=['group_id', code_to_num_sp_conversion[species]])[lambda d: d[code_to_num_sp_conversion[species]] != '*']
    #sp_ogs_table=sp_ogs_table.set_index('group_id')
    
    # Get information on OG-age
    og_age = pd.read_csv(ogs_age, sep="\t", header=None)
    og_age.columns = ["group_id","Phylostratum"]
    
    # Combine age information with OGs to get number of OGs in the given species for each phylostratum
    sp_ogs_age = pd.merge(sp_ogs_table,og_age, on="group_id")
    age_ogs_sp = sp_ogs_age.groupby('Phylostratum')['group_id'].nunique().reset_index()
    # Separate numeric and non-numeric Phylostratum values
    numeric_mask = pd.to_numeric(age_ogs_sp['Phylostratum'], errors='coerce').notna()
    # Convert numeric Phylostratum to int
    numeric_part = age_ogs_sp[numeric_mask].copy()
    numeric_part['Phylostratum'] = numeric_part['Phylostratum'].astype(int)
    numeric_part = numeric_part.sort_values(by='Phylostratum')
    # Non-numeric part (goes last)
    non_numeric_part = age_ogs_sp[~numeric_mask]
    # Combine: numeric first, then non-numeric
    age_ogs_sp_sorted = pd.concat([numeric_part, non_numeric_part], ignore_index=True)
    # Convert numeric count to int
    age_ogs_sp_sorted['group_id'] = age_ogs_sp_sorted['group_id'].astype(int)
    # Convert Phylostratum to str
    age_ogs_sp_sorted['Phylostratum'] = age_ogs_sp_sorted['Phylostratum'].astype(str)
    
    # Add extra info to phylostrata". Only for NNAJ1
    replacements = {
    "974": "Opisthokonta - N974",
    "976": "Holozoa - N976",
    "978": "Filozoa - N978",
    "981": "Coanozoa - N981",
    "983": "Metazoa - N983",
    "984": "Myriazoa - N984",
    "986": "Parahoxozoa - N986",
    "989": "Planulozoa - N989",
    "992": "Bilateria - N992",
    "996": "Nephrozoa - N996",
    "1002": "Protostomia - N1002",
    "1008": "Spiralia - N1008",
    "1017": "Platytrochozoa - N1017",
    "1032": "Lophotrochozoa - N1032",
    "1103": "Annelida - N1103",
    "1286": "Pleistoannelida - N1286",
    "1632": "Clitellata - N1632",
    "1908": "Crassiclitellata - N1908",
    "1925": "Lumbricina - N1925",
    "1938": "Hormogastirdae - N1938",
    "NNAJ1": "Norana najaformis - NNAJ1",
    "1048": "N1048",
    "1071": "N1071",
    "1141": "N1141",
    "1181": "N1181",
    "1229": "N1229",
    "1343": "N1343",
    "1403": "N1403",
    "1463": "N1463",
    "1517": "N1517",
    "1572": "N1572",
    "1684": "N1684",
    "1737": "N1737",
    "1777": "N1777",
    "1813": "N1813",
    "1840": "N1840",
    "1865": "N1865",
    "1887": "N1887"
    }
    
    age_ogs_sp_sorted['Phylostratum'] = age_ogs_sp_sorted['Phylostratum'].replace(replacements)
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    fig.subplots_adjust(hspace=0.05)  # adjust space between axes
    
    # plot the same data on both axes
    ax1.bar(age_ogs_sp_sorted['Phylostratum'], age_ogs_sp_sorted['group_id'], color='skyblue')
    ax2.bar(age_ogs_sp_sorted['Phylostratum'], age_ogs_sp_sorted['group_id'], color='skyblue')
    
    # zoom-in / limit the view to different portions of the data
    ax1.set_ylim(4700, 5000)  # outliers only
    ax2.set_ylim(0, 1100)  # most of the data
    
    # hide the spines between ax and ax2
    ax1.spines.bottom.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    
    # Now, let's turn towards the cut-out slanted lines.
    # We create line objects in axes coordinates, in which (0,0), (0,1),
    # (1,0), and (1,1) are the four corners of the axes.
    # The slanted lines themselves are markers at those locations, such that the
    # lines keep their angle and position, independent of the axes size or scale
    # Finally, we need to disable clipping.
    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                  linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
    
    # Custom x-tick labels: add "N" before each
    #xtick_labels = ['N' + str(ps) for ps in age_ogs_sp_sorted['Phylostratum']]
    xtick_labels = list(age_ogs_sp_sorted['Phylostratum'])
    ax2.set_xticks(range(len(age_ogs_sp_sorted['Phylostratum'])))
    ax2.set_xticklabels(xtick_labels, rotation=45)
    
    for label in ax2.get_xticklabels():
        label.set_horizontalalignment('right')
    
    # Labeling
    ax2.set_xlabel('Phylostratum')  # x-axis label for ax2
    ax2.set_ylabel('Number of OGs gained')  # y-axis label for ax2
    ax1.set_title('OG Count per Phylostratum')
    plt.tight_layout()
    
    # Save as SVG
    plt.savefig(species+"phylostratum_OGs_gained.svg", format='svg')
    
    #Old version
    # Create bar plot
    #plt.figure(figsize=(10, 6))
    #plt.bar(age_ogs_sp_sorted['Phylostratum'], age_ogs_sp_sorted['group_id'], color='skyblue')
    
    # Custom x-tick labels: add "N" before each
    #xtick_labels = ['N' + str(ps) for ps in age_ogs_sp_sorted['Phylostratum']]
    #plt.xticks(ticks=range(len(age_ogs_sp_sorted['Phylostratum'])), labels=xtick_labels, rotation=45)
    
    # Labeling
    #plt.xlabel('Phylostratum')
    #plt.ylabel('Number of OGs gained')
    #plt.title('OG Count per Phylostratum')
    #plt.tight_layout()
    
    # Save as SVG
    #plt.savefig(species+"phylostratum_OGs_gained.svg", format='svg')
    ```
    

![NNAJ1phylostratum_OGs_gained.svg](images/NNAJ1phylostratum_OGs_gained.svg)

### Gene loss

The idea is to look at the putative node where genes have been lost based on the presence/absence patterns of OGs. That way, if none of the children species of a specific node has the gene, we can assume a gene loss at that node.

- Code for reading species tree and OGs-species information
    
    ```python
    from ete3 import PhyloTree
    from itertools import compress
    import numpy as np
    import pickle
    
    og_file="ortholog_groups/ortholog_counts_per_species.stats.tsv"
    species_tree="metazoa_sp_tree_node_names.nwk"
    conversion_file="species_conversion.txt"
    
    #Check that names from species tree match the ones from orthogroups file if conversion file not used. If not match, raise an error
    if not conversion_file:
    	with open(og_file,"rt") as og_count:
    		tree_sps=set(PhyloTree(species_tree).get_leaves())
    		og_sps=set(og_count.readline().strip().split("\t")[1:])
    		if len(tree_sps.symmetric_difference(og_sps)) !=0:
    			raise Exception("Species tree and output of SonicParanoid do not contain the same species. Please, run again including the correct conversion file.")
    		else:
    			species_conversion={key: value for key, value in zip(og_sps,tree_sps)}
    else:
    	species_conversion={}
    	with open(conversion_file,"rt") as conv:
    		line=conv.readline().strip()
    		while line:
    			sp_sonic,sp_tree=line.split("\t")
    			species_conversion[sp_sonic]=[sp_tree]
    			line=conv.readline().strip()
    			
    #Obtain file saying which species are present at each OG (add all species where not 0 value).
    
    OG_species={}
    with open(og_file,"rt") as og_count:
    	line=og_count.readline().strip()
    	species=[]
    	for sp_file in line.split("\t")[1:-1]:
    		species.append(species_conversion[sp_file])
    	line=og_count.readline().strip()
    	species=np.array(species)
    	while line: #Convert line into True, False if values > 0 and multiply to species list to obtain species present
    		values=line.split("\t")
    		og,sp_count=values[0],values[1:-1]
    		sp_bool=[int(x)>0 for x in sp_count]
    		species_og=species[sp_bool].tolist()
    		OG_species[int(og)]=species_og
    		line=og_count.readline().strip()
    
    #Get number of species per OG and fix structure of OG_species
    sp_per_og=[]
    for og in OG_species.keys():
    	sp_per_og.append(len(OG_species[og]))
    	OG_species[og]=[item for sublist in OG_species[og] for item in sublist] # In the original OG_species species were lists inside a list
    
    #Load species tree
    t=PhyloTree(species_tree, format=8)
    ```
    

The code used to calculate them is based on this [script](https://github.com/fmarletaz/comp_genomics/blob/master/ortho_stats_dino.ipynb):

- Code for calculating gene losses of one OG
    
    ```python
    def get_losses_OG(OG_species,sp_tree):
    	if len(OG_species)>1:
    		phtyp=sp_tree.get_common_ancestor(OG_species)
    		ndesc=len([l.name for l in phtyp.get_leaves()])
    		if ndesc == len(OG_species):
    			return None
    		for leaf in sp_tree:
    			pv='1' if leaf.name in OG_species else '0'
    			leaf.add_features(presence=pv)
    		lost=[node.name for node in phtyp.get_monophyletic(values=['0'], target_attr="presence")]
    		return ndesc,phtyp.name,lost
    	else:
    		return None
    ```
    

The Jupyter Notebook this code is based upon also contains a part to calculate gene expansions per species. As stated in the paper (https://www.nature.com/articles/s41559-020-01327-6#Sec4), “Gene expansions were computed for each species using a hypergeometric test against the median gene number per species for a given family.”. However, several problems with that approach preclude us from using it: 1) It is calculated at the species level, meaning that we will not get information on the lineage. 2) Even if we modified the code to calculate the enrichment at the lineage level, we would have the problem of the different taxonomic representation. If in one lineage we have more species than in others, the expected median number of genes will be biased towards the lineage that contains more species. This will be true for both species and lineage gene expansions. 3) Depending on the gene family (less taxonomic representation), how do we know if it has been an expansion on the species/lineage of interest or a contraction in the others?

- Code for getting losses for all OGs and plot on the tree (name of root added manually)
    
    ```python
    #Get all losses for all OGs
    OG_losses = [get_losses_OG(OG_species[OG],t) for OG in OG_species]
    with open("OG_losses.pickle","wb") as f:
            pickle.dump(OG_losses,f,pickle.HIGHEST_PROTOCOL)
    
    #Get at each node of the species tree gene losses
    OG_losses_per_node = {}
    for i in range(len(OG_losses)):
    	if OG_losses[i] is None:
    		continue
    	else:
    		for node_lost in OG_losses[i][2]: #Name of node
    			node = t.search_nodes(name=node_lost)[0] #Node object obtained from node name
    			try:
    				node.num_loss = node.num_loss+1
    			except AttributeError:
    				node.add_feature("num_loss",1)
    			OG_losses_per_node[node_lost]=node.num_loss
    
    #Save tree with losses
    t.write(format=8,outfile="metazoa_sp_tree_gene_loss.nwk", features=["num_loss"])
    ```
    

![OGs_losses.png](images/OGs_losses.png)

Let’s check OGs gain/loss at the terrestrialization nodes (for nodes corresponding to just one species, let’s include also the genes not assigned to OGs -parenthesis-):

- Code for counting genes not assigned to OGs in those species
    
    ```bash
    cd /mnt/netapp2/Store_csbyegim/Metazoa_OGs/ortholog_groups
    ```
    
    ```python
    unassigned_genes_file="not_assigned_genes.ortholog_groups.tsv"
    sp_conversion_file="../species_conversion.txt"
    terr_species=["TRLO1","ETES1","RVAR1","PPUN2","ANGR1","OVER1","PELE1","CYCO1","ECRY2","HRPE1","PHEI1"]
    
    #Read unassigned genes
    unassigned_genes = {}
    with open(unassigned_genes_file, "rt") as f:
        while True:
            line = f.readline()  # Read a line from the file
            if not line:  # Break the loop if end of file is reached
                break
            line = line.strip()  # Remove leading and trailing whitespace
            if not line:  # Skip empty lines
                continue
            if line[0] == "#":  # Start of a new species section
                sp = line[1:]
                unassigned_genes[sp] = []
            else:  # Add gene to the current species
                unassigned_genes[sp].append(line)
    
    #Read conversion file
    sp_conversion={}
    with open(sp_conversion_file,"rt") as f:
    	line=f.readline().strip()
    	while line:
    		sp_file,sp_code=line.split("\t")
    		sp_conversion[sp_code]=sp_file
    		line=f.readline().strip()
    
    #Get list of unasigned genes for each of the species of interest
    for sp in terr_species:
    	print(sp)
    	len(unassigned_genes[sp_conversion[sp]])
    ```
    
- Results with Old SonicParanoid OGs
    
    
    | Terrestrialization event | Node | # Gains | # Loss |
    | --- | --- | --- | --- |
    | Hexapoda | 1353 | 296 | 2012 |
    | Anomura (Coenobitidae) | 1815 | 536 | 4579 |
    | Brachyura (Gecarcinidae) | 1913 | 84 | 4738 |
    | Isopoda | 1592 | 384 | 11164 |
    | TRLO1 (Amphipoda) | TRLO1 | 9 (+6317) | 5790 |
    | Myriapoda | 1117 | 493 | 60067 |
    | Arachnida | 1119 | 645 | 3081 |
    | Onychophora | 1057 | 2545 | 89347 |
    | ETES1 (tardigrade) | ETES1 (tardigrade) | 4 (+8778) | 1181 |
    | RVAR1 (tardigrade) | RVAR1 (tardigrade) | 19 (+7500) | 1072 |
    | Macrobiotidae + Richtersiidae | 1082 | 295 | 962 |
    | PPUN2 (nematode) | PPUN2 | 11 (+7268) | 1944 |
    | ANGR1 (nematode) | ANGR1 | 7 (+8567) | 3344 |
    | Rhabditida | 1236 | 534 | 2980 |
    | Stylomatophora | 1801 | 167 | 7202 |
    | Ellobiida (mollusca) | 1830 | 1011 | 7065 |
    | Onchidiidae (systelommatophora) | OVER1 | 11 (+3262) | 4297 |
    | Pomatiidae (caenogastropoda) | PELE1 | 10 (+4685) | 4501 |
    | DRAW1+PELO1+Crassi | 1860 | 319 | 2611 |
    | PELO1+Crass | 1882 | 236 | 1877 |
    | Crassiclitellata | 1903 | 552 | 2284 |
    | CYCO1 (leech) | CYCO1 (leech) | 17 (+2141) | 6640 |
    | ECRY2 (enchytreid) | ECRY2 (enchytreid) | 21 (+1975) | 4118 |
    | HRPE1 (annelid) | HRPE1 | 19 (+4402) | 3710 |
    | PHEI1 (annelid) | PHEI1 (annelid) | 7 (+3249) | 12304 |
    | Acteonemertidae (nemertea) | 1391 | 445 | 3649 |
    | Geoplanidae | 1612 | 348 | 2541 |
    | Tetrapoda | 1143 | 231 | 1189 |

| Terrestrialization event | Node | # Gains | # Loss |
| --- | --- | --- | --- |
| Hexapoda | 1357 | 305 | 1959 |
| Anomura (Coenobitidae) | 1820 | 537 | 4516 |
| Brachyura (Gecarcinidae) | 1918 | 85 | 4714 |
| Isopoda | 1597 | 380 | 7703 |
| TRLO1 (Amphipoda) | TRLO1 | 16 (+6319) | 5811 |
| Myriapoda | 1119 | 461 | 61232 |
| Arachnida | 1121 | 630 | 3018 |
| Onychophora | 1060 | 2341 | 88650 |
| ETES1 (tardigrade) | ETES1 (tardigrade) | 3 (+8775) | 1663 |
| RVAR1 (tardigrade) | RVAR1 (tardigrade) | 18 (+7486) | 1308 |
| Macrobiotidae + Richtersiidae | 1084 | 290 | 1094 |
| PPUN2 (nematode) | PPUN2 | 16 (+7187) | 2181 |
| ANGR1 (nematode) | ANGR1 | 9 (+8520) | 4067 |
| Rhabditida | 1239 | 525 | 2874 |
| Stylomatophora | 1806 | 158 | 6954 |
| Ellobiida (mollusca) | 1835 | 1026 | 6779 |
| Onchidiidae (systelommatophora) | OVER1 | 10 (+3277) | 4649 |
| Pomatiidae (caenogastropoda) | PELE1 | 3 (+4690) | 4727 |
| DRAW1+PELO1+Crassi | 1865 | 323 | 2409 |
| PELO1+Crass | 1887 | 278 | 1532 |
| Crassiclitellata | 1908 | 527 | 1758 |
| CYCO1 (leech) | CYCO1 (leech) | 18 (+1912) | 9939 |
| ECRY2 (enchytreid) | ECRY2 (enchytreid) | 30 (+1976) | 4670 |
| HRPE1 (annelid) | HRPE1 | 16 (+4391) | 4229 |
| PHEI1 (annelid) | PHEI1 (annelid) | 8 (+3235) | 13316 |
| Acteonemertidae (nemertea) | 1395 | 197 | 3622 |
| Geoplanidae | 1617 | 352 | 2806 |
| Tetrapoda | 1145 | 232 | 1128 |