#!/usr/bin/venv/python

'''
Script created by Gemma I. Martinez Redondo to calculate the gene losses per node in the species tree
Usage:
	python ogs_losses.py
'''

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
