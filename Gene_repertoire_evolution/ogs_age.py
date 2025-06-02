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
outtree="metazoa_ogs.nwk"
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
	OG_age[OG]=node

'''Check
count=0
for node in t.traverse():
	try:
		count=count+node.OG_num
	except AttributeError:
		continue
'''

#Don't understand why this is not working when loading after (something related to ete3??) pickle5
with open("OG_age.pickle","wb") as f:
	pickle.dump(OG_age,f,pickle.HIGHEST_PROTOCOL)

with open(outtree.split(".")[0]+".pickle","wb") as f:
        pickle.dump(t,f,pickle.HIGHEST_PROTOCOL)

t.write(format=8,outfile=outtree, features=["OG_num"])


'''
#Tree visualization (not working in CESGA)

from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace

def layout(node):
    if node.is_leaf():
        # Add node name to laef nodes
        N = AttrFace("name", fsize=14, fgcolor="black")
        faces.add_face_to_node(N, node, 0)
    if "OG_num" in node.features:
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"
        C = CircleFace(radius=node.OG_num, color="RoyalBlue", style="sphere")
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")

ts = TreeStyle()
ts.layout_fn=layout
ts.mode="c"
t.show(tree_style=ts)
'''
