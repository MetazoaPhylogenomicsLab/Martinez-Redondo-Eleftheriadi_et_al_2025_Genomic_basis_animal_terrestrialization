#!/usr/bin/env/python

#Read headers of file
with open("filtered_transposed_Species_OGs_habitat_phylum.tsv/00.part","rt") as f:
	data_columns=set(f.readline().strip().split("\t"))

#Read list of OGs to remove
feats_to_remove=[]
with open("features_eq0_first_it.txt","rt") as f:
	line=f.readline().strip()
	while line:
		feats_to_remove.append(line)
		line=f.readline().strip()

#Get columns to keep
columns_to_keep=data_columns-set(feats_to_remove)

#Add 5 (columns not OGs) to the number of OGs to keep and
OGs_final=[]
for col in columns_to_keep:
	if col.isnumeric():
		OGs_final.append(str(int(col)+5))
	else:
		continue

#Print final number of columns
print(",".join(OGs_final))
