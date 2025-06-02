# Phylostratigraphy of Hub-containing orthogroups

---

Keeping the 17 species of interest out of 973 the resulted matrix:

SonicParanoid_invertebrates_experiments.mod.txt (check our FigShare)

Using the following scripts we mapped the hub genes of the statistically significant modules per experiment for each species in the OGs inferred using SonicParanoid.

```python
import pandas as pd
import os
import glob
from collections import defaultdict

SPECIES = {
    "Annelida": ["Arenicola_marina", "Eisenia_andrei", "Hirudo_medicinalis"],
    "Arthropoda": ["Porcellio_laevis", "Ligia_oceanica"],
    "Mollusca": ["Phorcus_turbinatus", "Physella_acuta", "Theba_pisana", "Siphonaria_pectinata"],
    "Nemertea": ["Leptonemertes_chalicophora", "Tetrastemma_longissimum", "Tetrastemma_melanocephalum"],
    "Nematoda": ["Caenorhabditis_elegans", "Litoditis_marina"],
    "Onychophora": ["Peripatoides_aurorbis"],
    "Platyhelminthes": ["Obama_nungara", "Schmidtea_mediterranea"]
}

SPECIES_CODES = {
    "Arenicola_marina": 'COCO', "Eisenia_andrei": 'EAND', "Hirudo_medicinalis": 'HMED', # Annelida
    "Porcellio_laevis": 'PLAE', "Ligia_oceanica": 'MISO',  # Arthropoda
    "Phorcus_turbinatus": 'PTUR', "Physella_acuta": 'PACU', "Theba_pisana": 'TPIS', "Siphonaria_pectinata": 'SPEC', # Mollusca
    "Leptonemertes_chalicophora": 'LEPN', "Tetrastemma_longissimum": 'TLON', "Tetrastemma_melanocephalum": 'TMEL', # Nemertea
    "Caenorhabditis_elegans": 'CELE', "Litoditis_marina": 'LMAR', # Nematoda
    "Peripatoides_aurorbis": 'PEAU', # Onychophora
    "Obama_nungara": 'ONUN', "Schmidtea_mediterranea": 'SMED'  # Platyhelminthes
}

BasePath = "./"
OutputDir = os.path.join(BasePath,"Merged_Experiments")
#os.makedirs(OutputDir, exist_ok=True)  # Create output directory if it doesn't exist

experiment_files = {}

for phylum, species_list in SPECIES.items():
    for species in species_list:
        species_path = os.path.join(BasePath, phylum, species, "WGCNA/3.Relate_modules_to_external_traits/6.HubGenes")

        if os.path.exists(species_path):
            hub_files = glob.glob(os.path.join(species_path, "*_HubGenes.txt"))

            for file_path in hub_files:
                filename = os.path.basename(file_path)
                experiment_id = filename.replace("_HubGenes.txt", "")

                if experiment_id not in experiment_files:
                    experiment_files[experiment_id] = []

                experiment_files[experiment_id].append((species, file_path))
               

for experiment_id, files in experiment_files.items():
    merged_data = []

    for species, file_path in files:
        try:
            experiment_species_df = pd.read_csv(file_path, sep="\t")
            experiment_species_df["Species"] = SPECIES_CODES[species]  # Add species column
            experiment_species_df["Experiment"] = experiment_id  # Add experiment column
            experiment_species_df = experiment_species_df.drop_duplicates(subset = ["Gene"])
            merged_data.append(experiment_species_df)

        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    if merged_data:
        merged_df = pd.concat(merged_data, ignore_index=True)

        if not merged_df.empty:  # Check if the dataframe is not empty
            merged_df = merged_df.drop_duplicates(subset=['Protein', 'Species'], keep='first') #keep first occurrence of a gene for a species.

        output_file = os.path.join(OutputDir, f"{experiment_id}_Merged.txt")
        merged_df.to_csv(output_file, sep="\t", index=False)
        print(f"Merged file saved: {output_file}")
    else:
        print(f"No data to merge for experiment {experiment_id}") # Handle cases where there is no data

# Function to find group IDs for hub proteins
def find_OG_per_species_per_experiment(hub_proteins, OGs_df):
    
###############################################################################
## Function that returns 2 dataframes summarised per experiment 
## (meaning duplications in OGs):
##  1. match_hubs2OGs_df with info per: 
##	          'Experiment', 'Node', 'Species', 'OG',  'Proteins'
##  2. summarized_hubs2OGs_df same as match_hubs2OGs_df but 
##       adding up the number of OGs and Proteins per experiment (groupby Node)
###############################################################################
    results = []
    
    for _, row in OGs_df.iterrows():
        group_id = row['group_id']
        node = row['Node']
        OG_proteins = [protein for protein in row[species_code].split(',') if '*' not in protein] if pd.notna(row[species_code]) else []
        
        common_proteins = list(set(OG_proteins) & set(hub_proteins))
        if common_proteins:
            results.append({
                'Experiment': experiment_id ,
                'Node': node,
                'Species': species_code,
                'OG': group_id ,
                'Proteins': ', '.join(common_proteins)
            })
    
    match_hubs2OGs_df = pd.DataFrame(results, columns=['Experiment', 'Node', 'Species', 'OG',  'Proteins']) # A "flat" df without merging any info
    
    summarized_hubs2OGs_df = match_hubs2OGs_df.groupby(['Node']).agg({
        'Experiment': lambda x: ', '.join(set(x)),
        'Species': lambda x: ', '.join(set(x)),
        'OG': lambda x: ', '.join(map(str, set(x))),  
        'Proteins': lambda x: ', '.join(map(str,set(x))), 
        }).reset_index()

        # Add summary counts
    summarized_hubs2OGs_df['Sum_of_OGs'] = summarized_hubs2OGs_df['OG'].apply(lambda x: len(x.split(',')))
    summarized_hubs2OGs_df['Sum_of_proteins'] = summarized_hubs2OGs_df['Proteins'].apply(lambda x: len(x.split(',')))
    
    return match_hubs2OGs_df, summarized_hubs2OGs_df        
                
        
# Define paths
SONIC_PARANOID_PATH = BasePath + "/SonicParanoid/"
OGs_FILE = os.path.join(SONIC_PARANOID_PATH, "SonicParanoid_invertebrates_experiments.mod.txt")

# Read the OGs data
OGs_df = pd.read_csv(OGs_FILE, sep="\t")
MergedDir = OutputDir

# Dictionary to store species as keys and list of proteins as values
species_experiment_dict = {}

# Process each merged experiment file
for file_name in os.listdir(MergedDir):
    if file_name.endswith("_Merged.txt"):
        experiment_id = file_name.split("_")[0]  # Extract experiment ID
        print(f"Processing experiment: {experiment_id}")
        
        file_path = os.path.join(MergedDir, file_name)
        try:
            # Read the file
            df = pd.read_csv(file_path, sep="\t")

            # Ensure there are enough columns
            if df.shape[1] >= 5:
                # Selecting 3rd (Protein), 4th (Species), and 5th (Experiment) columns - Removing transcript and module columns 
                hubs_per_experiment = df.iloc[:, [2, 3, 4]]  
                
                # Drop rows where the protein column (3rd column) is NaN
                hubs_per_experiment = hubs_per_experiment.dropna(subset=[hubs_per_experiment.columns[0]])  

                # Iterate through rows and populate the dictionary
                for row in hubs_per_experiment.itertuples(index=False, name=None):
                    protein, species, _ = row  # Unpack first two columns, ignore the third
                
                    if species not in species_experiment_dict:
                        species_experiment_dict[species] = {}
                   
                    if experiment_id not in species_experiment_dict[species]:
                        species_experiment_dict[species][experiment_id] = set()  # Use set to avoid duplicates

                    # Add the protein to the set
                    species_experiment_dict[species][experiment_id].add(protein)

            else:
                print(f"Skipping {file_name}: Not enough columns")
        except Exception as e:
            print(f"Error processing {file_name}: {e}")

# Convert sets to lists for final output
species_experiment_dict = {
    species: {exp: list(proteins) for exp, proteins in experiments.items()}
    for species, experiments in species_experiment_dict.items()
}

all_match_dfs = []
all_summarized_dfs = []

for species, experiments in species_experiment_dict.items():
    species_code = species 
    
    for experiment_id, hub_proteins in experiments.items():
        print(f"Processing: Species={species_code}, Experiment={experiment_id}, {len(hub_proteins)}")

       
        match_hubs2OGs_df, summarized_hubs2OGs_df = find_OG_per_species_per_experiment(hub_proteins, OGs_df)

        
        all_match_dfs.append(match_hubs2OGs_df)
        all_summarized_dfs.append(summarized_hubs2OGs_df)

final_match_df = pd.concat(all_match_dfs, ignore_index=True)
final_summarized_df = pd.concat(all_summarized_dfs, ignore_index=True)

final_match_df.to_csv(os.path.join(OutputDir, "AllSpecies_Match_Hubs2OGs_per_experiment.txt"), sep='\t', index=False)
final_summarized_df.to_csv(os.path.join(OutputDir, "AllSpecies_Match_Hubs2OGs_per_experiment_Summarised_OGs_and_Proteins_by_NODE.txt"), sep='\t', index=False)
  
```

We will now use the file generated from the above script AllSpecies_Match_Hubs2OGs_per_experiment_Summarised_OGs_and_Proteins_by_NODE.txt

to map to the nodes of interest the percentage of the hub-containing OGs per species and visualize the Phylostratigraphy in iTOL 

```python
import pandas as pd
import os
import glob
import numpy as np
from collections import defaultdict

SPECIES = {
    "Annelida": ["Arenicola_marina", "Eisenia_andrei", "Hirudo_medicinalis"],
    "Arthropoda": ["Porcellio_laevis", "Ligia_oceanica"],
    "Mollusca": ["Phorcus_turbinatus", "Physella_acuta", "Theba_pisana", "Siphonaria_pectinata"],
    "Nemertea": ["Leptonemertes_chalicophora", "Tetrastemma_longissimum", "Tetrastemma_melanocephalum"],
    "Nematoda": ["Caenorhabditis_elegans", "Litoditis_marina"],
    "Onychophora": ["Peripatoides_aurorbis"],
    "Platyhelminthes": ["Obama_nungara", "Schmidtea_mediterranea"]
}

SPECIES_CODES = {
    "Arenicola_marina": 'COCO', "Eisenia_andrei": 'EAND', "Hirudo_medicinalis": 'HMED', # Annelida
    "Porcellio_laevis": 'PLAE', "Ligia_oceanica": 'MISO',  # Arthropoda
    "Phorcus_turbinatus": 'PTUR', "Physella_acuta": 'PACU', "Theba_pisana": 'TPIS', "Siphonaria_pectinata": 'SPEC', # Mollusca
    "Leptonemertes_chalicophora": 'LEPN', "Tetrastemma_longissimum": 'TLON', "Tetrastemma_melanocephalum": 'TMEL', # Nemertea
    "Caenorhabditis_elegans": 'CELE', "Litoditis_marina": 'LMAR', # Nematoda
    "Peripatoides_aurorbis": 'PEAU', # Onychophora
    "Obama_nungara": 'ONUN', "Schmidtea_mediterranea": 'SMED'  # Platyhelminthes
}

BasePath = "./"
OutputDir = os.path.join(BasePath,"Merged_Experiments")

final_match_df = pd.read_csv(os.path.join(OutputDir, "AllSpecies_Match_Hubs2OGs_per_experiment.txt"), sep="\t")
final_summarized_df = pd.read_csv(os.path.join(OutputDir, "AllSpecies_Match_Hubs2OGs_per_experiment_Summarised_OGs_and_Proteins_by_NODE.txt"), sep="\t")
)

def summarize_species_info(final_summarized_df):
    
    ###########################################################################
    ###    Function that summarises OGs, Proteins and #OGs and #Proteins   ####
    ###               per Node per species per experiment                  ####
    ## Cols: Node Experiment Species OG Proteins  Sum_of_OGs Sum_of_proteins ##
    ###########################################################################
    
    # Select relevant columns
    final_summarized_df = final_summarized_df[['Node', 'Experiment', 'Species', 'OG', 'Proteins']]
     
     # Get all unique species
    unique_species = sorted(final_summarized_df['Species'].unique())
    print(unique_species)
    # Pivot table to group by Node and Experiment, ensuring all species are represented
    grouped_byNode_bySpecies_df = final_summarized_df.groupby(['Node', 'Species']).agg({
                                        'OG': lambda x: ', '.join(sorted(set(og.strip() for val in x for og in val.split(',') if og.strip()))),  
                                        'Proteins': lambda x: ', '.join(sorted(set([item.strip() for sublist in x for item in sublist.split(',')]))) 
                                        }).reset_index()
    
    grouped_byNode_bySpecies_df['Sum_of_OGs'] = grouped_byNode_bySpecies_df['OG'].apply(lambda x: len(x.split(',')))
    grouped_byNode_bySpecies_df['Sum_of_proteins'] = grouped_byNode_bySpecies_df['Proteins'].apply(lambda x: len(x.split(',')))
    
    
    pivot_df = grouped_byNode_bySpecies_df.pivot_table(
        index=['Node'], 
        columns='Species', 
        values='Sum_of_OGs', 
        aggfunc='sum', 
        fill_value=0  # Fill missing species with 0
    ).reset_index()
    
    # Ensure all species in unique_species are present in pivot_df
    for species in unique_species:
        if species not in pivot_df.columns:
            pivot_df[species] = 0  # Add missing species columns with 0
    
    # Create a new Species column with species names or '0' for missing ones
    pivot_df['Species'] = pivot_df[unique_species].apply(
        lambda row: ','.join([species if row[species] > 0 else '0' for species in unique_species]), 
        axis=1
    )
    
    # Create a new Sum_of_OGs column joining OG counts as strings
    pivot_df['Sum_of_OGs'] = pivot_df[unique_species].apply(lambda row: ','.join(map(str, row)), axis=1)
    
    # Keep only required columns
    final_df = pivot_df[['Node', 'Species', 'Sum_of_OGs']]
    
    
    # Sort by Experiment first, then by Node
    final_df = final_df.sort_values(by=['Node'])
    
    return unique_species, final_df
    
    
final_summarized_df_noCONT = final_summarized_df[final_summarized_df["Experiment"] != "CONT"]

# Run function for summarising
species_order, AllSpecies_nodes_df_noCONT = summarize_species_info(final_summarized_df_noCONT)
AllSpecies_nodes_df_noCONT.head()

species_OG_sums_noCONT = defaultdict(int)
# Iterate over each row
for _, row in AllSpecies_nodes_df_noCONT.iterrows():
    species_list = row['Species'].split(',')  # List of species
    og_counts = list(map(int, row['Sum_of_OGs'].split(',')))  # Convert OG counts to integers
    
    # Sum OGs for each species
    for species, og in zip(species_list, og_counts):
        if species != '0':  # Ignore '0' placeholders
            species_OG_sums_noCONT[species] += og

species_summary_df_noCONT = pd.DataFrame(list(species_OG_sums_noCONT.items()), columns=['Species', 'Total_OGs'])
# Sort species by total OGs in descending order
species_summary_df_noCONT = species_summary_df_noCONT.sort_values(by='Total_OGs', ascending=False).reset_index(drop=True)

# Keep only OG counts, merge TETR and normalise 
og_counts_df = AllSpecies_nodes_df_noCONT["Sum_of_OGs"]        
og_counts_df = og_counts_df.str.split(',', expand=True)

# Name columns with species name based on species order calculated before 
og_counts_df.columns = species_order

# Convert all to float otherwise merges str
og_counts_df = og_counts_df.astype(float)
# Create the new column TETR
og_counts_df['TETR'] = og_counts_df['TLON'] + og_counts_df['TMEL']

new_species_order = [col for col in species_order if col not in ['TLON', 'TMEL']] + ['TETR']

# Remove columns TLON and TMEL
og_counts_df.drop(columns=['TLON', 'TMEL'], inplace=True)

TETR_value = species_summary_df_noCONT[species_summary_df_noCONT['Species'].isin(['TMEL', 'TLON'])]['Total_OGs'].sum()
species_summary_df_noCONT_TETR = pd.concat([species_summary_df_noCONT, pd.DataFrame([{'Species': 'TETR', 'Total_OGs': TETR_value}])], ignore_index=True)
species_summary_df_noCONT_TETR = species_summary_df_noCONT_TETR[~species_summary_df_noCONT_TETR['Species'].isin(['TMEL', 'TLON'])].reset_index(drop=True)

# Now we need to normalise og_counts_df with species_summary_df_noCONT_TETR

# Let's convert back to dict the df with species total OG counts 
species_OG_totals_noCONT = dict(zip(species_summary_df_noCONT_TETR['Species'], species_summary_df_noCONT_TETR['Total_OGs']))

normalized_df = og_counts_df.copy()

for species in normalized_df.columns:
    if species in species_OG_totals_noCONT:
        normalized_df[species] = round((normalized_df[species] / species_OG_totals_noCONT[species]) * 100, 2)
    else:
        print(f"Warning: {species} not found in species_OG_totals")

# Now normalized_df contains all OGs as percentages of species totals
```

## Calculating proportions of hub-containing orthogroups per species

---

For the phylostratigraphy (Figure 3B) we kept the percentages of the hub-containing OGs per species in the following Phylostrata:

- All Phylostrata from Root (N974) till the Protostomia (N1002)
- The sum from Protostomia (N1002) till the phylum node of each species
- The phylum node
- The sum from the phylum node to the species node
- The species node

```python
#(together with the script above)
# Final table for the iTOL 
merged_df = pd.concat([AllSpecies_nodes_df_noCONT['Node'].reset_index(drop=True), normalized_df.reset_index(drop=True)], axis=1)
# Add an N before each node name (to merge with the tree node names)
merged_df['Node'] = merged_df['Node'].apply(lambda x: f'N{x}')    

# Define the Nodes to be extracted from the tree
Nodes = {
    "Opisthokonta": "N974",  # Root
    "Holozoa": "N976",
    "Filozoa": "N978",
    "Coanozoa": "N981",
    "Metazoa": "N983",
    "Myriazoa": "N984",
    "Parahoxozoa": "N986",
    "Planulozoa": "N989",
    "Bilateria": "N992",
    "Nephrozoa": "N996",
    "Protostomia": "N1002",
    "Phylum": ["N1039", "N1047", "N1060", "N1061", "N1095", "N1102", "N1103"],
    "Species": ["NLEPN", "NTMEL", "NTLON", "NEAND", "NHMED", "NCOCO", "NCELE",
                "NLMAR", "NTPIS", "NPACU", "NPTUR", "NSPEC1", "NPEAU", "NPLAE",
                "NMISO", "NONUN", "NSMED"]
}

Phylum2Species = { 
    "N1039": ["CELE", "LMAR"],
    "N1047":["PTUR", "SPEC1","TPIS", "PACU"],
    "N1060":[ "PEAU"],
    "N1061": ["PLAE","MISO"],
    "N1095": ["ONUN", "SMED"],
    "N1102": ["LEPN", "TMEL", "TLON"],
    "N1103": ["EAND", "HMED", "COCO"]
    }

REFERENCE_NODE = "N1002"

from Bio import Phylo

tree_path = os.path.join(BasePath, "SonicParanoid", "pruned_tree.nwk")

def find_named_clade(tree, name):
    for clade in tree.find_clades():
        if clade.name == name:
            return clade
    return None

def extract_path_between(tree, start_clade, end_clade):
    full_path = tree.get_path(end_clade)
    full_path.append(end_clade)

    if start_clade not in full_path:
        return None
    start_index = full_path.index(start_clade)
    return full_path[start_index:]

def main(newick_file):
    tree = Phylo.read(newick_file, "newick")
    ref_clade = find_named_clade(tree, REFERENCE_NODE)
    if ref_clade is None:
        raise ValueError(f"Reference node '{REFERENCE_NODE}' not found.")

    records = []

    for phylum_node, species_list in Phylum2Species.items():
        phylum_clade = find_named_clade(tree, phylum_node)
        if phylum_clade is None:
            print(f"Phylum node '{phylum_node}' not found.")
            continue

        for species_code in species_list:
            species_tip = find_named_clade(tree, species_code)
            if species_tip is None:
                print(f"Species tip '{species_code}' not found.")
                continue

            path1 = extract_path_between(tree, ref_clade, phylum_clade)
            path2 = extract_path_between(tree, phylum_clade, species_tip)

            # Remove both ends: start and end node from each path
            if path1:
                path1 = path1[1:-2]
            if path2:
                path2 = path2[1:-2]

            def path_to_names(path):
                return [clade.name if clade.name else f"[internal_{id(clade)}]" for clade in path] if path else None

            records.append({
                "Species": species_code,
                "Phylum_Node": phylum_node,
                "Path_N1002_to_Phylum": path_to_names(path1),
                "Path_Phylum_to_Species": path_to_names(path2)
            })

    df = pd.DataFrame(records)
    return df 
    
if __name__ == "__main__":
    newick_path = tree_path
    df_paths = main(newick_path)
    print(df_paths.head())  # Preview

# Identify TMEL and TLON rows
df_tmeltlon = df_paths[df_paths['Species'].isin(['TMEL', 'TLON'])]

# Check if they share the same paths (i.e., identical path lists)
#if not df_tmeltlon.empty and df_tmeltlon['Path_N1002_to_Phylum'].nunique() == 1 and df_tmeltlon['Path_Phylum_to_Species'].nunique() == 1:
    # Keep one of them and rename species to 'TETR'
merged_row = df_tmeltlon.iloc[0].copy()
merged_row['Species'] = 'TETR'

# Drop TMEL and TLON from the original df
df_paths = df_paths[~df_paths['Species'].isin(['TMEL', 'TLON'])]

# Append the merged row
df_paths = df_paths.append(merged_row, ignore_index=True)

# Function to calculate sum for each path segment
def calculate_path_sum(path, species, iTOL_df):
    
    # Fix species name if needed
    if species == "SPEC1":
        species = "SPEC"

    normalized_path = [f"N{node}" if not str(node).startswith("N") else node for node in path]
    total_sum = 0
    for node in normalized_path:
        if node in iTOL_df['Node'].values:
            total_sum += iTOL_df.loc[iTOL_df['Node'] == node, species].values[0]
    return total_sum

# Add sum columns for Path_N1002_to_Phylum and Path_Phylum_to_Species
df_paths['Path_N1002_to_Phylum_Sum'] = df_paths.apply(
    lambda row: calculate_path_sum(row['Path_N1002_to_Phylum'], row['Species'], iTOL_df), axis=1)

df_paths['Path_Phylum_to_Species_Sum'] = df_paths.apply(
    lambda row: calculate_path_sum(row['Path_Phylum_to_Species'], row['Species'], iTOL_df), axis=1)

print(df_paths[['Species', 'Path_N1002_to_Phylum_Sum', 'Path_Phylum_to_Species_Sum']])

#Create the iTOL barplot input 

# Flatten the node IDs
node_ids = [v for val in Nodes.values() for v in (val if isinstance(val, list) else [val])]

# Filter rows and drop columns
filtered_df = (
    iTOL_df[iTOL_df['Node'].isin(node_ids)]
    .drop(columns=['Position', 'Radius'])
    .reset_index(drop=True)
)

# Reverse map: node -> group
node_to_group = {}
for group, values in Nodes.items():
    if isinstance(values, list):
        for v in values:
            node_to_group[v] = group
    else:
        node_to_group[values] = group

# Set Node as index to prep for summing
nodes_info_df = filtered_df.set_index("Node")

# Collapse Phylum and Species
collapsed_rows = []
for group in ["Phylum", "Species"]:
    members = Nodes[group]
    existing_members = [m for m in members if m in nodes_info_df.index]
    
    if existing_members:
        summed = nodes_info_df.loc[existing_members].sum()
        summed.name = group
        collapsed_rows.append(summed)

# Safely drop only existing Phylum/Species nodes
all_nodes_to_remove = [m for m in (Nodes["Phylum"] + Nodes["Species"]) if m in nodes_info_df.index]
nodes_info_df = nodes_info_df.drop(all_nodes_to_remove)

# Append collapsed rows
nodes_info_df = pd.concat([nodes_info_df, pd.DataFrame(collapsed_rows)])

# Reorder according to the original Nodes dict
final_order = []
for group in Nodes:
    if group in ["Phylum", "Species"]:
        final_order.append(group)
    else:
        val = Nodes[group]
        if val in nodes_info_df.index:
            final_order.append(val)

# Keep only available rows in final_order
final_order = [n for n in final_order if n in nodes_info_df.index]
nodes_info_df = nodes_info_df.loc[final_order].reset_index().rename(columns={"index": "Node"})

# Assign Group column
nodes_info_df["Group"] = nodes_info_df["Node"].map(node_to_group)
nodes_info_df.loc[nodes_info_df["Node"] == "Phylum", "Group"] = "Phylum"
nodes_info_df.loc[nodes_info_df["Node"] == "Species", "Group"] = "Species"

# Remove the Group column if it exists
if "Group" in nodes_info_df.columns:
    nodes_info_df = nodes_info_df.drop(columns=["Group"])

# Transpose the DataFrame
transposed_df = nodes_info_df.set_index("Node").T

# Reverse the column order for the mirroring effect in iTOL
transposed_df = transposed_df[transposed_df.columns[::-1]]

df_paths['Species'] = df_paths['Species'].replace({'SPEC1': 'SPEC'})

# Remove columns if they already exist to prevent duplicates
#transposed_df = transposed_df.drop(columns=['Path_N1002_to_Phylum_Sum', 'Path_Phylum_to_Species_Sum'], errors='ignore')

# Set index on df_paths
df_paths_indexed = df_paths.set_index('Species')

# Merge on index
merged_df = transposed_df.merge(
    df_paths_indexed[['Path_N1002_to_Phylum_Sum', 'Path_Phylum_to_Species_Sum']],
    left_index=True, right_index=True,
    how='left'
)

# Reorder columns: insert before/after 'Phylum'
cols = [col for col in merged_df.columns if col not in ['Path_N1002_to_Phylum_Sum', 'Path_Phylum_to_Species_Sum']]

# Locate where to insert relative to 'Phylum'
phylum_idx = cols.index('Phylum')

# Insert the two new columns around 'Phylum'
new_order = (
    cols[:phylum_idx] +
    ['Path_N1002_to_Phylum_Sum', 'Phylum', 'Path_Phylum_to_Species_Sum'] +
    cols[phylum_idx+1:]
)

# Reorder columns
final_df = merged_df[new_order]
final_df.to_csv(os.path.join(OutputDir, "BarPlot_input_ALLSPECIES_normalized_100perc.csv"),sep=",", index=True)  # index=True keeps the species names as row index

```

The csv file is the following:

[BarPlot_input_ALLSPECIES_normalized_100perc.csv](Phylostratigraphy%20of%20Hub-containing%20orthogroups%20201f4c09425a8076acfbd5a0a5ea7d64/BarPlot_input_ALLSPECIES_normalized_100perc.csv)

The percentages of this file were incorporated in the iTOL annotation file to be plotted as a barplot: 

[All_Species_iTOL_tree_normalised_barplot.txt](Phylostratigraphy%20of%20Hub-containing%20orthogroups%20201f4c09425a8076acfbd5a0a5ea7d64/All_Species_iTOL_tree_normalised_barplot.txt)

<aside>
💡

Calculate the same using Hub Genes matched with Proteomics

</aside>

```python
proteomics_filepath = os.path.join(BasePath, "Proteomics")

# Remove Control samples 
final_match_df_noCONT = final_match_df[final_match_df["Experiment"]!="CONT"]

# Step 1: Load filtered_accessions files by species
valid_proteins_by_species = {}

for file in glob.glob(os.path.join(proteomics_filepath,"*filtered_accessions*.txt")):
    # Extract the base filename
    base_filename = os.path.basename(file)
    # Extract species name from the base filename
    species = base_filename.split("_filtered_accessions")[0]
    with open(file, 'r') as f:
        proteins = set(line.strip() for line in f if line.strip())
        valid_proteins_by_species.setdefault(species, set()).update(proteins)

# Step 2: Define matching function
def get_matching_proteins(row):
    species = row['Species']
    protein_list = [p.strip() for p in str(row['Proteins']).split(',')]
    valid_proteins = valid_proteins_by_species.get(species, set())
    matched = [p for p in protein_list if p in valid_proteins]
    return ','.join(matched) if matched else None

# Step 3: Apply matching function and filter rows
final_match_df_noCONT['Proteomics'] = final_match_df_noCONT.apply(get_matching_proteins, axis=1)
final_match_df_noCONT = final_match_df_noCONT[final_match_df_noCONT['Proteomics'].notna()]

def summarize_species_info_proteomics(final_summarized_df):
    """
    ###########################################################################
    ###             Same as before but using column "Proteomics"           ####
    ###    Function that summarises OGs, Proteins and #OGs and #Proteins   ####
    ###               per Node per species per experiment                  ####
    ## Cols: Node Experiment Species OG Proteins  Sum_of_OGs Sum_of_proteins ##
    ###########################################################################
    
    """

    # Select relevant columns
    final_summarized_df = final_summarized_df[['Node', 'Experiment', 'Species', 'OG', 'Proteomics']]

    # Ensure 'OG' and 'Proteomics' columns are of string type
    final_summarized_df['OG'] = final_summarized_df['OG'].astype(str)
    final_summarized_df['Proteomics'] = final_summarized_df['Proteomics'].astype(str)

    # Get all unique species
    unique_species = sorted(final_summarized_df['Species'].unique())

    # Group by Node and Species, aggregating OGs and Proteomics
    grouped_byNode_bySpecies_df = final_summarized_df.groupby(['Node', 'Species']).agg({
        'OG': lambda x: ','.join(sorted(set(
            item.strip() for sublist in x for item in str(sublist).split(',')
        ))),
        'Proteomics': lambda x: ','.join(sorted(set(
            item.strip() for sublist in x for item in str(sublist).split(',')
        )))
    }).reset_index()

    # Calculate counts
    grouped_byNode_bySpecies_df['Sum_of_OGs'] = grouped_byNode_bySpecies_df['OG'].apply(
        lambda x: len(x.split(',')) if x else 0
    )
    grouped_byNode_bySpecies_df['Sum_of_proteins'] = grouped_byNode_bySpecies_df['Proteomics'].apply(
        lambda x: len(x.split(',')) if x else 0
    )

    # Pivot table to have species as columns
    pivot_df = grouped_byNode_bySpecies_df.pivot_table(
        index='Node',
        columns='Species',
        values='Sum_of_OGs',
        aggfunc='sum',
        fill_value=0
    ).reset_index()

    # Ensure all species are present in pivot_df
    for species in unique_species:
        if species not in pivot_df.columns:
            pivot_df[species] = 0

    # Create 'Species' column indicating presence or absence
    pivot_df['Species'] = pivot_df[unique_species].apply(
        lambda row: ','.join([species if row[species] > 0 else '0' for species in unique_species]),
        axis=1
    )

    # Create 'Sum_of_OGs' column with counts per species
    pivot_df['Sum_of_OGs'] = pivot_df[unique_species].apply(
        lambda row: ','.join(map(str, row)),
        axis=1
    )

    # Final DataFrame with required columns
    final_df = pivot_df[['Node', 'Species', 'Sum_of_OGs']]

    # Sort by Node
    final_df = final_df.sort_values(by='Node')

    return unique_species, final_df

# Run function
species_order_proteomics, AllSpecies_nodes_df_noCONT_proteomics = summarize_species_info_proteomics(final_match_df_noCONT)

species_OG_sums_noCONT_proteomics = defaultdict(int)

for _, row in AllSpecies_nodes_df_noCONT_proteomics.iterrows():
    species_list = row['Species'].split(',')  # List of species
    og_counts = list(map(int, row['Sum_of_OGs'].split(',')))  # Convert OG counts to integers
    
    # Sum OGs for each species
    for species, og in zip(species_list, og_counts):
        if species != '0':  # Ignore '0' placeholders
            species_OG_sums_noCONT_proteomics[species] += og

# Keep only OG counts and normalise 
og_counts_df = AllSpecies_nodes_df_noCONT_proteomics["Sum_of_OGs"]        
og_counts_df = og_counts_df.str.split(',', expand=True)

# Name columns with species name based on species order calculated before 
og_counts_df.columns = species_order_proteomics

normalized_df_proteomics = og_counts_df.copy()

# Convert all values to numeric 
normalized_df_proteomics = normalized_df_proteomics.apply(pd.to_numeric, errors='coerce')

# Perform normalization
for species in normalized_df_proteomics.columns:
    if species in species_OG_sums_noCONT_proteomics:
        normalized_df_proteomics[species] = round((normalized_df_proteomics[species] / species_OG_sums_noCONT_proteomics[species]) * 100, 2)
    else:
        print(f"Warning: {species} not found in species_OG_totals")

# Now normalized_df contains all OGs as percentages of species totals

# Let's create the final table for the iTOL 
merged_df = pd.concat([AllSpecies_nodes_df_noCONT_proteomics['Node'].reset_index(drop=True), normalized_df_proteomics.reset_index(drop=True)], axis=1)
# Add an N before each node name (to merge with the tree node names)
merged_df['Node'] = merged_df['Node'].apply(lambda x: f'N{x}')
```

```python
# Calculate matrix for Barplot iTOL 
# Define the Nodes to be extracted from the tree
Nodes = {
    "Opisthokonta": "N974",  # Root
    "Holozoa": "N976",
    "Filozoa": "N978",
    "Coanozoa": "N981",
    "Metazoa": "N983",
    "Myriazoa": "N984",
    "Parahoxozoa": "N986",
    "Planulozoa": "N989",
    "Bilateria": "N992",
    "Nephrozoa": "N996",
    "Protostomia": "N1002",
    "Phylum": ["N1039", "N1047", "N1061", "N1095", "N1103","N1102"],
    "Species": ["NLEPN", "NEAND", "NHMED", "NCOCO", "NCELE",
                "NLMAR", "NTPIS", "NPACU", "NPTUR", "NSPEC1", "NPLAE",
                "NMISO", "NONUN", "NSMED"]
}

Phylum2Species = { 
    "N1039": ["CELE", "LMAR"],
    "N1047":["PTUR", "SPEC1","TPIS", "PACU"],
    "N1061": ["PLAE","MISO"],
    "N1095": ["ONUN", "SMED"],
    "N1102": ["LEPN"],
    "N1103": ["EAND", "HMED", "COCO"]
    }

REFERENCE_NODE = "N1002"
tree_path = os.path.join(BasePath, "SonicParanoid", "pruned_tree.nwk")

def main(newick_file):
    tree = Phylo.read(newick_file, "newick")
    ref_clade = find_named_clade(tree, REFERENCE_NODE)
    if ref_clade is None:
        raise ValueError(f"Reference node '{REFERENCE_NODE}' not found.")

    records = []

    for phylum_node, species_list in Phylum2Species.items():
        phylum_clade = find_named_clade(tree, phylum_node)
        if phylum_clade is None:
            print(f"Phylum node '{phylum_node}' not found.")
            continue

        for species_code in species_list:
            species_tip = find_named_clade(tree, species_code)
            if species_tip is None:
                print(f"Species tip '{species_code}' not found.")
                continue

            path1 = extract_path_between(tree, ref_clade, phylum_clade)
            path2 = extract_path_between(tree, phylum_clade, species_tip)

            # Remove both ends: start and end node from each path
            if path1:
                path1 = path1[1:-2]
            if path2:
                path2 = path2[1:-2]

            def path_to_names(path):
                return [clade.name if clade.name else f"[internal_{id(clade)}]" for clade in path] if path else None

            records.append({
                "Species": species_code,
                "Phylum_Node": phylum_node,
                "Path_N1002_to_Phylum": path_to_names(path1),
                "Path_Phylum_to_Species": path_to_names(path2)
            })

    df = pd.DataFrame(records)
    return df

if __name__ == "__main__":
    newick_path = tree_path
    df_paths = main(newick_path)
    print(df_paths.head())  # Preview

# Function to calculate sum for each path segment
def calculate_path_sum(path, species, iTOL_df_proteomics):
    
    # Fix species name if needed
    if species == "SPEC1":
        species = "SPEC"
        

    normalized_path = [f"N{node}" if not str(node).startswith("N") else node for node in path]
    total_sum = 0
    for node in normalized_path:
        if node in iTOL_df_proteomics['Node'].values:
            total_sum += iTOL_df_proteomics.loc[iTOL_df_proteomics
                                                ['Node'] == node, species].values[0]
    return total_sum

# Add sum columns for Path_N1002_to_Phylum and Path_Phylum_to_Species
df_paths['Path_N1002_to_Phylum_Sum'] = df_paths.apply(
    lambda row: calculate_path_sum(row['Path_N1002_to_Phylum'], row['Species'], iTOL_df_proteomics), axis=1)

df_paths['Path_Phylum_to_Species_Sum'] = df_paths.apply(
    lambda row: calculate_path_sum(row['Path_Phylum_to_Species'], row['Species'], iTOL_df_proteomics), axis=1)

print(df_paths[['Species', 'Path_N1002_to_Phylum_Sum', 'Path_Phylum_to_Species_Sum']])

#Create the iTOL barplot input 
node_ids = [v for val in Nodes.values() for v in (val if isinstance(val, list) else [val])]

# Filter rows and drop columns
filtered_df = (
    iTOL_df_proteomics[iTOL_df_proteomics['Node'].isin(node_ids)]
    .drop(columns=['Position', 'Radius'])
    .reset_index(drop=True)
)

# Reverse map: node -> group
node_to_group = {}
for group, values in Nodes.items():
    if isinstance(values, list):
        for v in values:
            node_to_group[v] = group
    else:
        node_to_group[values] = group

# Set Node as index to prep for summing
nodes_info_df = filtered_df.set_index("Node")

# Collapse Phylum and Species
collapsed_rows = []
for group in ["Phylum", "Species"]:
    members = Nodes[group]
    existing_members = [m for m in members if m in nodes_info_df.index]
    
    if existing_members:
        summed = nodes_info_df.loc[existing_members].sum()
        summed.name = group
        collapsed_rows.append(summed)

# Safely drop only existing Phylum/Species nodes
all_nodes_to_remove = [m for m in (Nodes["Phylum"] + Nodes["Species"]) if m in nodes_info_df.index]
nodes_info_df = nodes_info_df.drop(all_nodes_to_remove)

# Append collapsed rows
nodes_info_df = pd.concat([nodes_info_df, pd.DataFrame(collapsed_rows)])

# Reorder according to the original Nodes dict
final_order = []
for group in Nodes:
    if group in ["Phylum", "Species"]:
        final_order.append(group)
    else:
        val = Nodes[group]
        if val in nodes_info_df.index:
            final_order.append(val)

# Keep only available rows in final_order
final_order = [n for n in final_order if n in nodes_info_df.index]
nodes_info_df = nodes_info_df.loc[final_order].reset_index().rename(columns={"index": "Node"})

# Assign Group column
nodes_info_df["Group"] = nodes_info_df["Node"].map(node_to_group)
nodes_info_df.loc[nodes_info_df["Node"] == "Phylum", "Group"] = "Phylum"
nodes_info_df.loc[nodes_info_df["Node"] == "Species", "Group"] = "Species"

# Remove the Group column if it exists
if "Group" in nodes_info_df.columns:
    nodes_info_df = nodes_info_df.drop(columns=["Group"])

# Transpose the DataFrame
transposed_df = nodes_info_df.set_index("Node").T

df_paths['Species'] = df_paths['Species'].replace({'SPEC1': 'SPEC'})

# Set index on df_paths
df_paths_indexed = df_paths.set_index('Species')

# Merge on index
merged_df = transposed_df.merge(
    df_paths_indexed[['Path_N1002_to_Phylum_Sum', 'Path_Phylum_to_Species_Sum']],
    left_index=True, right_index=True,
    how='left'
)

# Reorder columns: insert before/after 'Phylum'

cols = [col for col in merged_df.columns if col not in ['Path_N1002_to_Phylum_Sum', 'Path_Phylum_to_Species_Sum']]

# Locate where to insert relative to 'Phylum'
phylum_idx = cols.index('Phylum')

# Insert the two new columns around 'Phylum'
new_order = (
    cols[:phylum_idx] +
    ['Path_N1002_to_Phylum_Sum', 'Phylum', 'Path_Phylum_to_Species_Sum'] +
    cols[phylum_idx+1:]
)

# Reorder columns
final_df = merged_df[new_order]
final_df.to_csv(os.path.join(OutputDir, "BarPlot_input_ALLSPECIES_normalized_100perc_PROTEOMICS.csv"),sep=",", index=True)  # index=True keeps the species names as row index
```

The csv file is the following:

[BarPlot_input_ALLSPECIES_normalized_100perc_PROTEOMICS.csv](Phylostratigraphy%20of%20Hub-containing%20orthogroups%20201f4c09425a8076acfbd5a0a5ea7d64/BarPlot_input_ALLSPECIES_normalized_100perc_PROTEOMICS.csv)

The percentages of this file were incorporated in the iTOL annotation file to be plotted as a barplot: 

[All_Species_iTOL_tree_normalised_barplot_PROTEOMICS.txt](Phylostratigraphy%20of%20Hub-containing%20orthogroups%20201f4c09425a8076acfbd5a0a5ea7d64/All_Species_iTOL_tree_normalised_barplot_PROTEOMICS.txt)

<aside>
💡

To check the proportions of hub-containing orthogroups across experimental conditions we calculated 

</aside>

## Calculating proportions of hub-containing orthogroups per species and per experiment

---

```python
"""
Script created by Klara Eleftheriadi to calculate proportions of hub-containg OGs 
per species and per experiment to be plotted in iTOL as barplot
"""
import pandas as pd
import os
import glob
import numpy as np
from collections import defaultdict

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)

SPECIES = {
    "Annelida": ["Arenicola_marina", "Eisenia_andrei", "Hirudo_medicinalis"],
    "Arthropoda": ["Porcellio_laevis", "Ligia_oceanica"],
    "Mollusca": ["Phorcus_turbinatus", "Physella_acuta", "Theba_pisana", "Siphonaria_pectinata"],
    "Nemertea": ["Leptonemertes_chalicophora", "Tetrastemma_longissimum", "Tetrastemma_melanocephalum"],
    "Nematoda": ["Caenorhabditis_elegans", "Litoditis_marina"],
    "Onychophora": ["Peripatoides_aurorbis"],
    "Platyhelminthes": ["Obama_nungara", "Schmidtea_mediterranea"]
}

SPECIES_CODES = {
    "Arenicola_marina": 'COCO', "Eisenia_andrei": 'EAND', "Hirudo_medicinalis": 'HMED', # Annelida
    "Porcellio_laevis": 'PLAE', "Ligia_oceanica": 'MISO',  # Arthropoda
    "Phorcus_turbinatus": 'PTUR', "Physella_acuta": 'PACU', "Theba_pisana": 'TPIS', "Siphonaria_pectinata": 'SPEC', # Mollusca
    "Leptonemertes_chalicophora": 'LEPN', "Tetrastemma_longissimum": 'TLON', "Tetrastemma_melanocephalum": 'TMEL', # Nemertea
    "Caenorhabditis_elegans": 'CELE', "Litoditis_marina": 'LMAR', # Nematoda
    "Peripatoides_aurorbis": 'PEAU', # Onychophora
    "Obama_nungara": 'ONUN', "Schmidtea_mediterranea": 'SMED'  # Platyhelminthes
}

BasePath = "./"
OutputDir = os.path.join("./")

final_match_df = pd.read_csv(os.path.join(OutputDir, "AllSpecies_Match_Hubs2OGs_per_experiment.txt"), sep="\t")
final_match_df_noCONT = final_match_df[final_match_df["Experiment"]!="CONT"]

def summarize_species_by_experiment(df):
    """
    Summarizes the number of unique OGs per species per experiment.
    Also returns:
      - total unique OGs per species across all experiments (for normalization),
      - a consistent experiment order.
    """

    df = df.copy()

    # Merge TLON and TMEL into TETR 
    df['Species'] = df['Species'].replace({'TLON': 'TETR', 'TMEL': 'TETR'})

    df['OG_list'] = df['OG'].astype(str).apply(lambda x: [og.strip() for og in x.split(',') if og.strip()])
    df_long = df.explode('OG_list')

    # Summarize unique OGs per species and experiment
    summary = (
        df_long.groupby(['Species', 'Experiment'])['OG_list']
        .nunique()
        .reset_index(name='Unique_OGs')
    )

    # Get consistent experiment order ---
    experiment_order = sorted(df['Experiment'].unique())
    species_summary_df = summary.pivot(index='Species', columns='Experiment', values='Unique_OGs').fillna(0).astype(int)

    # Ensure all experiments are present as columns 
    for exp in experiment_order:
        if exp not in species_summary_df.columns:
            species_summary_df[exp] = 0

    species_summary_df = species_summary_df[experiment_order].reset_index()

    # Calculate total unique OGs per species (across all experiments) 
    species_unique_og_counts = (
        df_long.groupby('Species')['OG_list']
        .nunique()
        .to_dict()
    )

    # Normalize counts to percentage by total unique OGs per species 
    normalized_df = species_summary_df.copy()

    for species in normalized_df['Species']:
        total_ogs = species_unique_og_counts.get(species, 1)  # avoid division by 0
        normalized_df.loc[normalized_df['Species'] == species, experiment_order] = (
            normalized_df.loc[normalized_df['Species'] == species, experiment_order] / total_ogs * 100
        )

    normalized_df[experiment_order] = normalized_df[experiment_order].round(2)

    return experiment_order, species_summary_df, species_unique_og_counts, normalized_df

experiment_order, per_Experiment_summary_df, og_totals,normalized_perExperiment_summary_df = summarize_species_by_experiment(final_match_df_noCONT)

normalized_perExperiment_summary_df.columns.name = None  
normalized_perExperiment_summary_df.reset_index(drop=True, inplace=True)  

normalized_perExperiment_summary_df.to_csv(os.path.join(OutputDir,"normalized_og_percentages.csv"),sep=",", index=False)
```

The csv file is the following:

[AllSpecies_AllExperiments_normalized_og_percentages_barplot_input.csv](Phylostratigraphy%20of%20Hub-containing%20orthogroups%20201f4c09425a8076acfbd5a0a5ea7d64/AllSpecies_AllExperiments_normalized_og_percentages_barplot_input.csv)

The iTOL file used for plotting:

[All_Species_iTOL_tree_normalised_PerEXPERIMENT_barplot.txt](Phylostratigraphy%20of%20Hub-containing%20orthogroups%20201f4c09425a8076acfbd5a0a5ea7d64/All_Species_iTOL_tree_normalised_PerEXPERIMENT_barplot.txt)

## Upset plot to check for common OGs

---

Next, we will check for common OGs across species through upset plot 

```r
### UpSet plots for shared OGs by Spp for each species

library(dplyr)
library(UpSetR)
library(ComplexHeatmap)
library(ggplot2)
library(ComplexUpset)

##import full table with "Experiment","Node",	"Species", "OG",	"Proteins" and 	"Phylum" data as a dataframe

table_full <- read.table(file = "AllSpecies_Match_Hubs2OGs_per_experiment_with_Phylum.txt", 
                         header = TRUE, sep = "\t")

## remove "control" lines

table_full<-dplyr::filter(table_full, Experiment!="CONT")

## change "TLON" and "TMEL" to "TETR"

# table_full$Species <- as.character(table_full$Species)

table_full$Species[ table_full$Species == "TLON" ] <- "TETR"

table_full$Species[ table_full$Species == "TMEL" ] <- "TETR"

## TO GET SHARED OGs BETWEEN SPECIES ##

## get OGs for each species

# get species names

spp <- unique(table_full$Species)
spp

## make a vector of uniques OGs for each Spp

for(i in seq_along(spp)) {
  print(spp[i])
  d <- unique(table_full$OG[table_full$Species == spp[i]])
  assign(as.character(spp[i]), d)
  rm(d)
}

# make a list with the vectors of unique OGs

lt <- list(COCO=COCO, EAND=EAND, HMED=HMED, PLAE=PLAE, MISO=MISO, PTUR=PTUR, PACU=PACU, TPIS=TPIS,
           SPEC=SPEC, TETR=TETR, CELE=CELE, LMAR=LMAR, ONUN=ONUN, SMED=SMED, LEPN=LEPN, PEAU=PEAU)

# convert the list into a matrix

lt_mtx <- list_to_matrix(lt)

# matrix to data.frame

lt_mtx_df<-as.data.frame(lt_mtx)

#ordered Vector of spp 

spps<-rev(c("CELE", "LMAR", "PEAU", "PLAE", "MISO", "ONUN",
            "SMED", "PTUR", "SPEC", "TPIS", "PACU", "COCO", 
            "EAND", "HMED", "LEPN", "TETR"))

dev.new()

pdf(paste0("size_20_16colors_sicnames_it.pdf"), width = 13, height = 6, bg = "white") 
#calculate intersection matrix and plot

ComplexUpset::upset(
  lt_mtx_df,
  spps, width_ratio=0.2, height_ratio=1, keep_empty_groups=TRUE, min_size=20,
  mode = "exclusive_intersection", sort_intersections_by=c('degree', 'cardinality'), sort_sets=FALSE,
  themes=upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=12, face="italic", color= "black")))),
  labeller=function(labels) { rev(c("C. elegans", "L. marina", "P. aurorbis", "P. laevis", "L. oceanica",
                                    "O. nungara", "S. mediterranea", "P. turbinatus", "S. pectinata", "T. pisana",
                                    "P. acuta", "A. marina", "E. andrei", "H. medicinalis", "L. chalicophora", "Tetrastemma"))},
  queries=list(
    upset_query(set='CELE', fill='#BB3E03'),
    upset_query(set='LMAR', fill='#005f73'),
    upset_query(set='PEAU', fill='#BB3E03'),
    upset_query(set='PLAE', fill='#BB3E03'),
    upset_query(set='MISO', fill='#005f73'),
    upset_query(set='ONUN', fill='#BB3E03'),
    upset_query(set='SMED', fill="#005f73"),
    upset_query(set='PTUR', fill="#005f73"),
    upset_query(set='SPEC', fill="#005f73"),
    upset_query(set='TPIS', fill='#BB3E03'),
    upset_query(set='PACU', fill="#005f73"),
    upset_query(set='COCO', fill="#005f73"),
    upset_query(set='EAND', fill='#BB3E03'),
    upset_query(set='HMED', fill="#005f73"),
    upset_query(set='LEPN', fill='#BB3E03'),
    upset_query(set='TETR', fill="#005f73")
  ),
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        bar_number_threshold=1,  # show all numbers on top of bars
        width=0.9, # reduce width of the bars
        text=list(
          vjust=-0.1,
          hjust=-0.1,
          angle=45, size = 3.5))
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.3)))
      # + scale_y_continuous(expand=c(0, 1))
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=9),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colorured
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3,
      stroke=0.6)
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.8), position='right')+ geom_text(aes(label=..count..), hjust=-0.2, stat='count', size = 3.5))
  + scale_y_continuous(expand=expansion(mult=c(0, 0.3)))
  + theme(
    axis.line.x=element_line(colour='black'),
    axis.ticks.x=element_line(),
  )
)+ patchwork::plot_layout(heights=c(0.2, 0.35))

dev.off()
```