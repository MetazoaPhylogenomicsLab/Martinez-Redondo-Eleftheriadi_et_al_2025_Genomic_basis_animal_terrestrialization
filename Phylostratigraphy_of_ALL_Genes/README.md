# Phylostratigraphy using genes from significant modules

#### Step 1: Extraction and aggregation of significant module genes

After extracting significant lists from WGCNA module–trait association outputs across all species and experiments, we summarise this info by experiment, species and phylum. (1.Summarise_all_genes.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 07:54:35 2025 
@author: klara_el

This script scans the WGCNA output folders for each species and retrieves gene
lists from significant modules associated with each experiment.
For every detected gene list, it records the gene identifiers together with the
corresponding experiment, species code, and phylum.
Input files should follow the naming pattern:
    <SPECIES_CODE>_*_significant_module_geneList_for_<EXPERIMENT>.txt

"""
import os
import re
#import pandas as pd

# Define the species mapping
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
    "Arenicola_marina": 'COCO', "Eisenia_andrei": 'EAND', "Hirudo_medicinalis": 'HMED',
    "Porcellio_laevis": 'PLAE', "Ligia_oceanica": 'MISO',
    "Phorcus_turbinatus": 'PTUR', "Physella_acuta": 'PACU', "Theba_pisana": 'TPIS', "Siphonaria_pectinata": 'SPEC',
    "Leptonemertes_chalicophora": 'LEPN', "Tetrastemma_longissimum": 'TLON', "Tetrastemma_melanocephalum": 'TMEL',
    "Caenorhabditis_elegans": 'CELE', "Litoditis_marina": 'LMAR',
    "Peripatoides_aurorbis": 'PEAU',
    "Obama_nungara": 'ONUN', "Schmidtea_mediterranea": 'SMED'
}

BASE_PATH = "/path/to/basepath/"

# Dictionary to hold gene lists at different levels
phylum_gene_dict = {}
all_species_gene_set = {}
species_gene_dict = {}

        
for phylum, species_list in SPECIES.items():
    phylum_gene_dict[phylum] = {}
    
    for species in species_list:
        species_code = SPECIES_CODES[species]
        phylum_gene_dict[phylum][species_code] = {}
        species_dir = os.path.join(BASE_PATH, phylum, species, "WGCNA/3.Relate_modules_to_external_traits/1.signifMods2genes")
        
        if not os.path.exists(species_dir):
            continue  # Skip if the folder doesn't exist
        
        for filename in os.listdir(species_dir):
            match = re.match(rf"{species_code}_.+_significant_module_geneList_for_(.+)\.txt", filename)
            if match:
                experiment = match.group(1)
                
                # Read genes from the file
                with open(os.path.join(species_dir, filename), "r") as f:
                    genes = {line.strip() for line in f if line.strip()}
                
                # Store genes per experiment per species
                if species_code not in species_gene_dict:
                    species_gene_dict[species_code] = {}
                if experiment not in species_gene_dict[species_code]:
                    species_gene_dict[species_code][experiment] = set()
                species_gene_dict[species_code][experiment].update(genes)
                
               # Ensure phylum-level uniqueness with species as a sub-level
                if experiment not in phylum_gene_dict[phylum][species_code]:
                    phylum_gene_dict[phylum][species_code][experiment] = set()
                phylum_gene_dict[phylum][species_code][experiment].update(genes)
                
                # Aggregate at all-species level (not used for now need to be changed if to be used)
                if experiment not in all_species_gene_set:
                    all_species_gene_set[experiment] = set()
                all_species_gene_set[experiment].update(genes)
                
# Make species-level genes unique across experiments
final_species_gene_dict = {}
for species_code, exp_dict in species_gene_dict.items():
    final_species_gene_dict[species_code] = {}
    for experiment, genes in exp_dict.items():
        final_species_gene_dict[species_code][experiment] = set(genes)
        species_output_file = os.path.join(BASE_PATH,"SonicParanoid/Merging_with_all_genes", f"{species_code}_{experiment}_genes.txt")
        with open(species_output_file, "w") as out:
            out.write("Gene\tExperiment\n")  # Write header
            for gene in sorted(genes):
                out.write(f"{gene}\t{experiment}\n")

# Write phylum-level files with species information in 3rd column
for phylum, species_dict in phylum_gene_dict.items():
    for species_code, exp_dict in species_dict.items():
        for experiment, genes in exp_dict.items():
            phylum_output_file = f"{phylum}_{experiment}_genes.txt"
            with open(phylum_output_file, "w") as out:
                out.write("Gene\tExperiment\tSpecies\n")  # Write header
                for gene in sorted(genes):
                    out.write(f"{gene}\t{experiment}\t{species_code}\n")

# Merge per phylum per exepriment
for phylum, species_dict in phylum_gene_dict.items():
    phylum_exp_dict = {}  # Collect all genes per experiment within the phylum
    
    for species_code, exp_dict in species_dict.items():
        for experiment, genes in exp_dict.items():
            if experiment not in phylum_exp_dict:
                phylum_exp_dict[experiment] = []  # Store all (gene, species) pairs
            
            for gene in genes:
                phylum_exp_dict[experiment].append((gene, species_code))
    
    # Now write each experiment file with all species' genes
    for experiment, gene_species_list in phylum_exp_dict.items():
        phylum_output_file = f"{phylum}_{experiment}_Allgenes.txt"
        with open(phylum_output_file, "w") as out:
            out.write("Gene\tExperiment\tSpecies\n")  # Write header
            for gene, species_code in sorted(gene_species_list):
                out.write(f"{gene}\t{experiment}\t{species_code}\n")

# Merge per phylum regardless of experiment 
# Write a single file per phylum with all genes from all species and experiments
for phylum, species_dict in phylum_gene_dict.items():
    phylum_gene_list = []  # Collect all (gene, experiment, species) entries
    
    for species_code, exp_dict in species_dict.items():
        for experiment, genes in exp_dict.items():
            for gene in genes:
                phylum_gene_list.append((gene, experiment, species_code))
    
    # Now write the single merged phylum file
    phylum_output_file = f"{phylum}_Allgenes.txt"
    with open(phylum_output_file, "w") as out:
        out.write("Gene\tExperiment\tSpecies\n")  # Write header
        for gene, experiment, species_code in sorted(phylum_gene_list):
            out.write(f"{gene}\t{experiment}\t{species_code}\n")

# Write the all-species file
#for experiment, genes in all_species_gene_set.items():
#    all_species_output_file = f"All_species_{experiment}_genes.txt"
#    with open(all_species_output_file, "w") as out:
#        out.write("\n".join(sorted(genes)))
```

#### **Step 2: Matching significant genes to translated protein identifiers**

The phylum-level gene tables generated above (`<phylum>_Allgenes.txt`) were used to retrieve the corresponding translated protein identifiers for each gene. For each species, gene names were matched against the headers of the longest isoforms proteome files. The matched protein IDs were added as a new column, producing one table per phylum. (2.MatchGeneNames2Proteins.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 10:09:32 2025
@author: klara_el
"""

import os
import glob
import pandas as pd

SPECIES = {
    "Annelida": ["Arenicola_marina", "Eisenia_andrei", "Hirudo_medicinalis"],
    "Arthropoda": ["Porcellio_laevis", "Ligia_oceanica"],
    "Mollusca": ["Phorcus_turbinatus", "Physella_acuta", "Theba_pisana", "Siphonaria_pectinata"],
    "Nemertea": ["Leptonemertes_chalicophora", "Tetrastemma_longissimum", "Tetrastemma_melanocephalum"],
    "Nematoda": ["Caenorhabditis_elegans", "Litoditis_marina"],
    "Onychophora": ["Peripatoides_aurorbis"],
    "Platyhelminthes": ["Obama_nungara", "Schmidtea_mediterranea"]
}

def load_proteome(proteome_file):
    """Load proteome headers from a .pep file."""
    protein_dict = {}
    with open(proteome_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                protein_id = line[1:].strip()  # Remove '>'
                gene_name = protein_id.split('.')[0]  # Extract gene name
                if gene_name not in protein_dict:
                    protein_dict[gene_name] = []
                protein_dict[gene_name].append(protein_id)
    return protein_dict

def process_phylum_file(phylum_file, proteome_dir, output_dir):
    # Process a single phylum file, checking for gene matches in the proteomes.
    df = pd.read_csv(phylum_file, sep='\t', header=0, names=['Gene', 'Experiment', 'Species'])
    phylum_name = os.path.basename(phylum_file).replace('_Allgenes.txt', '')
    #print(phylum_name)
    # Dictionary to store proteome data
    proteome_cache = {}
    
    matched_proteins = []
    for _, row in df.iterrows():
        gene, species = row['Gene'], row['Species']
        #print(species)
        proteome_file = os.path.join(proteome_dir, f"{species}_transdecoder.filtered.longiso.mod.pep")
        #print(proteome_file)
        # Load proteome only if not already cached
        if species not in proteome_cache:
            if os.path.exists(proteome_file):
                proteome_cache[species] = load_proteome(proteome_file)
            else:
                proteome_cache[species] = {}  # Empty if file doesn't exist
        
        # Check for matching proteins
        matched = proteome_cache[species].get(gene, ['NA'])
        matched_proteins.append(",".join(matched))  # Store all matching proteins
    
    # Add Protein column and save output
    df['Protein'] = matched_proteins
    output_file = os.path.join(output_dir, f"{phylum_name}_AllGenes_with_Proteins.txt")
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed {phylum_file}, saved results to {output_file}")

# Directories
phyla_dir = "/SonicParanoid/Merging_with_all_genes/"
proteome_dir = "/longest_isoforms/"
output_dir = "/SonicParanoid/Merging_with_all_genes/new"  # Save results in the same directory

# Process all phylum files
for phylum in SPECIES.keys():
    print(phylum)
    phylum_file = os.path.join(phyla_dir, f"{phylum}_Allgenes.txt")
    print(phylum_file)
    if os.path.exists(phylum_file):
        process_phylum_file(phylum_file, proteome_dir, output_dir)

```

#### Step 3: Assigning significant genes to orthogroups

The following script uses the protein-annotated phylum-level gene files generated in the previous step (<phylum>_AllGenes_with_Proteins.txt) and matches the translated protein identifiers to orthogroups from the SonicParanoid output table (3.AllModsGenes2OGs_perPhylum.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 08:54:12 2025

@author: klara_el
"""

import pandas as pd
import os

# Define paths
#BasePath = "/home/klara/drive/PhD/Gene_CoExpression_Networks/Invertebrates"
BasePath = "/path/to/basepath"
SONIC_PARANOID_PATH = BasePath + "/SonicParanoid/"
OGs_FILE = os.path.join(SONIC_PARANOID_PATH, "SonicParanoid_invertebrates_experiments.mod.txt")

# Read the OGs data
OGs_df = pd.read_csv(OGs_FILE, sep="\t")
protein_data_filepath = os.path.join(BasePath, "SonicParanoid/Merging_with_all_genes/new/") #Files resulted from script 2 

# Phyla with respective Species names
SPECIES = {
    "Annelida": ["Arenicola_marina", "Eisenia_andrei", "Hirudo_medicinalis"],
    "Arthropoda": ["Porcellio_laevis", "Ligia_oceanica"],
    "Mollusca": ["Phorcus_turbinatus", "Physella_acuta", "Theba_pisana", "Siphonaria_pectinata"],
    "Nemertea": ["Leptonemertes_chalicophora", "Tetrastemma_longissimum", "Tetrastemma_melanocephalum"],
    "Nematoda": ["Caenorhabditis_elegans", "Litoditis_marina"],
    "Onychophora": ["Peripatoides_aurorbis"],
    "Platyhelminthes": ["Obama_nungara", "Schmidtea_mediterranea"]
}

# Species and their respective codes
SPECIES_CODES = {
    "Arenicola_marina": 'COCO', "Eisenia_andrei": 'EAND', "Hirudo_medicinalis": 'HMED', # Annelida
    "Porcellio_laevis": 'PLAE', "Ligia_oceanica": 'MISO',  # Arthropoda
    "Phorcus_turbinatus": 'PTUR', "Physella_acuta": 'PACU', "Theba_pisana": 'TPIS', "Siphonaria_pectinata": 'SPEC', # Mollusca
    "Leptonemertes_chalicophora": 'LEPN', "Tetrastemma_longissimum": 'TLON', "Tetrastemma_melanocephalum": 'TMEL', # Nemertea
    "Caenorhabditis_elegans": 'CELE', "Litoditis_marina": 'LMAR', # Nematoda
    "Peripatoides_aurorbis": 'PEAU', # Onychophora
    "Obama_nungara": 'ONUN', "Schmidtea_mediterranea": 'SMED'  # Platyhelminthes
}
results = []
    
def match_Genes_with_OGs(species_proteins,species_code, OGs_df):
# Keep species-specific OGs
    OGs_df_species = OGs_df.iloc[:, :4].join(OGs_df[species_code]) 
    
    for _, row in OGs_df_species.iterrows():
        group_id = row['group_id']
        node = row['Node']
        OG_proteins = [protein for protein in row[species_code].split(',') if '*' not in protein] if pd.notna(row[species_code]) else []
        
        common_proteins = list(set(OG_proteins) & set(species_proteins))
        if common_proteins:
            results.append({
                'Node': node,
                'OG': group_id,
                'Species': species_code,
                'Proteins': ', '.join(common_proteins)
            })
            
    match_Genes2OGs_df = pd.DataFrame(results, columns=['Node', 'OG', 'Species', 'Proteins']) # A "flat" df without merging any info

    return match_Genes2OGs_df

species_proteins_dict = {}  # Dictionary to store species-specific proteins
for phylum, species_list in SPECIES.items():
    phylum_file = os.path.join(protein_data_filepath, f"{phylum}_AllGenes_with_Proteins.txt")
    
    if not os.path.exists(phylum_file):
        print(f"File not found: {phylum_file}")
        continue
    
    # Read per-phylum data
    phylum_df = pd.read_csv(phylum_file, sep="\t")

    # Drop NA values in Protein column
    phylum_df = phylum_df.dropna(subset=["Protein"])

    # Group proteins by species
    for species in species_list:
        print(species)
        species_code = SPECIES_CODES.get(species, None)
        if species_code:
            species_proteins = phylum_df.loc[phylum_df["Species"] == species_code, "Protein"].tolist()
            species_proteins_dict[species_code] = species_proteins

            Genes2OGs_df = match_Genes_with_OGs(species_proteins,species_code,OGs_df)

output_file = os.path.join(protein_data_filepath, "All_Species_GenesAndProteins_to_OGs_mapping_new.txt")

Genes2OGs_df.to_csv(output_file, sep="\t", index=False)

OGs = Genes2OGs_df["OG"].drop_duplicates().reset_index()

# Save unique OGs in one-column list to check with top130
og_output_file = os.path.join(protein_data_filepath, "OGs_unique_list_from_ALL_genes.txt") 
Genes2OGs_df["OG"].drop_duplicates().to_csv(og_output_file, sep="\t", index=False, header=False)
print(f"OG list saved to: {og_output_file}")
```

#### Step 4: Phylostratigraphic profiling and iTOL barplot preparation

Now we need to perform module-associated genes phylostratigraphy and generate iTOL barplot input files. This script uses all significant module-associated genes from the WGCNA networks, mapped to translated proteins and SonicParanoid orthogroups, to summarize the phylogenetic origin of the gene sets across species and experiments. The output tables can be used in the iTOL formatted txt files for visualisation (Phylostratigraphy_Barplots_ALLGenes.py).

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: klara_el
"""
import os
import glob
from collections import defaultdict

import pandas as pd
from Bio import Phylo

BASE = "/path/to/basepath"

SONIC_PARANOID_DIR = os.path.join(BASE, "SonicParanoid")
OGS_FILE = os.path.join(SONIC_PARANOID_DIR, "SonicParanoid_invertebrates_experiments.mod.txt") 
TREE_FILE = os.path.join(SONIC_PARANOID_DIR, "pruned_tree.nwk") # With the 17 species and the internal nodes

# INPUT a folder containing multiple *.tsv/*.txt/*.csv files 
# with the same 4 columns 
# INCLUDING ALL GENES FROM SIGNIFICANT MODULES PER PHYLUM AND EXPERIMENT
ALL_GENES_DIR = os.path.join(SONIC_PARANOID_DIR, "Merging_with_all_genes/Phylostratigraphy_ALLGENES")

OUTDIR = os.path.join(ALL_GENES_DIR, "AllGenes_Phylo")
os.makedirs(OUTDIR, exist_ok=True)

REMOVE_CONTROL = True
CONTROL_NAME = "CONT"

# TMEL + TLON merged to TETR after mapping
TETR_MEMBERS = ["TLON", "TMEL"]  # keep this order irrelevant

# Fix for consistency 
SPECIES_RENAME = {"SPEC1": "SPEC"} 

# Reference node for Protostomia path calculations
REFERENCE_NODE = "N1002"

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
    "N1047": ["PTUR", "SPEC1", "TPIS", "PACU"],
    "N1060": ["PEAU"],
    "N1061": ["PLAE", "MISO"],
    "N1095": ["ONUN", "SMED"],
    "N1102": ["LEPN", "TMEL", "TLON"],
    "N1103": ["EAND", "HMED", "COCO"]
}

def load_allgenes_longtable(folder: str) -> pd.DataFrame:
    dfs = []
    
    if folder and os.path.isdir(folder):
        for fp in glob.glob(os.path.join(folder, "*")):
            if fp.endswith((".tsv", ".txt", ".csv")):
                dfs.append(pd.read_csv(fp, sep=None, engine="python"))

    if not dfs:
        raise FileNotFoundError(
            "No ALL-GENES input found. Provide ALL_GENES_FILE or ALL_GENES_DIR with readable files."
        )

    df = pd.concat(dfs, ignore_index=True)

    # Robust column renaming based on required semantics
    # Accepts any capitalization
    lower_map = {c.lower(): c for c in df.columns}
    required = ["gene", "experiment", "species", "protein"]
    missing = [r for r in required if r not in lower_map]
    if missing:
        raise ValueError(f"ALL-GENES table missing columns: {missing}. Found columns: {list(df.columns)}")

    df = df.rename(columns={
        lower_map["gene"]: "Gene",
        lower_map["experiment"]: "Experiment",
        lower_map["species"]: "Species",
        lower_map["protein"]: "Protein"
    })

    # Basic cleaning
    df = df.dropna(subset=["Protein", "Species", "Experiment"])
    df["Protein"] = df["Protein"].astype(str).str.strip()
    df["Species"] = df["Species"].astype(str).str.strip()
    df["Experiment"] = df["Experiment"].astype(str).str.strip()

    # Remove true duplicates within the same species+experiment (same protein repeated)
    df = df.drop_duplicates(subset=["Protein", "Species", "Experiment"])

    if REMOVE_CONTROL:
        df = df[df["Experiment"] != CONTROL_NAME].copy()

    return df

def build_species_experiment_dict(df_long: pd.DataFrame):
    """
    Returns dict: species -> experiment -> list(proteins)
    """
    d = defaultdict(lambda: defaultdict(set))
    for p, s, e in zip(df_long["Protein"], df_long["Species"], df_long["Experiment"]):
        d[s][e].add(p)
    return {s: {e: list(ps) for e, ps in exps.items()} for s, exps in d.items()}

# MAP PROTEINS TO OGs PER SPECIES PER EXPERIMENT
def map_proteins_to_OGs_for_species_experiment(proteins, OGs_df, species_code, experiment_id):
    """
    Returns:
      match_df: one row per OG hit (keeps Experiment)
      summary_df_node_species: collapsed across experiments later, but here we already summarize to Node×Species
                              by collapsing duplicates *within the current experiment mapping*.
    """
    if species_code not in OGs_df.columns:
        raise KeyError(
            f"Species column '{species_code}' not found in OG matrix. "
            f"Available species cols: {[c for c in OGs_df.columns if c not in ['group_id','Node','group_size','sp_in_grp','seed_ortholog_cnt']]}"
        )

    protset = set(proteins)
    rows = []

    for _, row in OGs_df.iterrows():
        if pd.isna(row[species_code]):
            continue

        group_id = row["group_id"]
        node = row["Node"]

        cell = str(row[species_code])
        og_proteins = [p.strip() for p in cell.split(",") if p.strip() and "*" not in p]

        common = sorted(protset.intersection(og_proteins))
        if common:
            rows.append({
                "Experiment": experiment_id,
                "Node": node,
                "Species": species_code,
                "OG": str(group_id),
                "Proteins": ",".join(common)
            })

    match_df = pd.DataFrame(rows, columns=["Experiment", "Node", "Species", "OG", "Proteins"])

    # Summarize within (Node, Species, Experiment) – still keeps experiment in string if needed
    # But for your main Figure3B-like barplot we will later collapse across experiments anyway.
    if match_df.empty:
        summary_df = pd.DataFrame(columns=["Node", "Species", "Experiment", "OG", "Proteins"])
        return match_df, summary_df

    summary_df = match_df.groupby(["Node", "Species"], as_index=False).agg({
        "Experiment": lambda x: ",".join(sorted(set(x))),
        "OG": lambda x: ",".join(sorted(set(map(str, x)))),
        "Proteins": lambda x: ",".join(sorted(set(map(str, x))))
    })

    return match_df, summary_df

def collapse_across_experiments_node_species(all_summary_df: pd.DataFrame) -> pd.DataFrame:
    """
    Takes concatenated summaries (Node×Species per experiment) and collapses across experiments:
    - OGs: unique OGs per Node×Species across ALL experiments
    - Proteins: unique proteins per Node×Species across ALL experiments (informative only)
    """
    if all_summary_df.empty:
        return pd.DataFrame(columns=["Node", "Species", "OG", "Proteins", "Sum_of_OGs", "Sum_of_proteins"])

    df = all_summary_df.copy()

    def uniq_join_csv(series):
        items = set()
        for v in series.dropna().astype(str):
            for it in v.split(","):
                it = it.strip()
                if it:
                    items.add(it)
        return ",".join(sorted(items))

    collapsed = df.groupby(["Node", "Species"], as_index=False).agg({
        "OG": uniq_join_csv,
        "Proteins": uniq_join_csv,
    })

    collapsed["Sum_of_OGs"] = collapsed["OG"].apply(lambda x: len([v for v in x.split(",") if v.strip()]))
    collapsed["Sum_of_proteins"] = collapsed["Proteins"].apply(lambda x: len([v for v in x.split(",") if v.strip()]))

    return collapsed

def build_node_species_matrix(collapsed_node_species_df: pd.DataFrame):

    counts = collapsed_node_species_df.pivot_table(
        index="Node",
        columns="Species",
        values="Sum_of_OGs",
        aggfunc="sum",
        fill_value=0
    ).astype(float)

    # TETR merge 
    if ("TLON" in counts.columns) and ("TMEL" in counts.columns):
        counts["TETR"] = counts["TLON"] + counts["TMEL"]
        counts = counts.drop(columns=["TLON", "TMEL"])

    # normalization denominator = total across nodes per species
    totals = counts.sum(axis=0).replace(0, 1.0)
    normalized = (counts / totals) * 100.0
    normalized = normalized.round(2)

    return counts, normalized

# TREE PATH SUMS (Protostomia->Phylum and Phylum->Species)
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

def build_paths_df(tree_file: str) -> pd.DataFrame:
    tree = Phylo.read(tree_file, "newick")
    ref_clade = find_named_clade(tree, REFERENCE_NODE)
    if ref_clade is None:
        raise ValueError(f"Reference node '{REFERENCE_NODE}' not found in tree.")

    records = []

    for phylum_node, species_list in Phylum2Species.items():
        phylum_clade = find_named_clade(tree, phylum_node)
        if phylum_clade is None:
            print(f"WARNING: Phylum node '{phylum_node}' not found in tree.")
            continue

        for species_code in species_list:
            species_tip = find_named_clade(tree, species_code)
            if species_tip is None:
                print(f"WARNING: Species tip '{species_code}' not found in tree.")
                continue

            path1 = extract_path_between(tree, ref_clade, phylum_clade)
            path2 = extract_path_between(tree, phylum_clade, species_tip)

            # mimic your trimming
            if path1:
                path1 = path1[1:-2]
            if path2:
                path2 = path2[1:-2]

            def path_to_names(path):
                return [c.name for c in path] if path else []

            records.append({
                "Species": species_code,
                "Phylum_Node": phylum_node,
                "Path_N1002_to_Phylum": path_to_names(path1),
                "Path_Phylum_to_Species": path_to_names(path2)
            })

    df = pd.DataFrame(records)

    # Merge TLON/TMEL into TETR 
    tetr_rows = df[df["Species"].isin(TETR_MEMBERS)]
    if not tetr_rows.empty:
        merged_row = tetr_rows.iloc[0].copy()
        merged_row["Species"] = "TETR"
        df = df[~df["Species"].isin(TETR_MEMBERS)]
        df = pd.concat([df, pd.DataFrame([merged_row])], ignore_index=True)

    # SPEC1 rename only for matching final iTOL columns
    df["Species"] = df["Species"].replace(SPECIES_RENAME)

    return df

def calculate_path_sum(path, species, iTOL_node_df):
    # SPEC1->SPEC
    species = SPECIES_RENAME.get(species, species)

    total_sum = 0.0
    for node in path:
        node = str(node)
        node = node if node.startswith("N") else f"N{node}"
        if node in set(iTOL_node_df["Node"]):
            total_sum += float(iTOL_node_df.loc[iTOL_node_df["Node"] == node, species].values[0])
    return total_sum

# iTOL BARPLOT MATRIX (collapse Phylum + Species, transpose, sum paths from N1002)
def build_itol_barplot_matrix(iTOL_node_df: pd.DataFrame, df_paths: pd.DataFrame) -> pd.DataFrame:
    # flatten node IDs from Nodes dict
    node_ids = [v for val in Nodes.values() for v in (val if isinstance(val, list) else [val])]

    filtered = iTOL_node_df[iTOL_node_df["Node"].isin(node_ids)].copy().reset_index(drop=True)

    # Drop Position/Radius if present (some of your intermediate files had them)
    for col in ["Position", "Radius"]:
        if col in filtered.columns:
            filtered = filtered.drop(columns=[col])

    nodes_info_df = filtered.set_index("Node")

    # Collapse Phylum and Species
    collapsed_rows = []
    for group in ["Phylum", "Species"]:
        members = Nodes[group]
        existing = [m for m in members if m in nodes_info_df.index]
        if existing:
            summed = nodes_info_df.loc[existing].sum(numeric_only=True)
            summed.name = group
            collapsed_rows.append(summed)

    # remove original members
    to_remove = [m for m in (Nodes["Phylum"] + Nodes["Species"]) if m in nodes_info_df.index]
    nodes_info_df = nodes_info_df.drop(to_remove, errors="ignore")

    # append collapsed
    if collapsed_rows:
        nodes_info_df = pd.concat([nodes_info_df, pd.DataFrame(collapsed_rows)], axis=0)

    # reorder rows by Nodes dict
    final_order = []
    for group, val in Nodes.items():
        if group in ["Phylum", "Species"]:
            final_order.append(group)
        else:
            if val in nodes_info_df.index:
                final_order.append(val)
    final_order = [x for x in final_order if x in nodes_info_df.index]

    nodes_info_df = nodes_info_df.loc[final_order].reset_index().rename(columns={"index": "Node"})

    # transpose => species rows, node columns
    transposed = nodes_info_df.set_index("Node").T

    # set index on df_paths
    df_paths_indexed = df_paths.set_index("Species")

    # merge path sums columns
    merged = transposed.merge(
        df_paths_indexed[["Path_N1002_to_Phylum_Sum", "Path_Phylum_to_Species_Sum"]],
        left_index=True, right_index=True, how="left"
    )

    # insert sums around 'Phylum' column
    cols = [c for c in merged.columns if c not in ["Path_N1002_to_Phylum_Sum", "Path_Phylum_to_Species_Sum"]]
    phylum_idx = cols.index("Phylum") if "Phylum" in cols else len(cols) // 2

    new_order = (
        cols[:phylum_idx] +
        ["Path_N1002_to_Phylum_Sum", "Phylum", "Path_Phylum_to_Species_Sum"] +
        cols[phylum_idx + 1:]
    )

    new_order = [c for c in new_order if c in merged.columns]
    final_df = merged[new_order]

    return final_df

# Per-species per-experiment normalized OG %
def normalize_per_experiment(match_df: pd.DataFrame) -> pd.DataFrame:
    
    df = match_df.copy()
    if df.empty:
        return pd.DataFrame()

    # explode OG list
    df["OG_list"] = df["OG"].astype(str).apply(lambda x: [og.strip() for og in x.split(",") if og.strip()])
    df_long = df.explode("OG_list")

    # unique OGs per species per experiment
    summary = (
        df_long.groupby(["Species", "Experiment"])["OG_list"]
        .nunique()
        .reset_index(name="Unique_OGs")
    )

    # total unique OGs per species across all experiments
    denom = (
        df_long.groupby("Species")["OG_list"]
        .nunique()
        .to_dict()
    )

    exp_order = sorted(df["Experiment"].unique())
    wide = summary.pivot(index="Species", columns="Experiment", values="Unique_OGs").fillna(0).astype(float)

    # merge TLON/TMEL -> TETR in wide
    if ("TLON" in wide.index) and ("TMEL" in wide.index):
        tetr_row = wide.loc["TLON"] + wide.loc["TMEL"]
        wide = wide.drop(index=["TLON", "TMEL"])
        wide.loc["TETR"] = tetr_row

        # merge denom too
        denom["TETR"] = denom.get("TLON", 0) + denom.get("TMEL", 0)
        denom.pop("TLON", None)
        denom.pop("TMEL", None)

    # normalize
    norm = wide.copy()
    for sp in norm.index:
        d = denom.get(sp, 1)
        if d == 0:
            d = 1
        norm.loc[sp, exp_order] = (norm.loc[sp, exp_order] / d) * 100.0

    norm = norm.round(2).reset_index()
    norm.columns.name = None
    return norm

OGs_df = pd.read_csv(OGS_FILE, sep="\t")
allgenes_df = load_allgenes_longtable(ALL_GENES_DIR) #DF reading all 17 species files together
print(f"ALL genes rows after cleaning (and CONT removal={REMOVE_CONTROL}): {allgenes_df.shape[0]}")

# Build a nested dict Species:exp1{genes}, exp2{genes} etc
species_experiment_dict = build_species_experiment_dict(allgenes_df) 

all_match_dfs = []
all_summary_dfs = []

for species_code, experiments in species_experiment_dict.items():
    for experiment_id, proteins in experiments.items():
        print(f"  Mapping Species={species_code} Experiment={experiment_id} Proteins={len(proteins)}") # Unique proteins per experiment 
        match_df, summary_df = map_proteins_to_OGs_for_species_experiment(
            proteins, OGs_df, species_code, experiment_id
        )
        if not match_df.empty:
            all_match_dfs.append(match_df)
        if not summary_df.empty:
            all_summary_dfs.append(summary_df)

if not all_match_dfs:
    raise ValueError("No matches found between ALL genes and OG matrix. Check Protein IDs and Species codes.")

match_df = pd.concat(all_match_dfs, ignore_index=True)
summary_df_per_exp_collapsed = pd.concat(all_summary_dfs, ignore_index=True) if all_summary_dfs else pd.DataFrame()

# Write match table
match_out = os.path.join(OUTDIR, "AllGenes_Match_Genes2OGs_per_experiment.txt")
match_df.to_csv(match_out, sep="\t", index=False)
print("Wrote:", match_out)

# Collapse across experiments for Node×Species (unique OGs across all experiments)
collapsed_node_species_df = collapse_across_experiments_node_species(summary_df_per_exp_collapsed)

summary_out = os.path.join(OUTDIR, "AllGenes_Summary_NodeSpecies_AllExperiments.txt")
collapsed_node_species_df.to_csv(summary_out, sep="\t", index=False)
print("Wrote:", summary_out)

counts_df, normalized_df = build_node_species_matrix(collapsed_node_species_df)

# Build iTOL node table: Node must have N prefix to match tree internal nodes
iTOL_node_df = normalized_df.copy()
iTOL_node_df.index = iTOL_node_df.index.map(lambda x: f"N{x}" if not str(x).startswith("N") else str(x))
iTOL_node_df = iTOL_node_df.reset_index().rename(columns={"index": "Node"})

# Apply SPEC1->SPEC rename to columns if needed
iTOL_node_df = iTOL_node_df.rename(columns=SPECIES_RENAME)

print("Building paths and computing path sums")
df_paths = build_paths_df(TREE_FILE)

# Compute path sums using iTOL_node_df (normalized % at each node)
df_paths["Path_N1002_to_Phylum_Sum"] = df_paths.apply(
    lambda r: calculate_path_sum(r["Path_N1002_to_Phylum"], r["Species"], iTOL_node_df),
    axis=1
)
df_paths["Path_Phylum_to_Species_Sum"] = df_paths.apply(
    lambda r: calculate_path_sum(r["Path_Phylum_to_Species"], r["Species"], iTOL_node_df),
    axis=1
)

final_barplot_df = build_itol_barplot_matrix(iTOL_node_df, df_paths)

barplot_out = os.path.join(OUTDIR, "BarPlot_input_ALLSPECIES_normalized_100perc_ALLGENES.csv")
final_barplot_df.to_csv(barplot_out, sep=",", index=True)
print("Wrote:", barplot_out)

# Per-experiment table
print("Per-species per-experiment normalized OG%")
norm_exp_df = normalize_per_experiment(match_df)
if not norm_exp_df.empty:
    exp_out = os.path.join(OUTDIR, "AllSpecies_AllExperiments_normalized_og_percentages_ALLGENES.csv")
    norm_exp_df.to_csv(exp_out, sep=",", index=False)
    print("Wrote:", exp_out)

print("DONE")
```