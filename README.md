# Genomic basis of animal terrestrialisation

## 1. Description
This repository contains the scripts and files needed to reproduce the analyses in: Martínez-Redondo, Eleftheriadi et al. (2025). *Lack of a universal genomic toolkit for animal terrestrialisation*. (under review)

Input proteomes for species for which an IsoSeq reference transcriptome was not generated, and their functional annotation with [FANTASIA](https://github.com/MetazoaPhylogenomicsLab/FANTASIA), can be found in [MATEdb2](https://github.com/MetazoaPhylogenomicsLab/MATEdb2).

## 2. Directories in this repository
- **Species_tree_habitat:** This directory contains the species tree, habitat and phylum information, and ancestral habitat reconstruction files and scripts.
- **Orthology_inference:** This directory contains scripts needed to infere orthogroups using SonicParanoid2, and the output files. Intermediate files are available upon request given their large size to be uploaded in most data repositories.
- **Intermediate_conversion_files:** This directory contains intermediate conversion files used in several analyses.
- **Gene_repertoire_evolution:** This directory contains the scripts and files needed to perform the gene repertoire evolutionary analyses (gene gain and loss).
- **GO_enrichments_terrestrialisation_events:** This directory contains the scripts and files needed to perform the GO enrichment analyses of OGs gained and lost at the independent terrestrialisation events, as well as output files.
- **RNA_data_preprocessing:** This directory contains the scripts needed to perform the IsoSeq data preprocessing and reference transcriptome assembly along with the short reads preprocessing and mapping to the reference transcriptome.
- **Gene_co-expression_networks:** This directory contains the scripts needed to perform the gene coexpression network analyses.
- **Phylostratigraphy_of_Hub_Genes:** This directory contains the scripts and files needed to perform the phylostratigraphic analyses of hub genes.
- **Semantic_similarity_ConstellatoR:** This directory contains the scripts and files needed to perform the functional convergence analyses using [pygosemsim](https://github.com/mojaie/pygosemsim/tree/master/pygosemsim) and [constellatoR](https://github.com/MetazoaPhylogenomicsLab/constellatoR).
- **Pathway_Functional_Convergence_of_Hub_Genes:** This directory contains the scripts needed to perform pathway analyses on hub genes (KO terms).
- **Machine_Learning:** This directory contains the scripts and files needed to perform the machine learning XGBoost classification model for habitat, and posterior analyses on the top OGs obtained.
