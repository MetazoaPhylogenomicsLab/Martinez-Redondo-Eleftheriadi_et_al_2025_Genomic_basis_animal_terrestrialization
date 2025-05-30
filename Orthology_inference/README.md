## 1. Description

Orthology inference (i.e. identification of gene families) was done using [SonicParanoid2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03298-4).

## 2. Workflow

### Obtaining OGs using SonicParanoid

After trying OrthoFinder (which failed spectacularly in random points after trying for months) and SwifthOrtho (which I could not parallelize), we used SonicParanoid2. SonicParanoid2 includes some machine learning analyses to speed up computation by excluding those pairwise diamond comparisons that are expected not to be informative.

![Overview of SonicParanoid2 as shown in the preprint
a, Graph-based orthology inference pipeline using a novel ML-based approach to substantially reduce execution time of homology searches; b, Domain-based orthology inference pipeline that compares domain architectures using methods from NLP; c, Pairwise orthologous tables predicted using the two pipelines combined to generate ortholog graphs from which the output ortholog groups are inferred.](images/Untitled.png)

Overview of SonicParanoid2 as shown in the preprint
a, Graph-based orthology inference pipeline using a novel ML-based approach to substantially reduce execution time of homology searches; b, Domain-based orthology inference pipeline that compares domain architectures using methods from NLP; c, Pairwise orthologous tables predicted using the two pipelines combined to generate ortholog graphs from which the output ortholog groups are inferred.

The author was contacted to solve some issues about memory.

Execution time was 10 days (first time ran). The following times (some datasets were removed and/or added) it took a bit less as the first step of all-vs-all pairwise alignments were already calculated.

<aside>
⚠️ The input files were obtained from the failed runs of OrthoFinder. Thus, the names of the files and sequences were kept from this. A conversion file was created to keep track of the relationship between the original species codes, the OrthoFinder numbers, and the numbers used by SonicParanoid2 (sorted alphabetically the files)

</aside>

In total, 973 species were included.

```bash
#!/bin/bash

WD=/data/users/martinezredg/SonicParanoid
INPUT_PATH=$WD/../all_metazoa_plus_Klara/
OUT_PATH=$WD/results/

source $WD/sonicparanoid/bin/activate

sonicparanoid -i $INPUT_PATH -o $OUT_PATH -t 48 -p Klara_terrestrialization_2 --project-id Klara_terrestrialization_3 --mode fast --overwrite-tables
```

- Scripts to create conversion files (from OrthoFinder/input files to the species code). The previous one from old SonicParanoid results was used as a base. Note that species from IsoSeq don’t have a number
    
    ```bash
    #Get list of input files for SonicParanoid
    ls $LUSTRE/all_metazoa_plus_Klara > files_used_in_SonicParanoid.txt
    
    #Get the species that were in the old run in the conversion file
    grep -f files_used_in_SonicParanoid.txt old_results_with_PGOT1_no_Klara/species_conversion.txt > species_conversion.txt
    
    #Get files not included in original species_conversion.txt
    cut -f1 old_results_with_PGOT1_no_Klara/species_conversion.txt > tmp.txt
    
    #Add new species to the conversion file
    grep -vf tmp.txt files_used_in_SonicParanoid.txt | awk -F'_' '{print $0 "\t" $1}'>> species_conversion.txt
    rm tmp.txt
    ```
    
    ```bash
    #Now add the species number (in file name) or code (for IsoSeq) to the first column. This will be the information that tell us in the gene names from which species it comes from
    awk '{
        if ($1 ~ /^Species[0-9]+\.fa$/) {
            match($1, /[0-9]+/, num);
            print num[0] "\t" $1 "\t" $2;
        } else {
            print $2 "\t" $1 "\t" $2;
        }
    }' species_conversion.txt > species_num_conversion.txt
    ```
    

Habitat (species_habitat.txt) and phylum (species_phylum.txt and species_phylum_simplified.txt) were manually modified.