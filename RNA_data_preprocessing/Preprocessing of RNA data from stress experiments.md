# Preprocessing of RNA data from stress experiments (species-specific workflow)

## Transcriptome assembly

Generate a reference transcriptome assembly using IsoSeq data

``` bash
SPECIES_CODE=CELE

#Step 1 ----->>>>> Primer removal and demultiplexing
lima CCS.${SPECIES_CODE}_iso-${SPECIES_CODE}_iso.bam primers.fasta ${SPECIES_CODE}.ccs.bam --isoseq --peek-guess

#Step 2 ----->>>>> Refine
isoseq3 refine ${SPECIES_CODE}.ccs.NEB_5p--NEB_Clontech_3p.bam primers.fasta ${SPECIES_CODE}.flnc.bam --require-polya

bam2fasta ${SPECIES_CODE}.flnc.bam

# Replace “/” in fasta headers with “_” for compatibility with the tools later 
# and add the SPECIES_CODE

sed '/^>/ {s/\//_/g; s/^>/>{SPECIES_CODE}_iso_/}' {SPECIES_CODE}_flnc.fasta > {SPECIES_CODE}_flnc.mod.fasta

# Clustering with cd-hit-est to obtain the de novo reference transcriptome 

cd-hit-est -i ${SPECIES_CODE}_flnc.mod.fasta -o ${SPECIES_CODE}.flnc.mod.cdhit99.fasta -c $similarity -n 9 -T 25 -M 10000
```

## Preprocessing of short RNA-seq reads

``` bash
BasePath=$(pwd) 

raw_data=$BasePath/../raw_data/ #path/to/fastq/files

#####################################################
######### Quality control in the raw reads ##########
#####################################################

#mkdir $BasePath/quals #Create this directory before running this script
#mkdir $BasePath/trimmed_data #Create this directory before running this script
#mkdir $BasePath/tquals #Create this directory before running this script

fastqc -t 30 -o $BasePath/quals $raw_data/*.fastq.gz

#####################################################
############## Trimming of adapters #################
#####################################################

Trimmo_Out=$BasePath/trimmed_data #path/to/output/of/trimmed/files/after/trimmomatic

cd $raw_data

for file in *1.fastq.gz
do
       root=${file%_1.fastq.gz}
       echo $root
    #Change the path to custom.fa (file with the illumiuna adapters)
       java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 30 \
$raw_data/$root\_1.fastq.gz \
$raw_data/$root\_2.fastq.gz \
$Trimmo_Out/$root\_trimmed_1_paired.fastq.gz \
$Trimmo_Out/$root\_trimmed_1_unpaired.fastq.gz \
$Trimmo_Out/$root\_trimmed_2_paired.fastq.gz \
$Trimmo_Out/$root\_trimmed_2_unpaired.fastq.gz \
ILLUMINACLIP:/path/to/custom.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:15 LEADING:10 TRAILING:10 MINLEN:75 AVGQUAL:30
done && \

#Quality control of trimmed reads
cd $Trimmo_Out
mkdir unpaired
mv *unpaired* unpaired/

fastqc -t 30 -o $BasePath/tquals $Trimmo_Out/*_paired.fastq.gz
```

## Quasi-mapping short RNA-seq reads to reference transcriptome

The quantification files resulted from the following script are used as input for WGCNA

``` bash

SPECIES_NAME=Caenorhabditis_elegans
SPECIES_CODE=CELE

salmon index -t /path/to/${SPECIES_NAME}/${SPECIES_CODE}.flnc.mod.cdhit99.fasta -i ${SPECIES_CODE}_index -k 31 && \

BasePath=/path/to/$SPECIES_NAME/trimmed_data && \

for file in ${BasePath}/*1.fq.gz; 
do
    samp=`basename ${file}`
    root=${samp%_1.fq.gz} 

    echo "Processing sample ${root}"
    salmon quant -i ${SPECIES_CODE}_index -l A \
            -1 ${BasePath}/${root}_1.fq.gz \ 
            -2 ${BasePath}/${root}_2.fq.gz \ 
            -p 40 --validateMappings -o ./quants/${root}_quant
done
```

## Obtain proteome using MATEdb2 pipeline

``` bash
SPECIES_CODE=CELE
BasePath=path/to/cdhit99/fasta/file

#Make TransDecoder directory
mkdir ./TransDecoder

#Execute TransDecoder.LongOrfs

TransDecoder.LongOrfs -t ${BasePath}/${SPECIES_CODE}.flnc.mod.cdhit99.fasta -O ./TransDecoder

ORFS=$( grep -c ">" ./TransDecoder/longest_orfs.pep )
PERC25=$( expr $ORFS / 4 )

#Execute TransDecoder.Predict

TransDecoder.Predict -t ../${SPECIES_CODE}.flnc.mod.cdhit99.fasta -O ./TransDecoder -T $PERC25

#Rename files

#mv ./${SPECIES_CODE}.mod.trinity.fasta.transdecoder.cds ./${SPECIES_CODE}_transdecoder.cds
#mv ./${SPECIES_CODE}.mod.trinity.fasta.transdecoder.pep ./${SPECIES_CODE}_transdecoder.pep
#mv ./${SPECIES_CODE}.mod.trinity.fasta.transdecoder.gff3 ./${SPECIES_CODE}_transdecoder.gff3

#Remove extra files

rm pipeliner*
rm -r ./TransDecoder*
rm ./*.bed

################################################################################
########################         Decontamination        ########################
################################################################################

diamond blastp --query ../2.TransDecoder/${SPECIES_CODE}_transdecoder.pep \
--db $SHARE/databases/DIAMOND_DB/nr.dmnd \
--sensitive --max-target-seqs 1 --evalue 1e-10 --threads 30 --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
--out ./${SPECIES_CODE}.diamond.blastp.out

cds=${SPECIES_CODE}_transdecoder.mod.cds
proteome=${SPECIES_CODE}_transdecoder.mod.pep

blastp_out=${SPECIES_CODE}.diamond.blastp.out

/blobtools2/blobtools create \
    --fasta $cds \
    --hits $blastp_out \
    --hits-cols 1=qseqid,2=staxids,3=bitscore,5=sseqid,10=sstart,11=send,14=evalue \
    --taxrule bestsum \
    --taxdump /mnt/netapp2/Store_csbye/software/taxdump \
    ./BlobDir

#Obtain list of contaminants

python /mnt/netapp2/Store_csbye/scripts/extract_phyla_for_blobtools.py ./BlobDir/bestsum_phylum.json | sed "s/', '/,/g" | tr -d "[]'" > ./contaminants.txt

PHYLA=$(cat ./contaminants.txt)

#Filter cds and protein file
#Filter cds file
/mnt/netapp2/Store_csbye/software/blobtoolkit/blobtools2/blobtools filter \
     --param bestsum_phylum--Keys="$PHYLA" \
     --taxrule bestsum \
     --fasta $cds \
     --summary STDOUT \
     --summary-rank kingdom \
     ./BlobDir > ./${SPECIES_CODE}_blobtools.summary

#Filter pep file
/mnt/netapp2/Store_csbye/software/blobtoolkit/blobtools2/blobtools filter \
     --param bestsum_phylum--Keys="$PHYLA"\
     --taxrule bestsum\
     --fasta $proteome \
     ./BlobDir
```

## Extract longest protein isoforms

``` python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 15:56:22 2023
This script has created to obtain the longest isoforms coming from isoseq assembly 
@author: klara_el
"""

from Bio import SeqIO
import sys

SPECIES_CODE = sys.argv[1]
pep_file = SPECIES_CODE + "_transdecoder.mod.filtered.pep"
long_iso_pep = SPECIES_CODE + "_iso_longiso.pep"

def filter_and_keep_longest(records):
    filtered_records = {}
    for record in records:
        header = record.description.split(".p")[0]
        if header not in filtered_records or len(record.seq) > len(filtered_records[header].seq):
            filtered_records[header] = record

    return filtered_records.values()

records = list(SeqIO.parse(pep_file, "fasta"))
filtered_records = filter_and_keep_longest(records)

with open(long_iso_pep, "w") as output_handle:
    SeqIO.write(filtered_records, output_handle, "fasta")
print(f"Filtered and saved to {long_iso_pep}")
```
