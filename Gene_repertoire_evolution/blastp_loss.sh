#!/bin/bash
#SBATCH --job-name=blastp_gene_loss   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH --time=1-0:0
#SBATCH --output=blastp_gene_loss_array_%A-%a.out    # Standard output and error log
#SBATCH -e blastp_gene_loss_array_%A-%a.err
#SBATCH --array=1-15
#SBATCH --mail-type=ALL
#SBATCH -c 24
#SBATCH --mem=20G

#Paths
WD=$STORE2/Metazoa_OGs/
OUT_PATH=$WD/go_enrichment/gene_loss/blast
SEQ_PATH=$OUT_PATH/FASTA_sequences
SP_FASTA_FILES=/mnt/lustre/scratch/nlsas/home/csic/bye/gim/all_metazoa_plus_Klara/

#Get terrestrial lineage
FILE=$WD/go_enrichment/terr_lineages.txt
TERR_LIN=$(awk -v num_line=$SLURM_ARRAY_TASK_ID '{if(NR==num_line) print $1}' $FILE)

#Get FASTA sequence of OGs lost in a lineage
if [[ ! -f "$SEQ_PATH/${TERR_LIN}_lost.fasta" ]]; then
    module load cesga/2020 seqkit/2.1.0
    while read -r SP; do
        SP_CODE=$(grep -w "$SP" "$WD/species_conversion.txt" | cut -f2)
        seqkit grep -v \
            -f "$OUT_PATH/../${SP_CODE}_subset_genes_${TERR_LIN}_lost.txt" \
            "$SP_FASTA_FILES/$SP" \
            -o "$SEQ_PATH/${SP_CODE}_${TERR_LIN}_lost.fasta"
    done < "$OUT_PATH/../${TERR_LIN}_sps.txt"
    cat "$SEQ_PATH"/*_"${TERR_LIN}_lost.fasta" > "$SEQ_PATH/${TERR_LIN}_lost.fasta"
fi

#Run BLAST against nr
module load blast
blastp -query $SEQ_PATH/${TERR_LIN}_lost.fasta -db $BLASTDB/nr_cluster_seq -evalue 1e-5 -out $OUT_PATH/${TERR_LIN}_lost_blastp.out -num_threads $SLURM_CPUS_PER_TASK -max_target_seqs 5 -outfmt "6 qseqid qlen slen qstart qend sstart send qframe sframe evalue length pident mismatch gaps qcovs stitle staxids sscinames"

