#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 50
#SBATCH --mem=1000G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gimartinezredondo@gmail.com
#SBATCH -o sp-OG_full_model_slurm.%j.out
#SBATCH -e sp-OG_full_model_slurm.%j.err
#SBATCH --time=7-0:0

source ml_venv/bin/activate

python execute_full_ML_species_OGs.py
