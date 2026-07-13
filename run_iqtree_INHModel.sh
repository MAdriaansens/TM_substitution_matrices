#!/bin/bash
#SBATCH --job-name=AMOA_Qmaker
#SBATCH --time=4:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=20GB          # Memory in MB
#SBATCH --cpus-per-task=10
#SBATCH --output=slurm_outputQQ/Qmaker_output%A.out
#SBATCH --error=slurm_outputQQ/Qmaker_error%A.err

Matrice_path=/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/Matrice_Qmaker/fl/iqtreeQmaker_fl/Qmatrice/Q.AMOA_FL_Matrice_GTR20_FO

module load IQ-TREE/3.1.1-foss-2023a
iqtree3 --seed 1 -T AUTO -s Pfam12942_vs_curated_AMOA_filter_curated_merged_final_linsialigned_Fl.fasta -m ${Matrice_path} -B 1000
