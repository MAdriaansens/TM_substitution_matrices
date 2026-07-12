#!/bin/bash
#SBATCH --job-name=AMOA_MFQmaker
#SBATCH --time=32:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=20GB          # Memory in MB
#SBATCH --cpus-per-task=15
#SBATCH --output=slurm_outputQQ/AAQmaker_output%A.out
#SBATCH --error=slurm_outputQQ/AAQmaker_error%A.err

HMMalign=/nesi/nobackup/uc04105/PDB_alpha/Matrice_generation/HMMalign/pMO
module load IQ-TREE/3.1.1-foss-2023a
scp ${HMMalign}/Pfam12942_vs_curated_AMOA_filter_curated_merged_final_linsialigned.fasta .
iqtree3 -nt AUTO -s Pfam12942_vs_curated_AMOA_filter_curated_merged_final_linsialigned.fasta  -m MFP
