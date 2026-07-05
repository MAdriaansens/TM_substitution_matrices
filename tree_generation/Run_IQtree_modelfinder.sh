#!/bin/bash
#SBATCH --job-name=AMOA_Qmaker
#SBATCH --time=2:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=4GB          # Memory in MB
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm_outputQQ/Qmaker_output%A.out
#SBATCH --error=slurm_outputQQ/Qmaker_error%A.err


module load IQ-TREE/3.1.1-foss-2023a

iqtree3 -nt AUTO -s PF12942_AMOA_filter_curated_manual.fasta -m MFP


