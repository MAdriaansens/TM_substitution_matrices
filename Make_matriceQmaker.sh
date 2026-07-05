#!/bin/bash
#SBATCH --job-name=Qmaker
#SBATCH --time=2:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=4GB          # Memory in MB
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_outputQQ/Qmaker_output%A.out
#SBATCH --error=slurm_outputQQ/Qmaker_error%A.err


module load IQ-TREE/3.1.1-foss-2023a

iqtree3 --seed 1 -T AUTO -s PF04744_pMOB_filter_curated_manual.fasta -m MFP -mset LG,WAG,JTT -cmax 4 --prefix train_pMOB
iqtree3 -seed 1 -T AUTO -s PF04744_pMOB_filter_curated_manual.fasta -te train_pMOB.treefile --init-model LG --model-joint GTR20+FO --prefix pMOB_GTR20_FO
#grep -A 21 "can be used as input for IQ-TREE" AMOA_GTR20_FO.iqtree | tail -n20 > Q.AMOA_tm

#https://iqtree.github.io/doc/Estimating-amino-acid-substitution-models
