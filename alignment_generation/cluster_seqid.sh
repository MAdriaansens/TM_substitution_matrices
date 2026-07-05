#!/bin/bash
#SBATCH --job-name=06450_Homology_Euk_direct
#SBATCH --time=2:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=5GB          # Memory in MB
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_output/06450_Homology_Euk_direct_output%A.out
#SBATCH --error=slurm_output/06450_Homology_Euk_direct_error%A.err


module load MMseqs2/15-6f452-gompi-2023a
#mmseqs easy-cluster  --min-seq-id 0.95 -c 0.0 HMMalign/Pfam04744_vspMOB_filter.fa HMMalign/Pfam04744_vspMOB_filter_clusterd.fa  temp
#rm -r temp

#mmseqs easy-cluster  --min-seq-id 0.95 -c 0.0 HMMalign/Pfam04896_vspMOC_filter.fa HMMalign/Pfam04896_vspMOC_filter_clusterd.fa  temp
#rm -r temp

mmseqs easy-cluster  --min-seq-id 0.95 -c 0.0 HMMalign/Pfam12942_vsArchaeal_MoA_filtered.fa HMMalign/Pfam12942_vsArchaeal_MoA_filtered_clusterd.fa  temp
rm -r temp

mmseqs easy-cluster  --min-seq-id 0.95 -c 0.0 HMMalign/Pfam02461_vsBacterial_MoA_filtered.fa HMMalign/Pfam02461_vsBacterial_MoA_filtered_clusterd.fa  temp
rm -r temp
