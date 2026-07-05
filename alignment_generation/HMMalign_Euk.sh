#!/bin/bash
#SBATCH --job-name=06450_Homology_Euk_direct
#SBATCH --time=2:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=2GB          # Memory in MB
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_output/06450_Homology_Euk_direct_output%A.out
#SBATCH --error=slurm_output/06450_Homology_Euk_direct_error%A.err

module load Python/3.11.6-foss-2023a
module load HMMER/3.3.2-GCC-12.3.0
HMMalign=/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign

HMMdir=/nesi/nobackup/uc04105/PDB_alpha/Matrice_generation/PFAM
Curated_seq=/nesi/nobackup/uc04105/PDB_alpha/Matrice_generation/CPA_inputseq

#CPA Bacteria
hmmalign --amino --trim -o ${HMMalign}/CPA_bacteria_vsPF00999.sthk ${HMMdir}/PF00999.hmm  ${Curated_seq}/Bacteria_passed_all_filters_OKT3_GTDBreps_alignedPF00999.fasta

python parse_stockholm_filter.py ${HMMalign}/CPA_bacteria_vsPF00999.sthk ${HMMalign}/CPA_bacteria_vsPF00999_filtered 0
#CPA Archaea
hmmalign --amino --trim -o ${HMMalign}/CPA_archaea_vsPF00999.sthk ${HMMdir}/PF00999.hmm  ${Curated_seq}/Archaea_passed_all_filters_7OKT_GTDBreps_alignedPF00999.fasta

python parse_stockholm_filter.py ${HMMalign}/CPA_archaea_vsPF00999.sthk ${HMMalign}/CPA_archaea_vsPF00999_filtered 0

#CPA Eukarya
hmmalign --amino --trim -o ${HMMalign}/CPA_eukarya_vsPF00999.sthk ${HMMdir}/PF00999.hmm  ${Curated_seq}/Eukarya_4nov_passed_all_setcpa.fasta

python parse_stockholm_filter.py ${HMMalign}/CPA_Eukarya_vsPF00999.sthk ${HMMalign}/CPA_Eukarya_vsPF00999_filtered 0
