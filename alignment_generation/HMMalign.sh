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

HMMdir=/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/PFAM
Curated_seq=/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/fasta_curated_sequences

#PF02461 BmoA
#hmmalign --amino --trim -o ${HMMalign}/Pfam02461_vsBacterial_MoA.sthk ${HMMdir}/PF02461.hmm  ${HMMdir}/protein-matching-BmoA_PF02461.fasta
#python parse_stockholm_filter.py ${HMMalign}/Pfam02461_vsBacterial_MoA.sthk ${HMMalign}/Pfam02461_vsBacterial_MoA_filtered 183

#hmmalign --amino --trim -o ${HMMalign}/Pfam02461_vsCurated_Bacterial_MoA.sthk ${HMMdir}/PF02461.hmm  ${Curated_seq}/pMOA_curated.fasta
#python parse_stockholm_filter.py  ${HMMalign}/Pfam02461_vsCurated_Bacterial_MoA.sthk  ${HMMalign}/Pfam02461_vsCurated_Bacterial_MoA_filtered 0

#PF12942 AmoA
#hmmalign --amino --trim -o ${HMMalign}/Pfam12942_vsArchaeal_MoA.sthk ${HMMdir}/PF12942.hmm ${HMMdir}/protein-matching-AMO_PF12942.fasta
#python parse_stockholm_filter.py ${HMMalign}/Pfam12942_vsArchaeal_MoA.sthk ${HMMalign}/Pfam12942_vsArchaeal_MoA_filtered 167

#hmmalign --amino --trim -o ${HMMalign}/Pfam12942_vsCurated_Archaeal_MoA.sthk ${HMMdir}/PF12942.hmm ${Curated_seq}/AMOA_curated.fasta
#python parse_stockholm_filter.py ${HMMalign}/Pfam12942_vsCurated_Archaeal_MoA.sthk ${HMMalign}/Pfam12942_vsCurated_Archaeal_MoA_filtered 0

#PF04744 moB
#hmmalign --amino --trim -o ${HMMalign}/Pfam04744_vspMOB.sthk ${HMMdir}/PF04744.hmm ${HMMdir}/protein-matching-moB_PF04744.fasta
#python parse_stockholm_filter.py ${HMMalign}/Pfam04744_vspMOB.sthk ${HMMalign}/Pfam04744_vspMOB_filter 266

#hmmalign --amino --trim -o ${HMMalign}/Pfam04744_vs_curated_pMOB.sthk ${HMMdir}/PF04744.hmm ${Curated_seq}/pMOB_curated.fasta
python parse_stockholm_filter.py ${HMMalign}/Pfam04744_vs_curated_pMOB.sthk  ${HMMalign}/Pfam04744_vs_curated_pMOB_filter 0

#PF04896 moC
#hmmalign --amino --trim -o ${HMMalign}/Pfam04896_vspMOC.sthk ${HMMdir}/PF04896.hmm /home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/PFAM/protein-matching-moC_PF04896.fasta
#python parse_stockholm_filter.py  ${HMMalign}/Pfam04896_vspMOC.sthk ${HMMalign}/Pfam04896_vspMOC_filter  172

hmmalign --amino --trim -o ${HMMalign}/Pfam04896_vs_curated_pMOC.sthk ${HMMdir}/PF04896.hmm ${Curated_seq}/pMOC_curated.fasta
python parse_stockholm_filter.py ${HMMalign}/Pfam04896_vs_curated_pMOC.sthk ${HMMalign}/Pfam04896_vs_curated_pMOC_filter 0

#PF00999
#hmmalign --amino --trim -o ${HMMalign}/Pfam0999_vs_curatedCPA.sthk ${HMMdir}/PF00999.hmm ${Curated_seq}/CPA_curated.fasta
#python parse_stockholm_filter.py ${HMMalign}/Pfam0999_vs_curatedCPA.sthk ${HMMalign}/Pfam0999_vs_curatedCPA_filter 0
