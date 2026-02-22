#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      MSA_tgroup
#SBATCH --time          48:00:00
#SBATCH --mem           50GB
#SBATCH --cpus-per-task 10
#SBATCH --error         slurm_msa/MSA_run_%A.err
#SBATCH --output        slurm_msa/MSA_run_%A.out
module load MAFFT/7.505-gimkl-2022a-with-extensions
STRING=">"
COUNT=1
for file in /nesi/nobackup/uc04105/PDB_alpha/tma_topdb/tgroup_sep/*.fasta; do
    # Count occurrences
    matches=$(grep -o "$STRING" "$file" | wc -l)
    # If matches equal count, remove
    if [ "$matches" -eq "$COUNT" ]; then
        echo "$file not used in MSA"
    else
        echo "$file is used in MSA"
        partname=${file%.fasta}
        BASENAME=$(basename "$partname")
        echo "$BASENAME"
        mafft --thread 10 --maxiterate 1000 --auto ${file} > ${BASENAME}_autoaligned.fasta
    fi
done
