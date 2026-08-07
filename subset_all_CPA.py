from Bio import SeqIO
from collections import Counter
gap_list = []
for record in SeqIO.parse('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/parsed_tm_CPA.fasta','fasta'):
    count =0
    for residue in record.seq:
        
        if residue == '-':
            gap_list.append(count)
        count = count + 1
print(gap_list)
print(Counter(gap_list))
my_counter=Counter(gap_list)
to_remove_list=[]
for item, count in my_counter.items():
    if count > 3:
        to_remove_list.append(item)
    else:
        pass
print(to_remove_list)
with open('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/TM/All_CPA_curated_merged_TM_PF00999_filtered_inspected_removed.fa', 'w') as inspected:
    for record in SeqIO.parse('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/FL/All_CPA_curated_merged_FL_PF00999_filtered.fa', 'fasta'):
        result = "".join([char for idx, char in enumerate(str(record.seq)) if idx not in to_remove_list])
        id_record = '>' + record.id + '_nonTM_removed' + '\n'
        Line = (id_record+result + '\n')
        inspected.write(Line)

#TM residues
with open('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/TM/All_Archaea_CPA_curated_merged_TM_PF00999_filtered_inspected_removed.fa', 'w') as inspectedA:
    with open('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/TM/All_Bacteria_CPA_curated_merged_TM_PF00999_filtered_inspected_removed.fa', 'w') as inspectedB:
        with open('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/TM/All_Eukarya_CPA_curated_merged_TM_PF00999_filtered_inspected_removed.fa', 'w') as inspectedE:
            for record in SeqIO.parse('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/TM/All_CPA_curated_merged_TM_PF00999_filtered_inspected_removed.fa', 'fasta'):
                if '>Arc226_' in record.id:
                    Line = '>' + record.id + '\n' + str(record.seq) + '\n'
                    inspectedA.write(Line)
                elif '>Bac226_' in record.id:
                    Line = '>' + record.id + '\n' + str(record.seq) + '\n'
                    inspectedB.write(Line)
                elif '>EukM6_' in record.id:
                    Line = '>' + record.id + '\n' + str(record.seq) + '\n'
                    inspectedE.write(Line)
#FL residues
with open('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/FL/All_Archaea_CPA_curated_merged_FL_PF00999_filtered.fa', 'w') as inspectedA:
    with open('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/FL/All_Bacteria_CPA_curated_merged_FL_PF00999_filtered.fa', 'w') as inspectedB:
        with open('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/FL/All_Eukarya_CPA_curated_merged_FL_PF00999_filtered.fa', 'w') as inspectedE:
            for record in SeqIO.parse('/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/HMMalign/CPA/HMMalign_CPA_Adriaansens/FL/All_CPA_curated_merged_FL_PF00999_filtered.fa', 'fasta'):
                if '>Arc226_' in record.id:
                    Line = '>' + record.id + '\n' + str(record.seq) + '\n'
                    inspectedA.write(Line)
                elif '>Bac226_' in record.id:
                    Line = '>' + record.id + '\n' + str(record.seq) + '\n'
                    inspectedB.write(Line)
                elif '>EukM6_' in record.id:
                    Line = '>' + record.id + '\n' + str(record.seq) + '\n'
                    inspectedE.write(Line)
