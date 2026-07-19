import random
from Bio import SeqIO

Protein_id_list = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/Matrice_Qmaker/MonoOxygenase/sequences/FL/Pfam02461_vsCurated_Bacterial_MoA_filtered_curated_final_linsi_FL.fa'
Output_dir = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/Matrice_generation/Matrice_Qmaker/MonoOxygenase/sequences/FL/random_100'
Id_list=[]
Id_dict={}
for record in SeqIO.parse('{}'.format(Protein_id_list), 'fasta'):
    Id_list.append(record.id)
    Id_dict[record.id] = record.seq
for i in range(100):
    if i == 0:
        unique_sample=Id_list[0:100]
        with open('{}/Pfam02461_BMOA_random_{}_fl.fasta'.format(Output_dir, i), 'w') as output:
            for entry in unique_sample:
                sequence = '>{}'.format(entry) + '\n' + Id_dict[entry] + '\n'
                output.write(str(sequence))
    elif i == 100:
        unique_sample_Id_list[-101:-1]
        with open('{}/Pfam02461_BMOA_random_{}_fl.fasta'.format(Output_dir, i), 'w') as output:
            for entry in unique_sample:
                sequence = '>{}'.format(entry) + '\n' + Id_dict[entry] + '\n'
                output.write(str(sequence))
    else:
        random.shuffle(Id_list)
        unique_sample=Id_list[0:100]
        with open('{}/Pfam02461_BMOA_random_{}_fl.fasta'.format(Output_dir, i), 'w') as output:
            for entry in unique_sample:
                sequence = '>{}'.format(entry) + '\n' + Id_dict[entry] + '\n'
                output.write(str(sequence))
