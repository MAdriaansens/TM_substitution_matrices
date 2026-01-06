import json
from Bio import SeqIO



header = 'Name' + '\t' + 'total_aa_count' + '\t' + 'A'+ '\t' + 'R' + '\t' + 'N'+ '\t' + 'D'+ '\t' + 'C'+ '\t' + \
'Q'+ '\t' + 'E'+ '\t' + 'G'+ '\t' + 'H'+ '\t' + 'I'+ '\t' + 'L'+ '\t' + 'K'+ '\t' + 'M'+ '\t' + \
'F'+ '\t' +'P'+ '\t' + 'S'+ '\t' + 'T'+ '\t' + 'W'+ '\t' + 'Y'+ '\t' + 'V' + '\n'

print(len(amino_acids_single))
fl_dir = '/nesi/nobackup/uc04105/PDB_alpha/fl_topdb'
atm_dir = '/nesi/nobackup/uc04105/PDB_alpha/tma_topdb'


#this function (in a laborius manner) calculates the relative frequency for each amino acid 
def return_line_of_frequencies(input_list, name):
    #the 0.000000000001, addition is done to ensure the calculation still work despite some amino acids possibly being absent
    total = len(input_list)
    if total != 0:
        A = str((input_list.count('A') + 0.000000000001)/total)
        R = str((input_list.count('R') + 0.000000000001)/total)
        N = str((input_list.count('N') + 0.000000000001)/total)
        D = str((input_list.count('D') + 0.000000000001)/total)
        C = str((input_list.count('C') + 0.000000000001)/total)
        Q = str((input_list.count('Q') + 0.000000000001)/total)
        E = str((input_list.count('E') + 0.000000000001)/total)
        G = str((input_list.count('G') + 0.000000000001)/total)
        H = str((input_list.count('H') + 0.000000000001)/total)
        I = str((input_list.count('I') + 0.000000000001)/total)
        L = str((input_list.count('L') + 0.000000000001)/total)
        K = str((input_list.count('K') + 0.000000000001)/total)
        M = str((input_list.count('M') + 0.000000000001)/total)
        F = str((input_list.count('F') + 0.000000000001)/total)
        P = str((input_list.count('P') + 0.000000000001)/total)
        S = str((input_list.count('S') + 0.000000000001)/total)
        T = str((input_list.count('T') + 0.000000000001)/total)
        W = str((input_list.count('W') + 0.000000000001)/total)
        Y = str((input_list.count('Y') + 0.000000000001)/total)
        V = str((input_list.count('V') + 0.000000000001)/total)
    else:
        A  =  R  = N  = D  =  C  =  Q  =  E  = G  =  H  =  I  = L  = K =  M  = F  =  P =  S =  T  = W =  Y  =  V  = 'NaN'
        
    line = name + '\t' + str(total) +  '\t' + A  + '\t' +  R  + '\t' +  N  + '\t' + D  + '\t' + \
    C  + '\t' +  Q  + '\t' +  E  + '\t' +  G  + '\t' +  H  + '\t' +  I  + '\t' +  L  + '\t' + \
    K  + '\t' +  M  + '\t' +  F  + '\t' +  P  + '\t' +  S  + '\t' +  T  + '\t' +  W  + '\t' +  Y  + '\t' +  V  + '\n'
    return(line)

#subset data in membrane localizations

with open('pdb_chain_id_mebrane_data.json', 'r', encoding='utf-8') as file:
        # Use json.load() to parse the file content into a Python object
        membrane_data = json.load(file)



#make a set list of all membrane types
membrane_type_list = []
for key in membrane_data.keys():
    membrane_type_list.append(membrane_data[key])

from collections import Counter
print(Counter(membrane_type_list))

#make a count dictionary
Count_dic = dict(Counter(membrane_type_list))

membrane_type_list = set(membrane_type_list)

#make an empty dictionairy
#this dictionairy will contain a list of amino acids, with each amino acid coming from a sequence from that membrane type
membrane_freqtm_library = {}
for type_membrane in membrane_type_list:
    membrane_freqtm_library[type_membrane] = []

#first make sure the files are in the dictionary, not all are due to annotation issues from topdb
for file in os.listdir(atm_dir):
    if file.split('_tm.fasta')[0] in membrane_data:
        #get membrane type per file
        membrane_type = (membrane_data[file.split('_tm.fasta')[0]])


        #get the sequence
        for record in SeqIO.parse('{}/{}'.format(atm_dir, file), 'fasta'):
            #should be one sequence
            for i in str(record.seq):
                membrane_freqtm_library[membrane_type].append(i)
        
        #open the right key in the membrane_freq_library dictionary
        

#a seperate one for fl 
membrane_freqfl_library = {}
for type_membrane in membrane_type_list:
    membrane_freqfl_library[type_membrane] = []

#first make sure the files are in the dictionary, not all are due to annotation issues from topdb
for file in os.listdir(fl_dir):
    if file.split('_fl.fasta')[0] in membrane_data:
        #get membrane type per file
        membrane_type = (membrane_data[file.split('_fl.fasta')[0]])


        #get the sequence
        for record in SeqIO.parse('{}/{}'.format(fl_dir, file), 'fasta'):
            #should be one sequence
            for i in str(record.seq):
                membrane_freqfl_library[membrane_type].append(i)
        
        #open the right key in the membrane_freq_library dictionary
fl_aa_list = []
#merge all the fl files

for record in SeqIO.parse('{}/all_fl_merged.faa'.format(fl_dir), 'fasta'):
    for aa in str(record.seq):
        fl_aa_list.append(aa)

#calculate the amino acid frequency
name = 'fl_tm_proteins'
fl_line =  return_line_of_frequencies(fl_aa_list, name)

atm_aa_list = []

#merge all the atm files

for record in SeqIO.parse('{}/all_atm_merged.faa'.format(atm_dir), 'fasta'):
    for aa in str(record.seq):
        atm_aa_list.append(aa)

#calculate the amino acid frequency

name = 'atm_tm_proteins'
atm_line = (return_line_of_frequencies(atm_aa_list, name))

import os
with open('relative_freq_aa_atm_fl.tsv', 'w') as X:
    X.write(header)
    X.write(fl_line)
    X.write(atm_line)

    
    for key in membrane_freqfl_library.keys():
        name = key + '_fl'
        input_list = membrane_freqfl_library[key]
        X.write(return_line_of_frequencies(input_list, name))
    for key in membrane_freqtm_library.keys():
        name = key + '_tm'
        input_list = membrane_freqtm_library[key]
        X.write(return_line_of_frequencies(input_list, name))
