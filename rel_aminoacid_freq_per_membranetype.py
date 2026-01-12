import json
from Bio import SeqIO


import os
header = 'Name' + '\t' + 'total_aa_count' + '\t' + 'type' + '\t' + 'A'+ '\t' + 'R' + '\t' + 'N'+ '\t' + 'D'+ '\t' + 'C'+ '\t' + \
'Q'+ '\t' + 'E'+ '\t' + 'G'+ '\t' + 'H'+ '\t' + 'I'+ '\t' + 'L'+ '\t' + 'K'+ '\t' + 'M'+ '\t' + \
'F'+ '\t' +'P'+ '\t' + 'S'+ '\t' + 'T'+ '\t' + 'W'+ '\t' + 'Y'+ '\t' + 'V' + '\n'

fl_dir = '/nesi/nobackup/uc04105/PDB_alpha/fl_topdb'
atm_dir = '/nesi/nobackup/uc04105/PDB_alpha/tma_topdb'


#this function (in a laborius manner) calculates the relative frequency for each amino acid 
def return_line_of_frequencies(input_list, name):
    #the 1, addition is done to ensure the calculation still work despite some amino acids possibly being absent
    total = len(input_list)
    if total != 0:
        A = str((input_list.count('A') + 1))
        R = str((input_list.count('R') + 1))
        N = str((input_list.count('N') + 1))
        D = str((input_list.count('D') + 1))
        C = str((input_list.count('C') + 1))
        Q = str((input_list.count('Q') + 1))
        E = str((input_list.count('E') + 1))
        G = str((input_list.count('G') + 1))
        H = str((input_list.count('H') + 1))
        I = str((input_list.count('I') + 1))
        L = str((input_list.count('L') + 1))
        K = str((input_list.count('K') + 1))
        M = str((input_list.count('M') + 1))
        F = str((input_list.count('F') + 1))
        P = str((input_list.count('P') + 1))
        S = str((input_list.count('S') + 1))
        T = str((input_list.count('T') + 1))
        W = str((input_list.count('W') + 1))
        Y = str((input_list.count('Y') + 1))
        V = str((input_list.count('V') + 1))
    else:
        A  =  R  = N  = D  =  C  =  Q  =  E  = G  =  H  =  I  = L  = K =  M  = F  =  P =  S =  T  = W =  Y  =  V  = 'NaN'
        
    line = name + '\t' + str(total) +  '\t' +  name.split('_')[-1] + '\t' + A  + '\t' +  R  + '\t' +  N  + '\t' + D  + '\t' + \
    C  + '\t' +  Q  + '\t' +  E  + '\t' +  G  + '\t' +  H  + '\t' +  I  + '\t' +  L  + '\t' + \
    K  + '\t' +  M  + '\t' +  F  + '\t' +  P  + '\t' +  S  + '\t' +  T  + '\t' +  W  + '\t' +  Y  + '\t' +  V  + '\n'
    return(line)

#subset data in membrane localizations

with open('pdb_chain_id_mebrane_data.json', 'r', encoding='utf-8') as file:
        # Use json.load() to parse the file content into a Python object
        membrane_data = json.load(file)

from collections import Counter

def subtract_with_counts(list_a, list_b):
    counts = Counter(list_b)
    result = []
    for item in list_a:
        if counts[item] > 0:
            counts[item] -= 1
        else:
            result.append(item)
    return result


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
        
        #open the right key in the membranefrom collections import Counter

fl_aa_list = []
#merge all the fl files

for record in SeqIO.parse('{}/all_fl_merged.faa'.format(fl_dir), 'fasta'):
    for aa in str(record.seq):
        fl_aa_list.append(aa)

#calculate the amino acid frequency
name = 'merged_fl'
fl_line =  return_line_of_frequencies(fl_aa_list, name)

atm_aa_list = []

#merge all the atm files

for record in SeqIO.parse('{}/all_atm_merged.faa'.format(atm_dir), 'fasta'):
    for aa in str(record.seq):
        atm_aa_list.append(aa)


#calculate the amino acid frequency

name = 'merged_tm'
atm_line = (return_line_of_frequencies(atm_aa_list, name))


#subtract lists from eachother
non_membrane_regions_aa_list = subtract_with_counts(fl_aa_list, atm_aa_list)
print(len(non_membrane_regions_aa_list))
print(len(fl_aa_list))
print(len(atm_aa_list))


name = 'merged_nontm'
non_mem_line = (return_line_of_frequencies(non_membrane_regions_aa_list, name))


import os
with open('relative_freq_aa_atm_fl.tsv', 'w') as X:
    X.write(header)
    X.write(fl_line)
    X.write(atm_line)
    X.write(non_mem_line)
    
    for key in membrane_freqfl_library.keys():
        name = key + '_fl'
        input_fl_list = membrane_freqfl_library[key]
        X.write(return_line_of_frequencies(input_fl_list, name))
        
        name = key + '_tm'
        input_atm_list = membrane_freqtm_library[key]
        X.write(return_line_of_frequencies(input_atm_list, name))

        name = key + '_nontm'
        non_tm_list = subtract_with_counts(input_fl_list, input_atm_list)
        X.write(return_line_of_frequencies(non_tm_list, name))
