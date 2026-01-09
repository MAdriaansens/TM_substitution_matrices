#goal of this code is to parse and categorize families.
#using the ecod database

ecod_info_dict = {}
with open('ecod.v293.domains.txt', 'r') as lines:
    for line in lines:
        if line[0] == '#':
            pass
        elif 'pdb_chain' in line:
            #this is to remove the header
            pass
        else:
            #sometimes  multiple components make up a pdb id so it is split up
            if '-' in line.split('\t')[5]:
                multimer_split = []
                multimer_split = line.split('\t')[5].split('-')
                for mer in multimer_split:
                    pdb_id = line.split('\t')[4] + '_' + mer
                    pdb_id = (line.split('\t')[4] + '_' + line.split('\t')[5])
                
                #ecod has its own family classification, firs part is architecture of structure (beta-barrel), second is structural similarity but little homology support. 
                #third is second but grouped based on topological connections, the last is the family (Often pfam). 
    
                #so for 2nmz_A, its family_id = 1.1.1.3, 1 1 (structural similarity). = cradle loop barrel, 
                #1 expected homology = RIFT-related, 1 topological connection = Acid protease, 3 family = RVP. 
                family_id = line.split('\t')[3]
                #I also append the structural architecture to this.
                Architecture = line.split('\t')[7]
                #a family id 
                f_name = line.split('\t')[11]
                info_list = (pdb_id,family_id, f_name, Architecture)
                ecod_info_dict[pdb_id] = info_list
            #if it just the idenvidual chain then we write it like this:      
            else:
                pdb_id = (line.split('\t')[4] + '_' + line.split('\t')[5])
                
                family_id = line.split('\t')[3]
                
                Architecture = line.split('\t')[7]
                
                f_name = line.split('\t')[11]
                #concatenate all the info and save into dictionairy.
                info_list = (pdb_id,family_id, f_name, Architecture)
                ecod_info_dict[pdb_id] = info_list
#this function (in a laborius manner) calculates the relative frequency for each amino acid 
def return_line_of_frequencies(input_list):
    #the 0.000000000001, addition is done to ensure the calculation still work despite some amino acids possibly being absent
    A = str(input_list.count('A'))
    R = str(input_list.count('R'))
    N = str(input_list.count('N'))
    D = str(input_list.count('D'))
    C = str(input_list.count('C'))
    Q = str(input_list.count('Q'))
    E = str(input_list.count('E'))
    G = str(input_list.count('G'))
    H = str(input_list.count('H'))
    I = str(input_list.count('I'))
    L = str(input_list.count('L'))
    K = str(input_list.count('K'))
    M = str(input_list.count('M'))
    F = str(input_list.count('F'))
    P = str(input_list.count('P')) 
    S = str(input_list.count('S'))
    T = str(input_list.count('T'))
    W = str(input_list.count('W'))
    Y = str(input_list.count('Y'))
    V = str(input_list.count('V'))

        
    line = A  + '\t' +  R  + '\t' +  N  + '\t' + D  + '\t' + \
    C  + '\t' +  Q  + '\t' +  E  + '\t' +  G  + '\t' +  H  + '\t' +  I  + '\t' +  L  + '\t' + \
    K  + '\t' +  M  + '\t' +  F  + '\t' +  P  + '\t' +  S  + '\t' +  T  + '\t' +  W  + '\t' +  Y  + '\t' +  V 
    return(line)


import os
from Bio import SeqIO
import json

writing_dir = '/nesi/nobackup/uc04105/PDB_alpha/tma_topdb/family_sep'
#parse the pdb ids present in the filtered pdbtm database

#read membrane data json file

membrane_data = '/nesi/nobackup/uc04105/PDB_alpha/pdb_chain_id_mebrane_data.json'
with open(membrane_data) as json_data:
    membrane_dict = json.load(json_data)
    json_data.close()

#write a fasta and a tsv output

with open('/nesi/nobackup/uc04105/PDB_alpha/Alpha-Helical_tm_PDBTM_ECOD_family.tsv', 'w') as tsv_out:
    header = 'pdb_id' + '\t' + 'family_name' + '\t' + 'Architecture' + '\t' + 'ecod_code' + '\t' + 'tm_seq' + '\t'+ 'length_atm' + '\t' + 'fl_seq' + '\t' + 'length_fl' + '\t' + 'membrane_type' '\t'  + 'tA'+ '\t' + 'tR' + '\t' + 'tN'+ '\t' + 'tD'+ '\t' + 'tC'+ '\t' + \
'tQ'+ '\t' + 'tE'+ '\t' + 'tG'+ '\t' + 'tH'+ '\t' + 'tI'+ '\t' + 'tL'+ '\t' + 'tK'+ '\t' + 'tM'+ '\t' + \
'tF'+ '\t' +'tP'+ '\t' + 'tS'+ '\t' + 'tT'+ '\t' + 'tW'+ '\t' + 'tY'+ '\t' + 'tV' + '\t' + 'fA'+ '\t' + 'fR' + '\t' + 'fN'+ '\t' + 'fD'+ '\t' + 'fC'+ '\t' + \
'fQ'+ '\t' + 'fE'+ '\t' + 'fG'+ '\t' + 'fH'+ '\t' + 'fI'+ '\t' + 'fL'+ '\t' + 'fK'+ '\t' + 'fM'+ '\t' + \
'fF'+ '\t' +'fP'+ '\t' + 'fS'+ '\t' + 'fT'+ '\t' + 'fW'+ '\t' + 'fY'+ '\t' + 'fV' +'\n'
    tsv_out.write(header)
    for pdb_chain_fasta in os.listdir('/nesi/nobackup/uc04105/PDB_alpha/tma_topdb'):  
        aa_list = []
        if pdb_chain_fasta.split('_tm.fasta')[0] in ecod_info_dict:
            #family_name will be used
            family_name = ecod_info_dict[pdb_chain_fasta.split('_tm.fasta')[0]][2].replace(' ', '_').replace('/', '_slash_')
            ecod_code = ecod_info_dict[pdb_chain_fasta.split('_tm.fasta')[0]][1].replace('.', '_')

            #annotation of the db is not complete some proteins do have a annotaiton but no family annotation, their family name will be their whole id
            if family_name == '':
                family_name =ecod_code
            #write a fasta file containin the tm for each pdb matching it
            with open('/nesi/nobackup/uc04105/PDB_alpha/tma_topdb/family_sep/{}_tm.fasta'.format(family_name), 'a') as output:
                for record in SeqIO.parse('/nesi/nobackup/uc04105/PDB_alpha/tma_topdb/{}'.format(pdb_chain_fasta), 'fasta'):
                    fasta_line = '>' + record.id + ' family_id: {} ecod_id {}:'.format(family_name, ecod_code) + '\n' + str(record.seq) + '\n'
                    tm_seq = str(record.seq)

                    #count the aa in a tm element
                    for aa in tm_seq:
                        aa_list.append(aa)
                    count_aa_tm =  return_line_of_frequencies(aa_list)
                    
                    output.write(fasta_line)
            output.close()
            
            with open('/nesi/nobackup/uc04105/PDB_alpha/fl_topdb/family_sep/{}_fl.fasta'.format(family_name), 'a') as output:
                for record in SeqIO.parse('/nesi/nobackup/uc04105/PDB_alpha/fl_topdb/{}'.format(pdb_chain_fasta.replace('tm', 'fl')), 'fasta'):
                    fasta_line = '>' + record.id + ' family_id: {} ecod_id: {}'.format(family_name, ecod_code) + '\n' + str(record.seq) + '\n'
                    fl_seq = str(record.seq)


                    for aa in fl_seq:
                        aa_list.append(aa)
                    count_aa_fl =  return_line_of_frequencies(aa_list)
                    
                    output.write(fasta_line)
                    
                output.write(fasta_line)
            output.close()
            Architecture = ecod_info_dict[pdb_chain_fasta.split('_tm.fasta')[0]][-1]
            
            if pdb_chain_fasta.split('_tm.fasta')[0] in membrane_dict:
                
                tsv_line = pdb_chain_fasta.split('_tm.fasta')[0] + '\t' + family_name + '\t' + Architecture + '\t' + ecod_code + '\t' + tm_seq + '\t' + str(len(tm_seq)) + '\t' + fl_seq + '\t' + str(len(fl_seq)) + '\t'  + membrane_dict[pdb_chain_fasta.split('_tm.fasta')[0]].replace(' ', '_') + '\t' + count_aa_tm + '\t' + count_aa_fl +'\n'
                print(tsv_line)
            else:
                pass
            tsv_out.write(tsv_line)
