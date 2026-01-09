import xml.etree.ElementTree as ET
#read and install the modules for .xml reading
tree = ET.parse('topdb_all.xml.1')
root = tree.getroot()

#output directories and output json for membrane types
writing_dir_fl = '/nesi/nobackup/uc04105/PDB_alpha/fl_topdb'
writing_dir_tm = '/nesi/nobackup/uc04105/PDB_alpha/tma_topdb'
membrane_type_dict = {}

#what this function does is that it reads the sequence and the topology data from xml and then parses out the TM sequence
def get_membrane_sequence(sequence,membrane_annotation):
    #keep in mind that the pdb starts counting at 1 but python counts at 0 and you would need to include the first residue of the count as well
    #so interval should be [begin-1: end]
    begin = int(membrane_annotation.split(" '")[1].split("',")[0])-1
    end = int(membrane_annotation.split(" '")[3].split("',")[0])
    part_sequence = sequence[begin:end]
    return(part_sequence)



print(len(root))
alpha_count = 0
total_count = 0
list_membr_type = []
id_list = []
#child of root is each entry in db
tm_count = 0

#the code cant tell apart chains yet


membrane_data = '/nesi/nobackup/uc04105/PDB_alpha/pdb_chain_id_mebrane_data.json'
with open('ecod_data.json') as json_data:
    ecod_info_dict = json.load(json_data)
    json_data.close()
beta_list =['beta meanders', 'beta complex topology', 'beta barrels', 'beta sandwiches' , 'beta duplicates or obligate multimers', '']


for child in root:
    membrane_protein_flag = False
    
    tm_sequence = ''
    
    attribute = str(child.attrib)
    
    #this makes sure only alpha polytopic and bitopic tm proteins get through
    if 'Alpha' not in attribute:
        pass
    else:
        alpha_count = alpha_count + 1
        if 'Sequence'  in (str((child[1]))):
            pass
            
            #everything where sequence is child[1] misses a pdb annotation. 
            #that means we cannot assign the correct family to them based on literature
            #furthermore I cannot always assign who identified this proteins since it missess all other data
        elif len(child) == 10:
            pass
            #only two match the >10 and both are antibodies so they will be exculded
        #this is in here to alarm you if any pdb entries are missed by the code
        elif len(child) > 8:
            print('present')
        elif len(child) == 8:
            pdb_id = str(child[1][0].attrib).split("'")[-2]
            Topology = (child[3])
            chain_id = str(child[1][0][2][0].attrib).split("'")[-2]
            pdb_id = pdb_id + '_' + chain_id

            #make sure it is not an ecod
            if pdb_id in ecod_info_dict:
                if ecod_info_dict[pdb_id][-1] not in beta_list:
    
                    #small edit to topology annotation here
                    Topology_annotation = (child[6][2])
        
                    sequence = str(child[2][0].text)
                    sequence = (sequence.replace(' ', '').replace('\n', ''))
                    
                    tm_count = tm_count + 1
                    
                    #from the membrane annotation we now retrieve the start and end sites of the membrane 
                    for entry in Topology:
                        membrane_type = entry.text.split(';')[0]
                        list_membr_type.append(membrane_type)
                        membrane_type_dict[pdb_id] = membrane_type
        
                    
                    
                    for topology_annot in Topology_annotation:
                        annotation = str(topology_annot.attrib)
                        if 'Membrane' in annotation:
                            
                            membrane_protein_flag = True
                            membrane_annotation = annotation
                            tm_sequence = tm_sequence + (get_membrane_sequence(sequence, membrane_annotation))
                            #print('>{}_membrane_type:{}_aTM_count:{}'.format(pdb_id, membrane_type, tm_count))
                    #print(len(tm_sequence), len(sequence), len(tm_sequence)/len(sequence))
        
                    if membrane_protein_flag == True:
                        #write fl file
                        id_list.append(pdb_id)
                        fasta_header_fl = '>{}_fl_membranetype_{}'.format(pdb_id, membrane_type) 
                        fasta_output = fasta_header_fl + '\n' + sequence + '\n'
                        with open('{}/{}_fl.fasta'.format(writing_dir_fl, pdb_id), 'w') as ffl: #Fasta Full Length (FFL)
                            ffl.write(fasta_output)
                        ffl.close()
        
                        #write atm
                        fasta_header_fl = '>{}_tm_membranetype_{}'.format(pdb_id, membrane_type) 
                        fasta_output = fasta_header_fl + '\n' + tm_sequence + '\n'
                        with open('{}/{}_tm.fasta'.format(writing_dir_tm, pdb_id), 'w') as atm: #Alpha tm (Atm)
                            atm.write(fasta_output)
                        atm.close()

        
        elif len(child) == 7:
            #7 are multimeric proteins, or homodimers
            pdb_id = str(child[1][0].attrib).split("'")[-2]
            Topology = (child[3])
            chain_id = str(child[1][0][2][0].attrib).split("'")[-2]
            pdb_id = pdb_id + '_' + chain_id
            if pdb_id in ecod_info_dict:
                if ecod_info_dict[pdb_id][-1] not in beta_list:

                    #small edit to topology annotation here
                    Topology_annotation = (child[5][2])
            
        
                    sequence = str(child[2][0].text)
                    sequence = (sequence.replace(' ', '').replace('\n', ''))
                    
                    tm_count = tm_count + 1
        
                    #from the membrane annotation we now retrieve the start and end sites of the membrane 
                    for entry in Topology:
                        membrane_type = entry.text.split(';')[0]
                        list_membr_type.append(membrane_type)
                        membrane_type_dict[pdb_id] = membrane_type
                
        
                    
                    for topology_annot in Topology_annotation:
                        annotation = str(topology_annot.attrib)
                        if 'Membrane' in annotation:
                            membrane_protein_flag = True
                            membrane_annotation = annotation
                            tm_sequence = tm_sequence + (get_membrane_sequence(sequence, membrane_annotation))
                            #print('>{}_membrane_type:{}_aTM_count:{}'.format(pdb_id, membrane_type, tm_count))
                    #print(len(tm_sequence), len(sequence), len(tm_sequence)/len(sequence))
        
                    if membrane_protein_flag == True:
                        #write fl file
                        id_list.append(pdb_id)
        
                        fasta_header_fl = '>{}_fl_membranetype_{}'.format(pdb_id, membrane_type) 
                        fasta_output = fasta_header_fl + '\n' + sequence + '\n'
                        with open('{}/{}_fl.fasta'.format(writing_dir_fl, pdb_id), 'w') as ffl: #Fasta Full Length (FFL)
                            ffl.write(fasta_output)
                        ffl.close()
        
                        #write atm
                        fasta_header_fl = '>{}_tm_membranetype_{}'.format(pdb_id, membrane_type) 
                        fasta_output = fasta_header_fl + '\n' + tm_sequence + '\n'
                        with open('{}/{}_tm.fasta'.format(writing_dir_tm, pdb_id), 'w') as atm: #Alpha tm (Atm)
                            atm.write(fasta_output)
                        atm.close()
                
        elif len(child) == 6:
            pdb_id = str(child[1][0].attrib).split("'")[-2]
            Topology = (child[3])
            chain_id = str(child[1][0][2][0].attrib).split("'")[-2]
            pdb_id = pdb_id + '_' + chain_id
            if pdb_id in ecod_info_dict:
                if ecod_info_dict[pdb_id][-1] not in beta_list:
    
                    Topology_annotation = (child[4][2])
                    #I have no clue why the chunk below seems to work, it is identical to the chunck below but without it the pdb_id isnt recognized
                    #but using this it parses pdbs with len(child[0]) = 4 and 5, without issues. 
                    sequence = str(child[2][0].text)
                    sequence = (sequence.replace(' ', '').replace('\n', ''))
        
                    tm_count = tm_count + 1
                    for entry in Topology:
                        membrane_type = entry.text.split(';')[0]
                        list_membr_type.append(membrane_type)
                        membrane_type_dict[pdb_id] = membrane_type
                            #from the membrane annotation we now retrieve the start and end sites of the membrane 
        
                                #write fl file
                    fasta_header_fl = '>{}_fl_membranetype_{}'.format(pdb_id, membrane_type) 
                    fasta_output = fasta_header_fl + '\n' + sequence + '\n'
                    
                    for topology_annot in Topology_annotation:
                        annotation = str(topology_annot.attrib)
                        if 'Membrane' in annotation:
                            membrane_protein_flag = True
                            membrane_annotation = annotation
                            tm_sequence = tm_sequence + (get_membrane_sequence(sequence, membrane_annotation))
                            #print('>{}_membrane_type:{}_aTM_count:{}'.format(pdb_id, membrane_type, tm_count))
                    #print(len(tm_sequence), len(sequence), len(tm_sequence)/len(sequence))
        
                    
                    if membrane_protein_flag == True:
                        #write fl file
                        id_list.append(pdb_id)
                        fasta_header_fl = '>{}_fl_membranetype_{}'.format(pdb_id, membrane_type) 
                        fasta_output = fasta_header_fl + '\n' + sequence + '\n'
                        with open('{}/{}_fl.fasta'.format(writing_dir_fl, pdb_id), 'w') as ffl: #Fasta Full Length (FFL)
                            ffl.write(fasta_output)
                        ffl.close()
        
                        #write atm
                        fasta_header_fl = '>{}_tm_membranetype_{}'.format(pdb_id, membrane_type) 
                        fasta_output = fasta_header_fl + '\n' + tm_sequence + '\n'
                        with open('{}/{}_tm.fasta'.format(writing_dir_tm, pdb_id), 'w') as atm: #Alpha tm (Atm)
                            atm.write(fasta_output)
                        atm.close()
                
            #if len(tm_sequence)/len(sequence) > 0.5:
            #    print(pdb_id, (len(tm_sequence)/len(sequence)), len(sequence), len(tm_sequence))
        #still need to fix 4, 7 and 10
        elif len(child) == 5:
            pdb_id = str(child[1][0].attrib).split("'")[-2]
            tm_count = tm_count + 1
            #no membrane localization is added the child len = 5
            Topology = (child[3])
            chain_id = str(child[1][0][2][0].attrib).split("'")[-2]
            pdb_id = pdb_id + '_' + chain_id
            if pdb_id in ecod_info_dict:
                if ecod_info_dict[pdb_id][-1] not in beta_list:
    
                    if (len(child[4])) <=  2:
                        #this  means no transmembrane elemtns
                        pass
                    else:
                        sequence = str(child[2][0].text)
                        sequence = (sequence.replace(' ', '').replace('\n', ''))
                      #  print(child[4][1].text)
                      #  print(child[4][0].text)
                        Topology_annotation = (child[4])
        
                        for topology_annot in Topology_annotation:
                            annotation = str(topology_annot.attrib)
                            if 'Membrane' in annotation:
                                membrane_protein_flag = True
                                membrane_annotation = annotation
                                tm_sequence = tm_sequence + (get_membrane_sequence(sequence, membrane_annotation))
                                #print('>{}_membrane_type:{}_aTM_count:{}'.format(pdb_id, membrane_type, tm_count))
                        #print(len(tm_sequence), len(sequence), len(tm_sequence)/len(sequence))
        
                        if membrane_protein_flag == True:
                            #write fl file
                            id_list.append(pdb_id)
                            fasta_header_fl = '>{}_fl'.format(pdb_id) 
                            fasta_output = fasta_header_fl + '\n' + sequence + '\n'
                            with open('{}/{}_fl.fasta'.format(writing_dir_fl, pdb_id), 'w') as ffl: #Fasta Full Length (FFL)
                                ffl.write(fasta_output)
                            ffl.close()
            
                            #write atm
                            fasta_header_fl = '>{}_tm'.format(pdb_id) 
                            fasta_output = fasta_header_fl + '\n' + tm_sequence + '\n'
                            with open('{}/{}_tm.fasta'.format(writing_dir_tm, pdb_id), 'w') as atm: #Alpha tm (Atm)
                                atm.write(fasta_output)
                            atm.close()
                
        #this is in here to alarm you if any pdb entries are missed by the code
        elif len(child) < 5:
            print('present')
        else:
            pass


       
    total_count = total_count + 1
print(total_count)
print(alpha_count) 

print(set(list_membr_type))
print(tm_count)
print(len(list_membr_type))

from collections import Counter
Counter(list_membr_type)
print(len(id_list))
print(len(set(id_list)))
import json
with open('pdb_chain_id_mebrane_data.json', 'w') as f:
    json.dump(membrane_type_dict, f)
