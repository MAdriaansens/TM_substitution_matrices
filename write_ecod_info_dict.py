#goal of this code is to parse and categorize families.
#using the ecod database
import json

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
                X_group=line.split('\t')[8]
                H_group=line.split('\t')[9]
                T_group=line.split('\t')[10]
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
                
                Architecture = line.split('\t')[8]
                X_group=line.split('\t')[9]
                H_group=line.split('\t')[10]
                T_group=line.split('\t')[11]
                
                f_name = line.split('\t')[12]
                #concatenate all the info and save into dictionairy.
                info_list = (pdb_id,family_id, f_name, Architecture)
                ecod_info_dict[pdb_id] = info_list
with open('ecod_data.json', 'w') as f:
    json.dump(ecod_info_dict, f)
