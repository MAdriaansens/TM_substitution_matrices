import xml.etree.ElementTree as ET
#read and install the modules for .xml reading
tree = ET.parse('topdb_all.xml.1')
root = tree.getroot()
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

#child of root is each entry in db
tm_count = 0

#the code cant tell apart chains yet

for child in root:
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
            #small edit to topology annotation here
            Topology_annotation = (child[6][2])

            sequence = str(child[2][0].text)
            sequence = (sequence.replace(' ', '').replace('\n', ''))
            tm_count = tm_count + 1
            for entry in Topology:
                membrane_type = entry.text.split(';')[0]
                list_membr_type.append(membrane_type)
        
            #from the membrane annotation we now retrieve the start and end sites of the membrane 
            for topology_annot in Topology_annotation:
                annotation = str(topology_annot.attrib)
                if 'Membrane' in annotation:
                    membrane_annotation = annotation
                    tm_sequence = tm_sequence + (get_membrane_sequence(sequence, membrane_annotation))
                    #print('>{}_membrane_type:{}_aTM_count:{}'.format(pdb_id, membrane_type, tm_count))
            #print(len(tm_sequence), len(sequence), len(tm_sequence)/len(sequence))
        elif len(child) == 7:
            #7 are multimeric proteins, or homodimers
            pdb_id = str(child[1][0].attrib).split("'")[-2]
            Topology = (child[3])

            #small edit to topology annotation here
            Topology_annotation = (child[5][2])


            sequence = str(child[2][0].text)
            sequence = (sequence.replace(' ', '').replace('\n', ''))
            tm_count = tm_count + 1
            for entry in Topology:
                membrane_type = entry.text.split(';')[0]
                list_membr_type.append(membrane_type)
        
            #from the membrane annotation we now retrieve the start and end sites of the membrane 
            for topology_annot in Topology_annotation:
                annotation = str(topology_annot.attrib)
                if 'Membrane' in annotation:
                    membrane_annotation = annotation
                    tm_sequence = tm_sequence + (get_membrane_sequence(sequence, membrane_annotation))
                    #print('>{}_membrane_type:{}_aTM_count:{}'.format(pdb_id, membrane_type, tm_count))
            #print(len(tm_sequence), len(sequence), len(tm_sequence)/len(sequence))
        elif len(child) == 6:
            pdb_id = str(child[1][0].attrib).split("'")[-2]
            Topology = (child[3])
            chain_id = str(child[1][0][2][0].attrib).split("'")[-2]
            pdb_id = pdb_id + '_' + chain_id
            print(pdb_id)

            Topology_annotation = (child[4][2])
            #I have no clue why the chunk below seems to work, it is identical to the chunck below but without it the pdb_id isnt recognized
            #but using this it parses pdbs with len(child[0]) = 4 and 5, without issues. 
            sequence = str(child[2][0].text)
            sequence = (sequence.replace(' ', '').replace('\n', ''))

            tm_count = tm_count + 1
            for entry in Topology:
                membrane_type = entry.text.split(';')[0]
                list_membr_type.append(membrane_type)
        
                    #from the membrane annotation we now retrieve the start and end sites of the membrane 
                for topology_annot in Topology_annotation:
                    annotation = str(topology_annot.attrib)
                    if 'Membrane' in annotation:
                        membrane_annotation = annotation
                        tm_sequence = tm_sequence + (get_membrane_sequence(sequence, membrane_annotation))
                    #print('>{}_membrane_type:{}_aTM_count:{}'.format(pdb_id, membrane_type, tm_count))
                    #print(tm_sequence)
            if len(tm_sequence)/len(sequence) > 0.5:
                print(pdb_id, (len(tm_sequence)/len(sequence)), len(sequence), len(tm_sequence))
        #still need to fix 4, 7 and 10
        elif len(child) == 5:
            pdb_id = str(child[1][0].attrib).split("'")[-2]
            tm_count = tm_count + 1
            #no membrane localization is added the child len = 5
            Topology = (child[3])

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
                        membrane_annotation = annotation
                        tm_sequence = tm_sequence + (get_membrane_sequence(sequence, membrane_annotation))
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
