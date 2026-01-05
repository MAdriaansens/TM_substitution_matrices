import xml.etree.ElementTree as ET
tree = ET.parse('topdb_all.xml.1')
root = tree.getroot()
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
for child in root:
    tm_count = 0
    tm_sequence = ''
    #child of root is each entry in db
    attribute = str(child.attrib)
    #this makes sure only alpha polytopic and bitopic tm proteins get through
    if 'Alpha' not in attribute:
        pass
    else:
        alpha_count = alpha_count + 1
        if 'Sequence'  in (str((child[1]))):
            pass
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! bug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elif len(child) != 6:
            pass
                        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! bug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else:
            sequence = str(child[2][0].text)
            sequence = (sequence.replace(' ', '').replace('\n', ''))
            pdb_id = str(child[1][0].attrib).split("'")[-2]
            Topology = (child[3])
            for entry in Topology:
                membrane_type = entry.text.split(';')[0]
                list_membr_type.append(membrane_type)
            for topology_annot in child[4][2]:
                annotation = str(topology_annot.attrib)
                if 'Membrane' in annotation:
                    tm_count =  tm_count + 1
                    membrane_annotation = annotation
                    tm_sequence = tm_sequence + (get_membrane_sequence(sequence, membrane_annotation))
            #print('>{}_membrane_type:{}_aTM_count:{}'.format(pdb_id, membrane_type, tm_count))
            #print(tm_sequence)
            total_count = total_count + 1
print(total_count)
print(alpha_count) 
print(set(list_membr_type))
