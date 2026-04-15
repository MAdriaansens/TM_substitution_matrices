import re
import sys

hit_file = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/protein-matching-AMO_PF12942_alignedPF12942.sthk'
outfile = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/protein-matching-AMO_PF12942_alignedPF12942_alignment_TMonly.faa'
# define alignment
alignment={}



# read in alignment

with open(hit_file) as hmmalignment:
    for line in hmmalignment:
        if(line[0]!="#"):
            splitLine=line.split()

            if(len(splitLine)==2):
                stringFilter = lambda text: re.sub('[a-z.*]', '', splitLine[1])
                filteredString=stringFilter(splitLine[1])

                if(splitLine[0] in alignment.keys()):
                    alignment[splitLine[0]]=alignment[splitLine[0]]+filteredString
                else:
                    alignment[splitLine[0]]=filteredString



# check for identical length

lengthDict={}

for header in alignment.keys():
    lengthDict[len(alignment[header])]=1



if len(lengthDict) > 1:
    sys.exit("ERROR! Aligned sequences differ by length")
#TM residues of Archaea AMOA after Q04507 through TMHMM2.0 on 15th April 2025
TM = []
TM.extend(list(range(28,50)))
TM.extend(list(range(64,87)))
TM.extend(list(range(91,114)))
TM.extend(list(range(123,143)))
TM.extend(list(range(147,170)))
TM.extend(list(range(218,241)))
print(TM)

with open(outfile, 'w') as out:
    for header in alignment.keys():
        TM_align = ''
        residue_num = 1
        for residue in alignment[header]:
            if residue_num in TM:
                TM_align += residue
            else:
                TM_align += '-'
            residue_num = residue_num + 1
            
        
        line = '>' + header.split('|')[0]+ '\n' + TM_align + '\n'
        print(header.split('|')[0])
        out.write(line)
