import re
import sys

hit_file = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/protein-matching-BacteriaMO_PF02461_alignedPF02461.sthk'
outfile='/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/protein-matching-BacteriaMO_PF02461_alignedPF02461_TM_alignments_only.fasta'
alignment={}

length=167

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
print(lengthDict)


if len(lengthDict) > 1:
    sys.exit("ERROR! Aligned sequences differ by length")


TM = []
TM.extend(list(range(15,42)))
TM.extend(list(range(52,74)))
TM.extend(list(range(77,102)))
TM.extend(list(range(113,155)))
TM.extend(list(range(187,194)))
TM.extend(list(range(204,207)))
TM.extend(list(range(213,234)))
print(TM)

with open(outfile, 'w') as out:
    for header in alignment.keys():
        if (len(alignment[header].replace('-', ''))) > length:
            TM_align = ''
            residue_num = 1
            for residue in alignment[header]:
                if residue_num in TM:
                    TM_align += residue
                else:
                    TM_align += '-'
                residue_num = residue_num + 1
            line = '>' + header.split('|')[0] + '\n' + TM_align + '\n'
            out.write(line)
        else:
            pass
