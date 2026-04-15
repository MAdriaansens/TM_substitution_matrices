import re
import sys

hit_file = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/AMOA_curated_alignedPF12942.sthk'
outfile='/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/protein-matching-AMO_PF12942_alignedPF12942_alignment_TMonly.faa'
alignment={}

length=128

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
TM.extend(list(range(9,11)))
TM.extend(list(range(27,47)))
TM.extend(list(range(50,71)))
TM.extend(list(range(80,100)))
TM.extend(list(range(105,123)))
print(TM)
print(alignment[header][9:11])
print(alignment[header][27:47])

print(alignment[header][50:71])
print(alignment[header][80:100])
print(alignment[header][105:123])

print(alignment[header])

entry = 'MSIFRTEEILKAAKMPPEAVHMSRLIDAVYFPILIILLVGTYHMHFMLLAGDWDFWMDWKDRQWWPVVTPIVGITYCSAIMYYLWVNYRQPFGATLCVVCLLIGEWLTRYWGFYWWSHYPINFVTPGIMLPGALMLDFTLYLTRNWLVTALVGGGFFGLLFYPGNWPIFGPTHLPIVVEGTLLSMADYMGHLYVRTGTPEYVRHIEQGSLRTFGGHTTVIAAFFSAFVSMLMFTVWWYLGKVYCTAFFYVKGKRGRIVHRNDVTAFGEEGFPEGIK'
residue_num=1
residue  = ''
for a in entry:
    if residue_num in TM:
        residue += a
    else:
        residue += '-'
    residue_num += 1
print(residue)





with open(outfile, 'w') as out:
    for header in alignment.keys():
        print(alignment[header])
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
