header = '>3RFR_4|Chains G[auth C], H[auth G], I[auth K]|PmoC|Methylocystis sp. M (51782)'
sequence_pmoC = 'MSSTTSTAAGAAAEVESVVDLRGMWIGLAVLNVFYLIVRIYEQVFGWRAGLDSFAPEFQTYWMSILWTEIPLELVSGLGLAGYLWKTRDRNVDAVAPREEMRRLVVLVQWLVVYGIAIYWGASFFTEQDGAWHMTVIRDTDFTPSHIIEFYMSYPIYSVIAVGAFFYAKTRIPYFAHGYSLAFLIVAIGPFMIIPNVGLNEWGHTFWFMEELFVAPLHWGFVFFGWMALGVFGVVLQILGRIHALIGKEGVALLTE'



pmoC_tm_seq_3RFR = 20*gap +  sequence_pmoC[20:46] +  11*gap + sequence_pmoC[57:84] + gap*23 + sequence_pmoC[107:133] + gap*7 + sequence_pmoC[140:167] + gap*12 + sequence_pmoC[179:198] + gap*20 + sequence_pmoC[218:235] + gap*21
#print(pmoC_tm_seq[218])
print(header)
print(pmoC_tm_seq_3RFR)
header = '>3CHX_3|Chains C, H[auth G], M[auth K]|PmoC|Methylosinus trichosporium (426)'
print(header)
sequence_pmoC = 'MSVTTETTAGAAAGSDAIVDLRGMWVGVAGLNIFYLIVRIYEQIYGWRAGLDSFAPEFQTYWLSILWTEIPLELVSGLALAGWLWKTRDRNVDAVAPREELRRHVVLVEWLVVYAVAIYWGASFFTEQDGTWHMTVIRDTDFTPSHIIEFYMSYPIYSIMAVGAFFYAKTRIPYFAHGFSLAFLIVAIGPFMIIPNVGLNEWGHTFWFMEELFVAPLHWGFVFFGWMALGVFGVVLQILMGVKRLIGKDCVAALVG'
pmoC_tm_seq_3CHX= gap*22 + sequence_pmoC[22:43] + gap*17 + sequence_pmoC[60:82] + gap*27+sequence_pmoC[109:130] +  gap*13 +sequence_pmoC[143:164] + gap*13 + sequence_pmoC[174:197] + gap*22 + sequence_pmoC[217:239] + gap*17
print(pmoC_tm_seq_3CHX)
import re
import sys

hit_file = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/pMOC_curated_alignedPF04896.sthk'
alignment={}
length=172
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

lengthDict={}

for header in alignment.keys():
    lengthDict[len(alignment[header])]=1
print(lengthDict)


if len(lengthDict) > 1:
    sys.exit("ERROR! Aligned sequences differ by length")
for header in alignment.keys():
     if (len(alignment[header].replace('-', ''))) > length:
         if 'CHX' in header:
             print(header)
             print(alignment[header])
             print(alignment[header][18:39])
             print(alignment[header][56:78])
             print(alignment[header][105:126])
             print(alignment[header][139:160])
             print(alignment[header][170:193])
             print(alignment[header][213:235])
         if '3RFR' in header:
             print(header)
             print(alignment[header])
             print(alignment[header][16:42])
             print(alignment[header][53:80])
             print(alignment[header][103:129])
             print(alignment[header][136:163])
             print(alignment[header][175:194])
             print(alignment[header][214:231])
TM = []
TM.extend(list(range(18,39)))
TM.extend(list(range(56,78)))
TM.extend(list(range(105,126)))
TM.extend(list(range(136,160)))
TM.extend(list(range(175,193)))
TM.extend(list(range(214,231)))
print(TM)

hit_file = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/protein-matching-moC_PF04896_alginedPF04896.sthk'
outfile='/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/protein-matching-moC_PF04896_alginedPF04896_TM_alignments_only.fasta'
alignment={}

length=172

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

lengthDict={}

for header in alignment.keys():
    lengthDict[len(alignment[header])]=1
print(lengthDict)



if len(lengthDict) > 1:
    sys.exit("ERROR! Aligned sequences differ by length")

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
