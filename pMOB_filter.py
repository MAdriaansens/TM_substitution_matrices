
header = '>3RFR_1|Chains A, J[auth E], K[auth I]|PmoB|Methylocystis sp. M (51782)'
sequence_pmoB = 'MKKLVKLAAFGAAAAVAATLGAIAPASAHGEKSQQAFLRMRTLNWYDVQWSKTTVNVNEEMILSGKVHVFSAWPQAVANPRVSFLNAGEPGPVLVRTAQFIGEQFAPRSVSLEIGKDYAFSINLRGRRAGRWHVHAQINVEGGGPIIGPGQWIEIKGDMKDFTDPVTLLDGSTVDLENYGISRIYAWHLPWLAVGAAWILFWFIRKGIIASYVRVAEGRPDDVIGDDDRRIGAIVLALTILATIVGYAVTNSTFPRTIPLQAGLQKPLTPIETEGTVGVGKEQVTTELNGGVYKVPGRELTINVKVKNGTSQPVRLGEYTAAGLRFLNPTVFTQKPDFPDYLLADRGLSNDDVIAPGESKEIVVKIQDARWDIERLSDLAYDTDSQVGGLLFFFTPDGKRFAAEIGGPVIPKFVAGDMP'
print(len(sequence_pmoB))

#3RFR says it is TM 185-204 / 231-249

gap = '-'

Sequence_pmoB = gap*180+ sequence_pmoB[180:205] + gap*18 + sequence_pmoB[223:250] + 169*gap

pmoB_tm_3RFR = '>3RFR_1|Chains A, J[auth E], K[auth I]|PmoB|Methylocystis sp. M (51782)_TMelements' + '\n' + Sequence_pmoB + '\n'
print(pmoB_tm_3RFR)

header ='>3CHX_1|Chains A, F[auth E], K[auth I]|PmoB|Methylosinus trichosporium (426)'

sequence_pmoB = 'HGEKSQQAFLRMRTLNWYDVKWSKTSLNVNESMVLSGKVHVFSAWPQAVANPKSSFLNAGEPGPVLVRTAQFIGEQFAPRSVSLEVGKDYAFSIDLKARRAGRWHVHAQINVEGGGPIIGPGQWIEIKGDMADFKDPVTLLDGTTVDLETYGIDRIYAWHFPWMIAAAAWILYWFFKKGIIASYLRISEGKDEEQIGDDDRRVGAIVLAVTILATIIGYAVTNSTFPRTIPLQAGLQKPLTPIIEEGTAGVGPHVVTAELKGGVYKVPGRELTIQVKVTNKTDEPLKLGEYTAAGLRFLNPDVFTTKPEFPDYLLADRGLSTDPTPLAPGETKTIEVKVQDARWDIERLSDLAYDTDSQIGGLLMFFSPSGKRYATEIGGPVIPKFVAGDMP'
pmoB_tm_3CHX = '{}'.format(header) + '\n' + gap*153 + gap*39+ sequence_pmoB[153:174] + gap*26 + sequence_pmoB[200:222] + gap*170 + '\n'
print(pmoB_tm_3CHX)
import re
import sys

hit_file = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/pMOB_curated_alginedPF04744.sthk'
alignment={}
length=266
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
             print(alignment[header][152:173])
             print(alignment[header][199:221])
         if '3RFR' in header:
             print(header)
             print(alignment[header][152:173])
             print(alignment[header][194:221])
TM = []
TM.extend(list(range(152,173)))
TM.extend(list(range(194,221)))

#consensus

hit_file = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/protein-matching-moB_PF04744_alginedPF04744.sthk'
outfile='/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/protein-matching-moB_PF04744_alginedPF04744_TM_alignments_only.fasta'
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
