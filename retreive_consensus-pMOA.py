import re
import sys

hit_file = '/home/mad149/00_nesi_projects/uc04105_nobackup/PDB_alpha/HMMalign/pMOA_curated_alignedPF02461.sthk'
# define alignment
alignment={}



# read in alignment

with open(hit_file) as hmmalignment:
    for line in hmmalignment:
        if(line[0]!="#"):
            splitLine=line.split()

            if(len(splitLine)==2):
                stringFilter = lambda text: re.sub('[.*]', '', splitLine[1])
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
sequence_pmoA = 'MFTSKSGGAIGPFHSVAEAAGCVKTTDWMFLTLLFLAVLGGYHIHFMLTAGDWDFWVDWKDRRMWPTVVPILGVTFAAAAQAFFWENFKLPFGATFAVSGLLIGEWINRYCNFWGWTYFPISLVFPSALVVPALWLDIIMLLSGSYVITAVVGSLGWGLLFYPNNWPAIAALHQATEQHGQLMSLADLVGFHFVRTSMPEYIRMVERGTLRTFGKEVVPVAAFFSGFVSMMVYFLWWFVGKWYSTTKVIQKI'
gap = '-'


pmoA_tm_seq_3CHX = gap*27+sequence_pmoA[27:53] + gap *11 + sequence_pmoA[64:85] + gap*4 + sequence_pmoA[89:113] + gap*12 + sequence_pmoA[125:168] + gap*32 + sequence_pmoA[199:218] + gap*7 + sequence_pmoA[225:247] + gap*5

print(pmoA_tm_seq_3CHX)

pmoA_tm_seq_3RFR = gap*24 +sequence_pmoA[24:55] + gap*17 + sequence_pmoA[62:115] + gap*5 + sequence_pmoA[120] + gap*3 + sequence_pmoA[123:169]+ gap*12 + sequence_pmoA[181:205] + gap*11 + sequence_pmoA[216:245] + gap*6 + sequence_pmoA[-2:]
print(pmoA_tm_seq_3RFR)

sequence_pmoA_3CHX2 = 'MFTSKSGGAIGPFHSVAEAAGCVKTTDWMFLTLLFLAVLGGYHIHFMLTAGDWDFWVDWKDRRMWPTVVPILGVTFAAAAQAFFWENFKLPFGATFAVSGLLIGEWINRYCNFWGWTYFPISLVFPSALVVPALWLDIIMLLSGSYVITAVVGSLGWGLLFYPNNWPAIAALHQATEQHGQLMSLADLVGFHFVRTSMPEYIRMVERGTLRTFGKEVVPVAAFFSGFVSMMVYFLWWFVGKWYSTTKVIQKI'
sequence_pmoA_3RFR = 'MSQSKSGGAVGPFNSVAEAAGCVATTDWMLLVLLFFAVLGGYHVHFMLTAGDWDFWVDWKDRRMWPTVLPILGVTFCAASQAFWWVNFRLPFGAVFAVLGLMIGEWINRYVNFWGWTYFPISLVFPSAMIVPAIWLDVILLLSGSYVITAVVGSLGWGLLFYPNNWPAIAAFHQATEQHGQLMTLADLIGLHFVRTSMPEYIRMVERGTLRTFGKDVVPVAAFFSGFVSMMVYFLWWFMGRWYSTTKRIEQI'
print(sequence_pmoA_3CHX2[14])
for header in alignment.keys():
    print(header)
    print(len(alignment[header]))
    if 'CHX' in header:
        #visually match to see what parts of the HMM match with TM elements of CHX3 pMOA
        print(alignment[header][213:235])
        print(alignment[header][187:206])
        print(alignment[header][113:154])
        print(alignment[header][77:101])
        print(alignment[header][52:73])
        print(alignment[header][15:41])
    elif 'RFR' in header:
        print(alignment[header][12:43])
        print(alignment[header][50:103])
        print(alignment[header][108])
        print(alignment[header][111:157])
        print(alignment[header][169:193])
        print(alignment[header][204:233])
        print(alignment[header][236])
#the consensus of these these two sequences is:
#1 15-41
#2 52-73
#3 77-101
#4 113-154
#5 187-193
#6 204-206
#7 213-233
