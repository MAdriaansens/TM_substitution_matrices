import sys
A = sys.argv[1]
Res_list = []
with open(A, 'r') as Alignment:
    for sequence in Alignment:
        if sequence[0] != '>':
            length = len(sequence)-1
            break
Alignment.close()
print(length)
#get length of the alignment


def return_list_residues_on_position(number):
    #test is a crude way of removing gaps of looping through a file and only taking the number i of the residue
    residue_list = []
    with open(A, 'r') as Alignment:
        for sequence in Alignment:
            if sequence[0] != '>':
                if sequence[number] != '-':
                    residue = sequence[number]

                    residue_list.append(residue)
    return(residue_list)

from collections import Counter
total_pairs = []
for number in range(length):
    #per aa in residue (on #i position) it will return a list of residues present on it
    #from this list it well then
    residue_list_on_position = return_list_residues_on_position(number)
    pairs = [(i, j) for i in residue_list_on_position for j in residue_list_on_position]
    total_pairs.append(pairs)

total_list = sum(total_pairs, [])
count_dict = dict(Counter(total_list).items())

print(count_dict)
print(len(count_dict))
all_pairs = len(total_list)
#only issue of this chunk is that it does not add DE or ED as one, which is fine but needs to be resolved later on

import numpy as np
results = list(count_dict.items())
Freq_dict = {}
control_list = []
for entry in results:
    First_let = entry[0][0]
    Second_let = entry[0][1]
    #this is for the E -> E and F -> F residues. ie the ones which 'do not change'
    if First_let == Second_let:
        combi = First_let + Second_let
        Freq_dict[combi] = entry[1]
    else:
        combi = First_let + Second_let
        combi = ''.join(sorted(combi))

        if combi in control_list:
            pass
        else:
            #for K -> L and L-> K
            Freq_dict[combi] = entry[1]*2
            control_list.append(combi)
print(Freq_dict)
print(control_list)


A_array = ['AA']
C_array = ['AC', 'CC']
D_array = ['AD', 'CD', 'DD']
E_array = ['AE', 'CE', 'DE', 'EE']
F_array = ['AF', 'CF', 'DF', 'EF', 'FF']
G_array = ['AG', 'CG', 'DG', 'EG', 'FG', 'GG']
H_array = ['AH', 'CH', 'DH', 'EH', 'FH', 'GH', 'HH']
I_array = ['AI', 'CI', 'DI', 'EI', 'FI', 'GI', 'HI', 'II']
K_array = ['AK', 'CK', 'DK', 'EK', 'FK', 'GK', 'HK', 'IK', 'KK']
L_array = ['AL', 'CL', 'DL', 'EL', 'FL', 'GL', 'HL', 'IL', 'KL', 'LL']
M_array = ['AM', 'CM', 'DM', 'EM', 'FM', 'GM', 'HM', 'IM', 'KM', 'LM', 'MM']
N_array = ['AN', 'CN', 'DN', 'EN', 'FN', 'GN', 'HN', 'IN', 'KN', 'LN', 'MN', 'NN']
P_array = ['AP', 'CP', 'DP', 'EP', 'FP', 'GP', 'HP', 'IP', 'KP', 'LP', 'MP', 'NP', 'PP']
Q_array = ['AQ', 'CQ', 'DQ', 'EQ', 'FQ', 'GQ', 'HQ', 'IQ', 'KQ', 'LQ', 'MQ', 'NQ', 'PQ', 'QQ']
R_array = ['AR', 'CR', 'DR', 'ER', 'FR', 'GR', 'HR', 'IR', 'KR', 'LR', 'MR', 'NR', 'PR', 'QR', 'RR']
S_array = ['AS', 'CS', 'DS', 'ES', 'FS', 'GS', 'HS', 'IS', 'KS', 'LS', 'MS', 'NS', 'PS', 'QS', 'RS', 'SS']
T_array = ['AT', 'CT', 'DT', 'ET', 'FT', 'GT', 'HT', 'IT', 'KT', 'LT', 'MT', 'NT', 'PT', 'QT', 'RT', 'ST', 'TT']
V_array = ['AV', 'CV', 'DV', 'EV', 'FV', 'GV', 'HV', 'IV', 'KV', 'LV', 'MV', 'NV', 'PV', 'QV', 'RV', 'SV', 'TV', 'VV']
W_aray =  ['AW', 'CW', 'DW', 'EW', 'FW', 'GW', 'HW', 'IW', 'KW', 'LW', 'MW', 'NW', 'PW', 'QW', 'RW', 'SW', 'TW', 'VW', 'WW']
Y_array = ['AY', 'CY', 'DY', 'EY', 'FY', 'GY', 'HY', 'IY', 'KY', 'LY', 'MY', 'NY', 'PY', 'QY', 'RY', 'SY', 'TY', 'VY', 'WY', 'YY']
array_list = [A_array, C_array, D_array, E_array, F_array, G_array, H_array, I_array, K_array,L_array, M_array, N_array, P_array, Q_array, R_array, S_array, T_array, V_array, W_aray, Y_array]

for array in array_list:
    number= 0
    for pair in array:
        if pair in Freq_dict.keys():
            array[number] = Freq_dict[pair]/all_pairs
        else:
            array[number] = 'NA' #0.1/all_pairs #will become 0.1 later on
        number =  number + 1

order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N','P','Q','R','S','T','V','W','Y']
num_array = 0
for array in array_list:
    print(order[num_array], array)
    num_array = num_array + 1
    if num_array == 20:
        print('   A    C    D    E    F   G    H   I   K   L    M    N   P   Q    R    S   T    V    W    Y')
