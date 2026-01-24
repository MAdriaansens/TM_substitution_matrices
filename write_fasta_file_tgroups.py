from Bio import SeqIO

writing_site_tm = '/nesi/nobackup/uc04105/PDB_alpha/tma_topdb/tgroup_sep'
writing_site_fl = '/nesi/nobackup/uc04105/PDB_alpha/fl_topdb/tgroup_sep'

with open("/nesi/nobackup/uc04105/PDB_alpha/Alpha-Helical_tm_PDBTM_ECOD_family.tsv", 'r') as individual_protein_data:
    for protein in individual_protein_data:
        Tgroup = protein.split('\t')[2].replace(' ', '_').replace('/', '_slash_')
        tm_seq = protein.split('\t')[7]
        fl_seq = protein.split('\t')[9]
        pdb_id = protein.split('\t')[0]
        with open('{}/{}_tm.fasta'.format(writing_site_tm, Tgroup), 'a') as tm_in:
            line = '>' + pdb_id + '_tm_Tgroup:{}'.format(Tgroup) + '\n' + tm_seq + '\n'
        tm_in.close()
        with open('{}/{}_fl.fasta'.format(writing_site_fl, Tgroup), 'a') as fl_in:
            line = '>' + pdb_id + '_fl_Tgroup:{}'.format(Tgroup) + '\n' + fl_seq + '\n'
        fl_in.close()
    
