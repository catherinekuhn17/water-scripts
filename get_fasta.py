import sys
from Bio import SeqIO
import os


pdb_ids_fn = 'pdb_id_list.txt'
outfil = 'all_fasta.fasta'
with open(pdb_ids_fn) as infile, open(outfil, 'w') as outfile:
    while True:
        pdb_id = infile.readline().rstrip()
        pdb_fn = os.path.join('all',pdb_id, pdb_id + '.pdb')
        with open(pdb_fn, 'r') as pdb_file:
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                outfile.write('>' + str(record.id)+"\n")
                outfile.write(str(record.seq)+"\n")
