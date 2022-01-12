from Bio.PDB import *
import os
import shutil

pdb_ids_fn = 'pdb_id_list.txt'
outfile = 'pdb_id_list_part_occupancy.txt'
with open(pdb_ids_fn) as infile, open(outfile, 'w') as outfile:
    while True:
        pdb_id = infile.readline().rstrip()
        pdb_fn = os.path.join('all',pdb_id, pdb_id + '.pdb')
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, pdb_fn)

        atoms = structure.get_atoms()
        partial_occupancy = False
        for a in atoms:
            if 'HOH' in str(a.get_parent()):
                if a.get_occupancy()<1:
                    partial_occupancy = True
        if partial_occupancy:
            outfile.write(pdb_id+"\n")
            os.mkdir(os.path.join('partial_occupancy',pdb_id))
            file_names = os.listdir(os.path.join('all',pdb_id))
            for file_name in file_names:
                shutil.move(os.path.join('all', pdb_id, file_name), os.path.join('partial_occupancy', pdb_id, file_name))
        



