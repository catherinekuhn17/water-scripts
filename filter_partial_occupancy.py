from Bio.PDB import *
import os
import shutil
base_folder = '/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/pdb_files/download2/'
pdb_ids_fn = '/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/pdb_files/download2/pdb_ids_all.txt'
outfile = '/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/pdb_files/download2/pdb_id_list_part_occupancy.txt'
with open(pdb_ids_fn) as infile, open(outfile, 'w') as outfile:
    while True:
        pdb_id = infile.readline().rstrip()
        pdb_fn = os.path.join(base_folder, 'all',pdb_id, pdb_id + '.pdb')
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, pdb_fn)

        atoms = structure.get_atoms()
        for a in atoms:
            if 'HOH' in str(a.get_parent()):
                if a.get_occupancy()<1:
<<<<<<< HEAD
                    partial_occupancy = True
        if partial_occupancy:
            outfile.write(pdb_id+"\n")
            os.mkdir(os.path.join(base_folder, 'partial_occupancy',pdb_id))
            file_names = os.listdir(os.path.join(base_folder, 'all',pdb_id))
            for file_name in file_names:
                shutil.copy(os.path.join(base_folder, 'all', pdb_id, file_name), os.path.join(base_folder,'partial_occupancy', pdb_id, file_name))
=======
                    outfile.write(pdb_id+"\n")
                    os.mkdir(os.path.join('partial_occupancy',pdb_id))
                    file_names = os.listdir(os.path.join('all',pdb_id))
                    for file_name in file_names:
                        shutil.move(os.path.join('all', pdb_id, file_name), os.path.join('partial_occupancy', pdb_id, file_name))
>>>>>>> 56507b3ae89b2e5ffc0748847aa7e42e6071d903
        



