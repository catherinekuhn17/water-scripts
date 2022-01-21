import pkg_resources  
import os.path
import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from qfit.structure import Structure
import glob
from multiprocessing import Pool

def parse_args():
    p = ArgumentParser(description=__doc__)
   # p.add_argument("structure", type=str,
   #                help="Path to PDB files with structures.")
    p.add_argument("--base")
    p.add_argument("--dist", type=float, default='3.0',  help="Distance between water residue and close residues")  
    p.add_argument("--pdb", help="Name of the input PDB.")
    args = p.parse_args()
    return args

# finds for ALL waters (full and partial)
def residues_close(fn, dist, pdb_n):
    print(fn)
    try:
        structure = Structure.fromfile(fn).reorder()
    except NotImplementedError:
        return
    resi_code_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 
                  'HIS','ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                  'THR', 'TRP', 'TYR', 'VAL']
    pdb = structure.extract('resn', 'HOH', '!=') 
    pdb = pdb.extract('e', 'H', '!=') # want no hydrogens
    water = structure.extract('resn', 'HOH', '==')
    water = water.extract('e', 'O', '==') # want no hydrogens
    wat_id, wat_b, wat_q, wat_chain, wat_altloc, pdb_chain, pdb_resi, pdb_atom, pdb_altloc, pdb_b, pdb_resn, pdb_dist, pdb_q,pdb_all_res_b = ([] for i in range(14))
    for chain in np.unique(water.chain): # itterate through each water chain
        tmp_water = water.extract('chain', chain, '==') 
        for i in tmp_water.resi: # itterate through each resi
            for alt in tmp_water.extract('resi', i, '==').altloc: # for each altloc in water
                if alt == '':
                    tmp_water_ext = tmp_water.extract(f'resi {i}') # extract this water
                else:
                    tmp_water_ext = tmp_water.extract(f'altloc {alt} and resi {i}')
                dist_pdb = np.linalg.norm(pdb.coor - tmp_water_ext.coor, axis=1) # find distance of water to protein
                closest_at = np.where(dist_pdb < dist)[0] # find where distance is < 5
                for r, n, c ,a, rn in zip(pdb.resi[closest_at], pdb.name[closest_at], pdb.chain[closest_at], pdb.altloc[closest_at], pdb.resn[closest_at]): # for each place where distance is less than 5
                    if rn in resi_code_list:
                        if a !='': 
                            pdb_ext = pdb.extract(f'name {n} and altloc {a} and resi {r} and chain {c}')
                            pdb_ext_res = pdb.extract(f'altloc {a} and resi {r} and chain {c}')
                        else:
                            pdb_ext = pdb.extract(f'name {n} and resi {r} and chain {c}')
                            pdb_ext_res = pdb.extract(f'resi {r} and chain {c}')
                        disty = np.linalg.norm(pdb_ext.coor - tmp_water_ext.coor, axis=1) # find distance for this (should be min)
                        pdb_dist.append(disty[0]) # append this distance
                        wat_id.append(i) # keep track of water ID
                        wat_altloc.append(tmp_water_ext.altloc) # keep track of water altloc
                        wat_b.append(tmp_water_ext.b[0]) # keep track of water b factor
                        wat_q.append(tmp_water_ext.q[0]) # keep track of water occupancy
                        wat_chain.append(chain)                 
                        pdb_chain.append(c)
                        pdb_resi.append(r)
                        pdb_atom.append(n)
                        pdb_altloc.append(a)
                        pdb_b.append(pdb_ext.b[0])
                        pdb_q.append(pdb_ext.q[0])
                        pdb_resn.append(pdb_ext.resn[0])
                        pdb_all_res_b.append(pdb_ext_res.b[0]) # we want the b of the whole res too

    df=pd.DataFrame()
    df=pd.DataFrame(list(zip(wat_id,wat_altloc, wat_chain, wat_b, wat_q, pdb_chain, pdb_resi, pdb_resn, 
                             pdb_atom, pdb_altloc, pdb_b, pdb_q, pdb_dist, pdb_all_res_b)),
                    columns=['wat_id', 'wat_altloc', 'wat_chain', 'wat_b', 'wat_q', 'prot_chain', 'prot_resi', 'prot_resn', 
                              'prot_atom', 'prot_altloc', 'prot_b', 'prot_q', 'prot_dist', 'prot_all_res_b'])
    df.to_csv('test/' + str(pdb_n) + '_' + str(dist) + '_test3.csv')


def main():
    args = parse_args()
    os.chdir(args.base)
    fns = []
    print(args.pdb)
    for file in glob.glob('*/*.pdb'): # get all pdb files in directory
        fns.append(file)
    print(fns)
    dist = np.ones(len(fns))*args.dist # create array of distances
    args_for_pool = list(zip(fns, dist, [f[0:4] for f in fns])) # setup args for pool
    print(args_for_pool)
    pool = Pool(processes=32)
    results = pool.starmap(residues_close, args_for_pool)

if __name__ == '__main__':
    main()
