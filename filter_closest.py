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
    p.add_argument("--base", help="base directory of PDB files")
    p.add_argument("--filenames", type=str, help="pattern for input dataframes")
    p.add_argument("--filter_type", type=str, default=None,help="how to filter unwanted distances")
    p.add_argument("--cutoff",type=float,default='2.5',  help="distance to start filtering at") 
    args = p.parse_args()
    return args

def find_closest(input_id, input_df, cutoff,filter_type=None):
    '''
    Function for finding the closest atom to a water.
    
    Parameters
    ----------
    filter type : how you want to filter out the "bad" data. Options are match_alt, remvove_alt, and None 
    cutoff : what distance under should be filtered
    
    returns 
    -------

    '''
    print(f'Starting {input_id}')
    
    wat_id, min_dist, wat_b, pdb_b, pdb_q, pdb_chain, pdb_resi, pdb_atom, pdb_resn, pdb_id, wat_q, prot_b, wat_altloc, pdb_altloc, wat_chain, pdb_name, pdb_in_d = ([] for i in range(17))
    
    for w_id in np.unique(input_df['wat_id']):
        
        tmp_df_a = input_df.copy()[input_df['wat_id'] == w_id]
        
        for alt in np.unique(tmp_df_a['wat_altloc']):
            tmp_df = tmp_df_a[tmp_df_a['wat_altloc']==alt]
            wat_id.append(w_id)
            min_d = min(tmp_df['prot_dist'])
            pdb_in_d.append(min_d)
            where = np.where(tmp_df['prot_dist'] == min_d)[0][0]
            w_q = list(tmp_df['wat_q'])[where]
            p_q = list(tmp_df['prot_q'])[where]
            new_df = tmp_df.copy()
            
            if w_q < 1 and p_q < 1:
                while min_d < cutoff:
                    if filter_type == 'remvove_alt':
                        where = np.where(new_df['prot_dist'] == min_d)[0][0]
                        if list(new_df['wat_q'])[where] < 1 and list(new_df['prot_q'])[where] < 1:
                            pdb_c = list(new_df['prot_chain'])[where]
                            pdb_r = list(new_df['prot_resi'])[where]
                            pdb_a = list(new_df['prot_altloc'])[where]
                            index_names = new_df[((new_df['prot_chain']==pdb_c)
                                 &(new_df['prot_resi']==pdb_r)
                                 &(new_df['prot_altloc']==pdb_a))].index
                            if len(index_names) == len(new_df['prot_dist']):
                                break
                            else:
                                new_df = new_df.drop(index_names)
                                min_d = min(new_df['prot_dist'])
                        else:
                            break
                    elif filter_type == 'match_alt':
                        pdb_a = list(tmp_df['wat_altloc'])[where]
                        new_df = tmp_df[(tmp_df['prot_altloc']==pdb_a)]
                        if len(new_df['prot_dist'])>0:
                            min_d = min(new_df['prot_dist'])
                    elif filter_type == None:
                        break                              
            
            min_dist.append(min_d)
            where = np.where(new_df['prot_dist'] == min_d)[0][0]
            pdb_name.append(input_id)
            wat_b.append(list(new_df['wat_b'])[where])
            wat_q.append(list(new_df['wat_q'])[where])
            wat_chain.append(list(new_df['wat_chain'])[where])
            wat_altloc.append(list(new_df['wat_altloc'])[where])
            pdb_b.append(list(new_df['prot_b'])[where])
            pdb_q.append(list(new_df['prot_q'])[where])
            pdb_chain.append(list(new_df['prot_chain'])[where])
            pdb_resi.append(list(new_df['prot_resi'])[where])
            pdb_resn.append(list(new_df['prot_resn'])[where])
            pdb_atom.append(list(new_df['prot_atom'])[where])
            pdb_altloc.append(list(new_df['prot_altloc'])[where])
            prot_b.append(list(new_df['prot_all_res_b'])[where])
   
    df = pd.DataFrame()
    df = pd.DataFrame(list(zip(wat_id, pdb_name, wat_altloc, wat_chain, wat_b, wat_q, pdb_chain, pdb_resi, pdb_resn, 
                             pdb_atom, pdb_altloc, pdb_b, pdb_q, pdb_dist, prot_b, pdb_in_d)),
                      columns=['wat_id', 'wat_altloc', 'wat_chain', 'wat_b', 'wat_q', 'pdb_chain', 'pdb_resi', 'pdb_resn', 
                               'psb_atom', 'pdb_altloc', 'pdb_b', 'ppdb_q', 'prot_dist', 'prot_b', 'pdb_in_d'])
    print(f'Finishing {input_id}')
    return df


def main():
    
    args = parse_args()
    dicty_full=dict()
    
    for file in glob.glob(args.filenames):
        dicty_full[file[0:4]] = pd.read_csv(file)
    args_for_pool = list(zip(dicty_full.keys(), dicty_full.values(), [args.cutoff for i in range(len(dicty_full.keys()))], [args.filter_type for i in range(len(dicty_full.keys()))])) # setup args for pool

    pool = Pool(processes=32)
    results = pool.starmap(find_closest, args_for_pool)

if __name__ == '__main__':
    main()




