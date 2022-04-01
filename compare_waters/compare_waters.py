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
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import matplotlib.pyplot as plt
import WATER_LOCS as WATER_LOCS
import WATER_LOCS_SS as WATER_LOCS_SS

# simple functions to obtrain the protein part of the structure and water part of the stricture
def read_protein(structure):
    pdb = structure.extract('resn', 'HOH', '!=') 
    pdb = pdb.extract('e', 'H', '!=')
    return pdb

def read_waters(structure):
    water = structure.extract('resn', 'HOH', '==')
    water = water.extract('e', 'O', '==') 
    return water

def get_chi1(res, protein):
    '''
    Function for finding X1 of a residue
    
    Parameters:
    ----------
    res : the residue you want to find the X1 of
    protein : the structure file of the whole protein (for cases of altlocs)
    
    Returns:
    --------
    the x1 symbol (g+ for 0-120, t for 120-240, g- for 240-360)
    
    '''
    
    chi1_dict = {
        0 : 'gplus',
        1 : 't',
        2 : 'gminus'
    }
    if res.type == 'rotamer-residue':
        if res.nchi > 0:
            try:   
                chi_val = res.get_chi(1)
            except IndexError: # for cases where the altloc doesnt have the residues that arent in an altloc
                try: # have to extract from rest or residue not included in altloc
                    chi_val = list(protein.extract(f'resi {res.resi[0]} and chain {res.chain[0]}').residues)[0].get_chi(1) 
                except IndexError: # if that doesn't work just exit 
                    return None
            if chi_val<0:
                chi_val+=360
            chi_val = chi1_dict[int(chi_val//120)]
        else:
            chi_val = 'NA'
        return chi_val

def get_dssp(fn):
    '''
    Function for finding SS of each residue of a protein (and also the ASA_
    
    Parameters:
    ----------
    fn : filename of pdb
    
    Returns:
    --------
    dict_dssp : a dictionary of each residues SS, formatted {chain : {resi : ss}}
    dict_asa : a dictionary of each residues ASA, formatted {chain : {resi : asa}}
    
    '''
    ss_converter = {
        'H' : 'H',
        'G' : 'H',
        'E' : 'E',
        'B' : 'E',
        'I' : 'C',
        'T' : 'C',
        'S' : 'C',
        '-' : 'NA'
    }

    p = PDBParser()
    structure = p.get_structure('id', fn)
    model = structure[0]
    dssp = DSSP(model, fn)
    
    chains = [list(dssp.keys())[i][0] for i in range(len(list(dssp.keys())))]
    resis = [list(dssp.keys())[i][1][1] for i in range(len(list(dssp.keys())))]
    ss = [ss_converter[dssp[list(dssp.keys())[i]][2]] for i in range(len(list(dssp.keys())))]
    asa = [dssp[list(dssp.keys())[i]][3] for i in range(len(list(dssp.keys())))]
    
    dict_dssp={}
    dict_asa = {}
    for c, r, s, a in zip(chains, resis, ss, asa):
        if c in dict_dssp.keys():
            dict_dssp[c][r] = s
            dict_asa[c][r] = a
        else:
            dict_dssp[c] = {}
            dict_asa[c] = {}
            dict_dssp[c][r] = s
            dict_asa[c][r] = a

    return dict_dssp, dict_asa

def find_closest_waters_id(protein, water, res, dist_cutoff):
    '''
    function for finding the closest waters to a specified residue
    
    Parameters
    ----------
    protein : the protein protion of the pdb (structure)
    water : the water portion of the pdb (structure)
    res : the residue you are searching for waters around
    dist_cutoff : the furthest from the protein you will search for waters
    
    Returns
    --------
    the indexes of the waters closest to the residue
    
    '''
    water_dict={}
    water_id_all = []
    for name, rc  in zip(res.name, res.coor):
        if res.altloc[0]!='': # dealing with altlocs, want to make sure we JUST have the waters with the same altloc
            water = water.extract(f'altloc {res.altloc[0]}')
        dist = [np.linalg.norm(water.coor - rc, axis=1)][0]
        idx_upper = np.where(dist < dist_cutoff)[0] # filter to waters within cutoff dist of res
        if name in pos_c:
            idx_lower = np.where(dist > 3.1)[0] # if dist to carbon greater than 3.1, consider it a class and ignore
        else:
            idx_lower = np.where(dist > 2.5)[0] # if dist to N,O, or S greater than 2.5, consider it a class and ignore
        water_id = list(set(idx_upper) & set(idx_lower))
        water_id = idx_upper
        water_id_all = np.append(water_id_all, water_id)
    wat_id_fin = np.unique(np.array(water_id_all).flatten())
    watnum = [water.resi[int(i)] for i in wat_id_fin] 
    wat_alt = [water.altloc[int(i)] for i in wat_id_fin]
    wat_chain = [water.chain[int(i)] for i in wat_id_fin]
    
    return watnum, wat_alt, wat_chain # save the water resi #, altloc, and chain

def trilaterate(P1,P2,P3,P4, r1,r2,r3,r4):
    '''
    function to do a trilaterate calculation
    
    Parameters
    ----------
    P1, P2, P3, P4 : coordinates of 4 atoms in structure 
    r1, r2, r3, r4 : distances of water to these 4 atoms
    
    Returns 
    -------
    xyz coordinate of this water in relation to the input residue
    '''
    temp1 = P2.flatten()- P1.flatten()
    e_x = temp1/np.linalg.norm(temp1)                              
    temp2 = P3.flatten() - P1.flatten()                                      
    i = np.dot(e_x,temp2)                                   
    temp3 = temp2 - i*e_x                               
    e_y = temp3/np.linalg.norm(temp3)
    e_z = np.cross(e_x,e_y)                                 
    d = np.linalg.norm(P2.flatten()-P1.flatten())                                      
    j = np.dot(e_y,temp2)                                   
    x = (r1*r1 - r2*r2 + d*d) / (2*d)                    
    y = (r1*r1 - r3*r3 -2*i*x + i*i + j*j) / (2*j)       
    temp4 = r1*r1 - x*x - y*y                            
    z = np.sqrt(max(0,temp4))  # look into this?                               
    p_12_a = P1.flatten() + x*e_x + y*e_y + z*e_z                  
    p_12_b = P1.flatten() + x*e_x + y*e_y - z*e_z 
    
    dist1=np.linalg.norm(P4-p_12_a)
    dist2=np.linalg.norm(P4-p_12_b)
    if np.abs(r4-dist1)<np.abs(r4-dist2):
        return p_12_a.reshape(3, 1).T[0]
    else: 
        return p_12_b.reshape(3, 1).T[0]
    #return p_12_a.reshape(3, 1).T[0], p_12_b.reshape(3, 1).T[0]
    
    
def rmsd_waters(pdb_id, protein, water, resname, x1, ss, water_locs, dist_cutoff, dssp_dict):
    '''
    function for finding the rmsd of watAA placed waters and actually placed waters, within a cutoff
    
    Parameters:
    ----------
    protein : protein structure object
    resname : residue name to look at
    x1 : x1 angle to look at
    ss : secondary structure (None if you do not care about ss)
    water_loc : dictionary of water locations from watAA
    dist_cutoff : furthest waters to pull from 
    dssp_dict : seconday structure dictionary (can be None)
    
    Returns
    -------
    r_rmsd : a dictionary formatted
            {(chain ID, residue ID, residue altloc) : 
              {water ID : (dist to closest watAA, water_q, watAA water number, watAA xyz coord}}
    '''
    res = protein.extract(f'resn {resname}')
     # dictionary for location of waters (will separate close and far waters
    r_rmsd_full = {}
    r_rmsd_part = {}
    if ss!=None:
        close_atoms = water_locs[resname][ss][x1] # if we DO care about SS
        res_list = list(res.residues)
        res_list_new=[]
        for r in res_list:
            if r.chain[0] in dssp_dict.keys():
                if r.resi[0] in dssp_dict[r.chain[0]].keys():
                    if dssp_dict[r.chain[0]][r.resi[0]] == ss: # only look at res with this SS
                        res_list_new.append(r)
                    
        res_list = res_list_new  
    else:
        close_atoms = water_locs[resname][x1] # if we don't care about SS
        res_list = list(res.residues)
    r_rmsd_full = {}
    r_rmsd_part = {}
    for r in res_list: # iterating through each residue with this SS, X1, and resname combo:
        if get_chi1(r, protein) == x1: # if x1 matches what we want
            f, p = find_waters(pdb_id, protein, water, r, x1, ss, water_locs, dist_cutoff, dssp_dict, close_atoms)
            r_rmsd_full[str((r.chain[0], r.resi[0], r.altloc[0], r.resn[0]))] = f
            r_rmsd_part[str((r.chain[0], r.resi[0], r.altloc[0], r.resn[0]))] = p
    return pdb_id, r_rmsd_full, r_rmsd_part

def find_waters(pdb_id, protein, water, r, x1, ss, water_locs, dist_cutoff, dssp_dict, close_atoms):
    ''' 
    function for finding the closert watAA placed water to the actually placed waters.
    
    Parameters:
    ----------
    protein : protein structure object
    resname : residue name to look at
    x1 : x1 angle to look at
    ss : secondary structure (None if you do not care about ss)
    water_loc : dictionary of water locations from watAA
    dist_cutoff : furthest waters to pull from 
    dssp_dict : seconday structure dictionary (can be None)
    close_atoms : the sets of atom distances for watAA waters
    
    
    Returns
    -------
    r_rmsd : a dictionary formatted
            {(chain ID, residue ID, residue altloc) : 
              {water ID : (dist to closest watAA, water_q, watAA water number, watAA xyz coord}}
    '''
    wat_loc={}
    w_rmsd_partial = {}
    w_rmsd_full = {}
    for idx in range(0, len(close_atoms)): # going through sets of the closest atoms
        
        atom = list(close_atoms[idx+1].keys())
        dist = list(close_atoms[idx+1].values()) # these are the distances that will be used in trilaterate
        
        P={} # this is the points that will be used in trilaterate
        for i,at in enumerate(atom): 
            try: 
                P[i] = r.extract(f'chain {r.chain[0]} and name {at}').coor[0] 
                # has to do with how altlocs sometimes mess things up :/ 
            except IndexError:
                if sum(r.altloc == '') == len(r.altloc): # if it's just the non-altloc atoms of the altloc
                    return w_rmsd_partial, w_rmsd_full
                if at not in protein.extract(f'chain {r.chain[0]} and name {at} and resi {r.resi[0]}').name:
                    return w_rmsd_partial, w_rmsd_full
                else:
                    P[i] = protein.extract(f'chain {r.chain[0]} and name {at} and resi {r.resi[0]}').coor[0]
        
        if len(list(P.values())) == 4: # only if we were able to extract all P's we need
            wat_loc[idx+1] = trilaterate(P[0], P[1], P[2], P[3],dist[0], dist[1], dist[2], dist[3])
    wat_idx, wat_alt, wat_chain = find_closest_waters_id(protein, water, r, dist_cutoff) # id's of waters closest to residue
    for w, alt, c in zip(wat_idx, wat_alt, wat_chain): # going through each water
        rmsd=10 # just an initial val
        dist = 10
        if alt!='':
            wat = water.extract(f'resi {w} and altloc {alt} and chain {c}') # extract waters
        else:
            wat = water.extract(f'resi {w} and chain {c}')
        #print(list(wat_loc.values()))
        water_coor = wat.coor[0]
        for idx, coord in wat_loc.items():
            rmsd_tmp =  np.sqrt(((water_coor[0]-coord[0])**2+(water_coor[1]-coord[1])**2+(water_coor[2]-coord[2])**2)/3)
            if rmsd_tmp < rmsd:
                rmsd = rmsd_tmp # calculate rmsd to watAA water placements
                idx_fin = idx # index of final water used

        for idx, coord in wat_loc.items():
            dist_tmp = np.linalg.norm(water_coor - coord)
            if dist_tmp < dist:
                dist = dist_tmp # calculate dist to watAA water placements
               # idx_fin = idx # index of final water used

        if wat.q[0] < 1:
            w_rmsd_partial[wat.resi[0]] = (rmsd, dist, idx_fin, list(wat_loc.values())) # output with format 
        else:
            w_rmsd_full[wat.resi[0]] = (rmsd, dist, idx_fin, list(wat_loc.values()))      

    return w_rmsd_full, w_rmsd_partial
 
def main:
        
