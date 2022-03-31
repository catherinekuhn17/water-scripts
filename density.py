import os
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
from sklearn.cluster import MeanShift
from sklearn.cluster import  estimate_bandwidth
from scipy.spatial.distance import cdist

def parse_args():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--out", help="directory to read from/write files out to")
    p.add_argument("--pt", help="percentile cuttoff for meanshift")
    p.add_argument("--length", help="max number of waters to use")
    args = p.parse_args()
    return args

def find_density(coords_all, pt, length):
    '''
    Function for finding density of water distribution around a set of 4 atoms, and
    also clustering density using Mean Shift clustering from sklearn
    Parameters
    ----------
    coords_all : dictionary
        a kind of confusing dictionary (sorry) formatted
        {4 atom set : tuple(center 4 at coord,
                            4 at coords, 
                            water coords, 
                            normalized b vals for water, 
                            occupancy for water)}
    pt : int
        percentile cutoff for what to do mean shift clustering on
    length : int
        maximum size of data to do analysis on
    Returns
    -------
    writes out to a bunch of npy files including:
        center_coor : center coordinates of clusters from mean shift
        spread : spread around each center coordinate
        dens_all : density value of every point
        coord_set_all : coordinate value of every point
        rel_b_list_all : normalized b factor of every point
        q_list_all : occupancy of every point
        labs : labels of points included in meanshift clustering
        cutoff_idx : indices of points included in meanshift clustering
    
    '''
    print('starting')
    # initialize dictionaries
    center_coor={} 
    spread={}
    dens_all={}
    coord_set_all={}
    rel_b_list_all = {}
    q_list_all={}
    labs={}
    cutoff_idx={}
    cutoff_idx_current={}
    for atom_set, atom_set_v in coords_all.items():
        center_coor[atom_set] = {}
        spread[atom_set] = {}
        dens_all[atom_set] = {}
        coord_set_all[atom_set] = {}
        rel_b_list_all[atom_set] = {}
        q_list_all[atom_set] = {}
        labs[atom_set]={}
        cutoff_idx[atom_set]={}
        cutoff_idx_current[atom_set]={}
        for dih_assig, coords in atom_set_v.items():
            center_coor_tmp=[]
            spread_tmp=[]
            # water coordinates
            wat_coords = np.array(list(coords_all[atom_set][dih_assig][2])).flatten().reshape(-1,3)
            # normalized b for waters
            rel_b_val = np.array(coords_all[atom_set][dih_assig][3])
            # occupancy
            q_val = np.array(coords_all[atom_set][dih_assig][4])
            # this is so we can subsample a bit
            sampling = max(1, int(len(wat_coords)/length))
            wat_coords = wat_coords[::sampling]
            rel_b_val = rel_b_val[::sampling]
            q_val = q_val[::sampling]
            # first, we fo KDE, to find a density value for each point
            # here I have bandwidth = 1 since that generally works and is the default param
            kde = KernelDensity(kernel='gaussian', bandwidth=1, rtol=1E-4, atol=1E-4).fit(wat_coords)
            density = kde.score_samples(wat_coords)
            # finding indices of points with density above a specified percentile
            idx = np.where(density>np.percentile(density, pt))[0]
            cutoff_idx[atom_set][dih_assig]={}
            # just so that we can see the indices of a bunch of percentile cutoffs
            for pt_i in np.arange(0, 100, 10):
                cutoff_idx[atom_set][dih_assig][pt_i] = np.where(density>np.percentile(density, pt_i))[0]
            if len(idx)>0:
                print(atom_set, dih_assig)
                # now doing meanshift on these points
                # the bin_seeding means that points are binned into grids. Increasing min_bin_freq
                # increases the # of points that need to be in a bin
                # cluster_all=False means we don't cluster everything
                # these params have not been optimized at all
                msc = MeanShift(bandwidth=1, bin_seeding=False)
                msc.fit(wat_coords[idx])
                cluster_centers = msc.cluster_centers_
                labels = msc.labels_
                cluster_label = np.unique(labels)
                n_clusters = len(cluster_label)
                # this is to find "spread" of points (or distance of furthest point within a
                # cluster to it's cluster center)
                print(n_clusters)
                print(cluster_label)
                for i, cc in enumerate(cluster_centers):
                    pos = np.argsort(cdist([cc], wat_coords[idx]))[0][0]
                    center_coor_tmp.append(wat_coords[idx][pos])
                    radius = max(cdist(wat_coords[idx][labels==i], [cluster_centers[i]]))
                    print(radius)
                    spread_tmp.append(radius)
                center_coor[atom_set][dih_assig] = center_coor_tmp # xyz of center of cluster
                spread[atom_set][dih_assig] = spread_tmp # "radius" of cluster
                dens_all[atom_set][dih_assig] = density # density val of all points
                coord_set_all[atom_set][dih_assig] = wat_coords # xyz of all points
                rel_b_list_all[atom_set][dih_assig] = rel_b_val # b factors of points
                q_list_all[atom_set][dih_assig] = q_val # occupancy of waters
                labs[atom_set][dih_assig] = labels # cluster each point belongs in
                cutoff_idx_current[atom_set][dih_assig] = idx # indices above percentile density
                
    np.save(f'center_coor_{length}_{pt}.npy', center_coor) 
    np.save(f'spread_{length}_{pt}.npy', spread) 
    np.save(f'density_vals_{length}_{pt}.npy', dens_all) 
    np.save(f'all_xyz_coords_{length}_{pt}.npy', coord_set_all) 
    np.save(f'rel_b_list_{length}_{pt}.npy', rel_b_list_all) 
    np.save(f'q_list_{length}_{pt}.npy', q_list_all) 
    np.save(f'labels_{length}_{pt}.npy', labs)
    np.save(f'cutoff_idx_{length}_{pt}.npy', cutoff_idx_current)
    np.save(f'cutoff_idx_{length}_all.npy', cutoff_idx)
    return

def main():
    args = parse_args()
    out_dir = args.out
    pt = args.pt
    length = args.length
    os.chdir(out_dir)
    coords_all = np.load('dih_info.npy',allow_pickle='TRUE').item()
    find_density(coords_all, pt, length)

if __name__ == '__main__':
    main()