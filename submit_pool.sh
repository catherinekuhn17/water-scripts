#!/bin/bash
#$ -l h_vmem=2G
#$ -l mem_free=2G
#$ -t 1-1
#$ -l h_rt=10:00:00
#$ -pe smp 1
#$ -R yes
#$ -V

#Stephanie Wankowicz
#Started: 19-09-05
#Last Editted: 19-09-20


#__________________SET PATHS________________________________________________#
source /wynton/group/fraser/swankowicz/phenix-installer-1.19.2-4158-intel-linux-2.6-x86_64-centos6/phenix-1.19.2-4158/phenix_env.sh
export PATH="/wynton/home/fraserlab/swankowicz/anaconda3/bin:$PATH"
source activate qfit

PDB_file=/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/pdb_files/download2/pdb_id_list.txt
base_dir='/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/pdb_files/download2/all2/test'
cd $base_dir
echo $base_dir
#________________________________________________Qfit Analysis________________________________________________#

#python /wynton/home/rotation/ckuhn/Desktop/Fraser_lab/water-scripts/residues_close_pool.py --base=$base_dir --pdb='*/*.pdb' --dist=5.0

python /wynton/home/rotation/ckuhn/Desktop/Fraser_lab/water-scripts/filter_closest.py --base=$base_dir --filenames='*pool.csv' --cutoff=2.5 
