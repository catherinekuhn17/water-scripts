#!/bin/bash
#$ -l h_vmem=500G
#$ -l mem_free=500G
#$ -t 1-1
#$ -l h_rt=10:00:00
#$ -pe smp 1
#$ -R yes
#$ -V
#$ -e /wynton/home/rotation/ckuhn/Desktop/Fraser_lab/water-scripts/4_atom_scripts/error/error\$TASK_ID.txt
#$ -e /wynton/home/rotation/ckuhn/Desktop/Fraser_lab/water-scripts/4_atom_scripts/out/\$TASK_ID.log
source= /wynton/home/rotation/ckuhn/anaconda3/envs
source activate water
base_dir='/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/water-scripts/4_atom_scripts'
out_dir='/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/water-scripts/4_atom_scripts/out'
file_dir='/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/top2018_cifs_full_filtered_hom30'
cd $base_dir
python /wynton/home/rotation/ckuhn/Desktop/Fraser_lab/water-scripts/4_atom_scripts/find_wat.py --pdb=$file_dir --out=$out_dir
