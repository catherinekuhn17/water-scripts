#!/bin/bash
#Stephanie Wankowicz
#04/25/2019

#this must be done before you submit to SGE since SGE cannot connect to the internet!

#________________________________________________INPUTS________________________________________________#
base_folder='/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/pdb_files/download3/' #base folder (where you want to put folders/pdb files

pdb_filelist=/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/pdb_files/download3/pdb_ids.txt
while read -r line; do
  PDB=$line
  cd $base_folder
  if [ -d "/$PDB" ]; then
    echo "Folder exists." 
  else
    mkdir $PDB
  fi
  #mkdir $PDB
  cd $PDB
  curl -L -O https://files.rcsb.org/download/${PDB}.pdb
  curl -L -O https://files.rcsb.org/download/${PDB}-sf.cif
  curl -L -O http://edmaps.rcsb.org/coefficients/${PDB}.mtz
done < $pdb_filelist
