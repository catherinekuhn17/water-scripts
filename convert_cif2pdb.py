import pymol2
import os

fn = '/wynton/home/rotation/ckuhn/Desktop/Fraser_lab/top2018_cifs_full_filtered_hom30' 

os.chdir(fn)
fns=[]
for fi in glob.glob("*/*/*"):
    fns.append(fi)
    
for infile in fns:
    with pymol2.PyMOL() as pymol:
         pymol.cmd.load(infile)
         pymol.cmd.save(infile.replace('.cif', '.pdb'))