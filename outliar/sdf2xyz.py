import os
import sys
sdf_path = sys.argv[1]
xyz_path = sys.argv[2]


from ase.io import read, write
atoms = read(sdf_path, format='sdf')
# os.makedirs('/pubhome/xtzhang/interaction/data/csd/csd_pairs/modified_tar/xyz', exist_ok=True)
write(  xyz_path, atoms)