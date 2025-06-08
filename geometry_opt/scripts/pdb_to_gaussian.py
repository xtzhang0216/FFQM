#!/usr/bin/env python
import os
import glob
from pathlib import Path

# Define paths
pdb_dir = "/pubhome/xtzhang/interaction/FF/mol2"
output_dir = "/pubhome/xtzhang/interaction/FF/gaussian_opt"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Get all PDB files
pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))

for pdb_file in pdb_files:
    # Get the base name without extension
    base_name = os.path.basename(pdb_file).split('.')[0]
    output_file = os.path.join(output_dir, f"{base_name}.gjf")
    
    # Read PDB file
    with open(pdb_file, 'r') as f:
        pdb_lines = f.readlines()
    
    # Extract atom coordinates and element symbols
    atoms = []
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Extract element symbol (columns 77-78 in PDB format)
            if len(line) >= 78:
                element = line[76:78].strip()
            else:
                # If element field is not available, use the atom name
                element = line[12:16].strip()[0]
            
            # Extract x, y, z coordinates
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            
            atoms.append((element, x, y, z))
    
    # Create Gaussian input file
    with open(output_file, 'w') as f:
        # Write header
        f.write(f"%chk={base_name}.chk\n")
        f.write("%mem=8GB\n")
        f.write("%nprocshared=32\n")
        f.write("#p B3LYP/6-311G** Opt=cartesian\n\n")
        
        # Write title
        f.write(f"Geometry optimization of {base_name} using B3LYP/6-311G**\n\n")
        
        # Write charge and multiplicity (assuming neutral, singlet)
        f.write("0 1\n")
        
        # Write coordinates
        for element, x, y, z in atoms:
            f.write(f"{element:<2} {x:12.6f} {y:12.6f} {z:12.6f}\n")
        
        # Write blank line at the end
        f.write("\n")
    
    print(f"Created Gaussian input file: {output_file}")

print("All Gaussian input files have been created.")
