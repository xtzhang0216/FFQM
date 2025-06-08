#%%
import os
import shutil
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import AllChem
import subprocess
import pandas as pd
import json
from pathlib import Path
import sys
from openmm.app import AmberPrmtopFile, AmberInpcrdFile
from openmm import System, Context, Platform
from openmm.unit import *
from openmm.app import NoCutoff,PDBFile
from openmm import VerletIntegrator
from openmm.app import ForceField
import datetime
import argparse

def calc_ene(name, xyz_file,key_file):
    # echo "e" | analyze.x xyz -k key and then get energy
    # subprocess.run(["echo", "'e'", "|", "analyze.x", xyz_file, "-k", key_file, ">", f"{name}.out"])
    os.system(f"echo 'e' | analyze.x {xyz_file} -k {key_file} > {name}.out")
    with open(f"{name}.out", "r") as f:
        for line in f:
            if "Total Potential Energy : " in line:
            
                    return float(line.split()[4])
    return None

def extract_inter(cmx_out):
    with open(cmx_out, "r") as f:
        for line in f:
            if "Intermolecular Energy : " in line:
                return float(line.split()[3])
    return None

def pdb2xyz(query_name, atom_types, ref_key):
    #  pdbxyz.x pdb_file -k ref_key
    pdb_file = f"{query_name}.pdb"
    subprocess.run(["pdbxyz.x", pdb_file, "-k", ref_key])
    
    # with open(ref_xyz, "r") as f:
    #     lines = f.readlines()
    #     atom_types = [line.split()[5] for line in lines[1:]]
    
    # update atom_type in xyz file
    with open(f"{query_name}.xyz", "r") as f:
        lines = f.readlines()
    
    for i in range(1, len(lines)):
        lines[i] = lines[i][:56] + atom_types[i-1] + lines[i][60:]
    
    with open(f"{query_name}.xyz", "w") as f:
        f.writelines(lines)

def make_cmx_xyz(fga_xyz, fgb_xyz):
    namea = fga_xyz.split(".")[0].split("_")[0]
    nameb = fgb_xyz.split(".")[0].split("_")[0]
    with open(fga_xyz, "r") as f:
        fga_lines = f.readlines()
    with open(fgb_xyz, "r") as f:
        fgb_lines = f.readlines()
    len_fga = int(fga_lines[0])
    len_fgb = int(fgb_lines[0])
    total_atoms = len_fga + len_fgb
    
    combined_lines = [f"{total_atoms:>6}  xxx\n"]
    combined_lines.extend(fga_lines[1:])
    
    # make sure the column index is right
    for i, line in enumerate(fgb_lines[1:], start=1):
        parts = line.split()
        parts[0] = str(int(parts[0]) + len_fga)
        n_bonds = len(parts) - 6
        for j in range(1, n_bonds + 1):
            parts[5 + j] = str(int(parts[5 + j]) + len_fga)
        bonds_str = "".join(f"{part:>6}" for part in parts[6:])
        bonds_str = bonds_str[1:]
        new_line = f"{parts[0]:>6}  {parts[1]:<3}{parts[2]:>14}{parts[3]:>14}{parts[4]:>14}{parts[5]:>6}{bonds_str}\n"
        combined_lines.append(new_line)
    
    with open(f"{namea}_{nameb}.xyz", "w") as f:
        f.writelines(combined_lines)


def make_cmx_key(fga_key, fgb_key):
    with open(fga_key, "r") as f:
        fga_lines = f.readlines()
    with open(fgb_key, "r") as f:
        fgb_lines = f.readlines()
    combined_lines = fga_lines + [line for line in fgb_lines if line not in fga_lines]
    
    with open("cmx.key", "w") as f:
        f.writelines(combined_lines)
    
def get_atom_type(xyz_file):
    with open(xyz_file, "r") as f:
        lines = f.readlines()
        atom_types = [line.split()[5] for line in lines[1:]]
    return atom_types

def check_nh_distance(fragA, fragB, fragA_name, idx):
    if fragA_name not in ["ACEM", "MIMM", "MIME", "MIMD"]:
        return True
    for atomA in fragA.GetAtoms():
        if atomA.GetSymbol() == 'N':
            for neighbor in atomA.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    posH = fragA.GetConformer().GetAtomPosition(neighbor.GetIdx())
                    for atomB in fragB.GetAtoms():
                        posB = fragB.GetConformer().GetAtomPosition(atomB.GetIdx())
                        dist = ((posH.x-posB.x)**2 + (posH.y-posB.y)**2 + (posH.z-posB.z)**2)**0.5
                        if dist < 1.5:
                            print(f"[Warning] {fragA_name}-N-H...atom distance < 1.5A at conformer {idx}, skip.")
                            return False
    return True

def main():
    parser = argparse.ArgumentParser(description='CHARMM Drude energy calculation')
    parser.add_argument('--xyz', type=str, default="ACEM_ACEM_00", help='xyz name, e.g. ACEM_HOH_00')
    parser.add_argument('--start', type=int, default=0, help='start index for sdf conformers')
    parser.add_argument('--end', type=int, default=100, help='end index for sdf conformers (exclusive)')
    parser.add_argument('--check_nh', type=bool, default=True, help='check N-H distance for specific molecules')
    args = parser.parse_args()

    xyz = args.xyz
    start = args.start
    end = args.end

    amoeba_para = "/pubhome/xtzhang/interaction/FF/amoeba/para"
    energy_json_path = Path('/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json')
    with open(energy_json_path, 'r') as energy_file:
        qm_energy_data = json.load(energy_file)
    ene = qm_energy_data[f"{xyz}.xyz"]
    results = {}
    with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
        lines = f.readlines()
    sdfpath = next((line for line in lines if xyz in line), None)
    sdfpath = sdfpath.strip()
    sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
    import tempfile
    import time
    WORKDIR = tempfile.mkdtemp(prefix=f"{xyz}_")
    total = len(sdf)
    not_rot = [k for k in ene.keys() if not k.startswith('rot')]
    for idx in range(start, end):
        if str(idx) not in not_rot:
            continue
        mol = sdf[idx]
        if mol is None:
            continue
        fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()
        fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if fga_name == "ACEH" or fgb_name == "ACEH":
            continue
        if fga_name == "ACET" and fraga.GetNumAtoms() != 7:
            continue
        if fgb_name == "ACET" and fragb.GetNumAtoms() != 7:
            continue
        if args.check_nh:
            if not check_nh_distance(fraga, fragb, fga_name, idx) or not check_nh_distance(fragb, fraga, fgb_name, idx):
                continue
        conf_dir = os.path.join(WORKDIR, f"conf_{idx}")
        print(f"Processing conformer {idx} / {total} at {conf_dir} at {datetime.datetime.now()}", flush=True)
        os.makedirs(conf_dir, exist_ok=True)
        os.chdir(conf_dir)
        Chem.MolToPDBFile(fraga, f"{fga_name}_a.pdb")
        Chem.MolToPDBFile(fragb, f"{fgb_name}_b.pdb")
        Chem.MolToPDBFile(mol, f"{fga_name}_{fgb_name}.pdb")
        shutil.copy(f"{amoeba_para}/{fga_name}.xyz", conf_dir)
        shutil.copy(f"{amoeba_para}/{fgb_name}.xyz", conf_dir)
        shutil.copy(f"{amoeba_para}/{fga_name}.key", conf_dir)
        shutil.copy(f"{amoeba_para}/{fgb_name}.key", conf_dir)
        make_cmx_key(f"{fga_name}.key", f"{fgb_name}.key")
        atom_type_a = get_atom_type(f"{fga_name}.xyz")
        atom_type_b = get_atom_type(f"{fgb_name}.xyz")
        assert len(atom_type_a) == len(fraga.GetAtoms())
        assert len(atom_type_b) == len(fragb.GetAtoms())
        pdb2xyz(f"{fga_name}_a", atom_type_a, f"{fga_name}.key")
        pdb2xyz(f"{fgb_name}_b", atom_type_b, f"{fgb_name}.key")
        make_cmx_xyz(f"{fga_name}_a.xyz", f"{fgb_name}_b.xyz")
        e_a = calc_ene(f"{fga_name}_a", f"{fga_name}_a.xyz", f"{fga_name}.key")
        e_b = calc_ene(f"{fgb_name}_b", f"{fgb_name}_b.xyz", f"{fgb_name}.key")
        e_complex = calc_ene(f"{fga_name}_{fgb_name}", f"{fga_name}_{fgb_name}.xyz", "cmx.key")
        supermolecular_energy = e_complex - (e_a + e_b)
        inter_ene = extract_inter(f"{fga_name}_{fgb_name}.out")
        qm_ene = ene[str(idx)]
        results[str(idx)] = [round(supermolecular_energy, 2), round(inter_ene, 2), round(qm_ene, 2)]
        os.chdir("../")
        shutil.rmtree(conf_dir)
    os.chdir("/")
    shutil.rmtree(WORKDIR)
    results_df = pd.DataFrame.from_dict(results, orient='index', columns=['Tinker_Supermolecular_Energy', 'Tinker_Intermolecular_Energy', 'QM_Energy'])
    results_df.index.name = 'Index'
    if fgb_name < fga_name:
        pair = f"{fgb_name}_{fga_name}"
    else:
        pair = f"{fga_name}_{fgb_name}"
    output_file = Path(f'/pubhome/xtzhang/interaction/FF/amoeba/data/{pair}.csv')
    results_df.to_csv(output_file, mode='a', header=not output_file.exists())

if __name__ == "__main__":
    main()
# ...existing code...
