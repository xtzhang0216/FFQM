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
#%%
sdfpath = "/pubhome/xtzhang/interaction/FF/outliar/tail/ACEM_NMA/amoeba/ACEM_nNMA_00_8372.sdf"
xyz = "ACET_PRPA_00"
start = 0
end = 100
amoeba_para = "/pubhome/xtzhang/interaction/FF/amoeba/para"
energy_json_path = Path('/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json')
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)
# for key in qm_energy_data.keys():
ene = qm_energy_data[f"{xyz}.xyz"]
results = {}

sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
# Settings
import tempfile
import time
WORKDIR = tempfile.mkdtemp(prefix=f"{xyz}_")
total = len(sdf)
mol = sdf[0]
idx=0
not_rot = [k for k in ene.keys() if not k.startswith('rot')]
for idx in range(start, end):
    mol = sdf[idx]

    # start_time = time.time()
    # with open(f"{WORKDIR}/log.txt", "a") as log_file:
    #     log_file.write(f"Processing conformer {idx + 1} / {total}\n")
    # if idx > 100:
    #     break
    if str(idx) not in not_rot:
        continue
    if mol is None:
        continue
    fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()
    # if fga_name == ref_fga_name and fgb_name == ref_fgb_name:
    #     pass
    # elif fga_name == ref_fgb_name and fgb_name == ref_fga_name:
    #     ref_fga_name, ref_fgb_name = fgb_name, fga_name
    #     ref_fga, ref_fgb = ref_fgb, ref_fga
    # else:
    #     print("Different conformer, skip")
    #     raise ValueError("Different conformer, skip")

    # 
    conf_dir = "/pubhome/xtzhang/interaction/FF/outliar/tail/ACEM_NMA/amoeba"

    print(f"Processing conformer {idx + 1} / {total} at {conf_dir} at {datetime.datetime.now()}", flush=True)
    os.makedirs(conf_dir, exist_ok=True)
    os.chdir(conf_dir)
    fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if fga_name == "ACEH" or fgb_name == "ACEH":
        continue
    if fga_name == "ACET" and fraga.GetNumAtoms() != 7:
        continue
    if fgb_name == "ACET" and fragb.GetNumAtoms() != 7:
        continue
    Chem.MolToPDBFile(fraga, f"{fga_name}_a.pdb")
    Chem.MolToPDBFile(fragb, f"{fgb_name}_b.pdb")
    Chem.MolToPDBFile(mol, f"{fga_name}_{fgb_name}.pdb")

    # cp para
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
    # pdb2xyz(f"{fga_name}_{fgb_name}", atom_type_a + atom_type_b, "cmx.key")
    e_a = calc_ene(f"{fga_name}_a", f"{fga_name}_a.xyz", f"{fga_name}.key")
    e_b = calc_ene(f"{fgb_name}_b", f"{fgb_name}_b.xyz", f"{fgb_name}.key")
    e_complex = calc_ene(f"{fga_name}_{fgb_name}", f"{fga_name}_{fgb_name}.xyz", "cmx.key")
    supermolecular_energy = e_complex - (e_a + e_b)
    inter_ene = extract_inter(f"{fga_name}_{fgb_name}.out")
    qm_ene = ene[str(idx)]
    # print(e_a, e_b, e_complex)
    # exit()
    results[str(idx)] = [round(supermolecular_energy, 2),round(inter_ene, 2), round(qm_ene, 2)]

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
# results_df.to_csv(output_file, mode='a', header=not output_file.exists())

# #%%
# rot_sdf_path = next((path for path in Path('/pubhome/lzeng/data/pair25/rot_split/').glob(f'*{xyz}*')), None)
# if rot_sdf_path:
#     rot_sdf = Chem.SDMolSupplier(str(rot_sdf_path), removeHs=False)
# else:
#     print(f"No rotated SDF file found for {xyz}")
#     rot_sdf = None

# rot_ene = {k: v for k, v in ene.items() if k.startswith('rot')}
# total = len(rot_sdf)
# if rot_ene:
#     for idx, mol in enumerate(rot_sdf):

#         # start_time = time.time()
#         # with open(f"{WORKDIR}/log.txt", "a") as log_file:
#         #     log_file.write(f"Processing conformer {idx + 1} / {total}\n")
#         # if idx > 100:
#         #     break
#         if f"rot-{idx}" not in rot_ene.keys():
#             continue
#         if mol is None:
#             continue
#         fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()
#         # if fga_name == ref_fga_name and fgb_name == ref_fgb_name:
#         #     pass
#         # elif fga_name == ref_fgb_name and fgb_name == ref_fga_name:
#         #     ref_fga_name, ref_fgb_name = fgb_name, fga_name
#         #     ref_fga, ref_fgb = ref_fgb, ref_fga
#         # else:
#         #     print("Different conformer, skip")
#         #     raise ValueError("Different conformer, skip")

#         # 
#         conf_dir = os.path.join(WORKDIR, f"conf_{idx}")

#         print(f"Processing conformer {idx + 1} / {total} at {conf_dir} at {datetime.datetime.now()}", flush=True)
#         os.makedirs(conf_dir, exist_ok=True)
#         os.chdir(conf_dir)

#         fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
#         Chem.MolToPDBFile(fraga, f"{fga_name}_a.pdb")
#         Chem.MolToPDBFile(fragb, f"{fgb_name}_b.pdb")
#         Chem.MolToPDBFile(mol, f"{fga_name}_{fgb_name}.pdb")

#         # cp para
#         shutil.copy(f"{amoeba_para}/{fga_name}.xyz", conf_dir)
#         shutil.copy(f"{amoeba_para}/{fgb_name}.xyz", conf_dir)
#         shutil.copy(f"{amoeba_para}/{fga_name}.key", conf_dir)
#         shutil.copy(f"{amoeba_para}/{fgb_name}.key", conf_dir)
#         make_cmx_key(f"{fga_name}.key", f"{fgb_name}.key")

#         atom_type_a = get_atom_type(f"{fga_name}.xyz")
#         atom_type_b = get_atom_type(f"{fgb_name}.xyz")
#         assert len(atom_type_a) == len(fraga.GetAtoms())
#         assert len(atom_type_b) == len(fragb.GetAtoms())
#         pdb2xyz(f"{fga_name}_a", atom_type_a, f"{fga_name}.key")
#         pdb2xyz(f"{fgb_name}_b", atom_type_b, f"{fgb_name}.key")
#         pdb2xyz(f"{fga_name}_{fgb_name}", atom_type_a + atom_type_b, "cmx.key")
#         e_a = calc_ene(f"{fga_name}_a", f"{fga_name}_a.xyz", f"{fga_name}.key")
#         e_b = calc_ene(f"{fgb_name}_b", f"{fgb_name}_b.xyz", f"{fgb_name}.key")
#         e_complex = calc_ene(f"{fga_name}_{fgb_name}", f"{fga_name}_{fgb_name}.xyz", "cmx.key")
#         interaction_energy = e_complex - (e_a + e_b)
#         qm_ene = rot_ene[f"rot-{idx}"]
#         results[f"rot-{idx}"] = [round(interaction_energy, 2), round(qm_ene, 2)]

#         os.chdir("../")
#     results_df = pd.DataFrame.from_dict(results, orient='index', columns=['FF_Energy', 'QM_Energy'])
#     results_df.index.name = 'Index'
#     if fgb_name[0] < fga_name[0]:
#         pair = f"{fgb_name}_{fga_name}"
#     else:
#         pair = f"{fga_name}_{fgb_name}"
#     output_file = Path(f'/pubhome/xtzhang/interaction/FF/amoeba/data/{pair}.csv')
#     results_df.to_csv(output_file, mode='a', header=not output_file.exists())
            