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

#%%
     
def mapping(frag, query,id):
    ref_mol_folder = "/pubhome/xtzhang/interaction/FF/amber/parameter"

    with open(os.path.join(ref_mol_folder, f"{frag}.mol2"), "r") as f:
        line1 = f.readlines()
    # Modify coordinates in line1 based on uery):
    crd_start=8
    for i in range(query.GetNumAtoms()):
        atom = query.GetAtomWithIdx(i)
        query_symbol = atom.GetSymbol()
        ref_symbol = line1[crd_start + i].split()[1][0]
        if query_symbol != ref_symbol:
            raise ValueError(f"Atom symbol mismatch: {query_symbol} != {ref_symbol}")
        crd = query.GetConformer().GetAtomPosition(i)
        line_idx = crd_start + i
        line1[line_idx] = line1[line_idx][:16] + f"{crd[0]:11.4f}" + line1[line_idx][27:]
        line1[line_idx] = line1[line_idx][:27] + f"{crd[1]:11.4f}" + line1[line_idx][38:]
        line1[line_idx] = line1[line_idx][:38] + f"{crd[2]:11.4f}" + line1[line_idx][49:]
    with open(f"{frag}_{id}.mol2", "w") as f:
        f.writelines(line1)

def calc_ene( crd_file):
    # Load the topology and coordinates
    # psf = ForceField("/pubhome/xtzhang/interaction/FF/charmm/charmm36.zxt.xml")
    psf = ForceField('/pubhome/xtzhang/interaction/FF/opls/param/opls.zxt.xml','tip3p.xml') 

    pdb = PDBFile(crd_file)
    for res in pdb.topology.residues():
        print(f"PDB Residue: {res.name}")
    system = psf.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )
    integrator = VerletIntegrator(0.001 * picoseconds)  # 步长设为极小值（不影响单点能量计算）
    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    state = context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    return potential_energy
import xml.etree.ElementTree as ET



def calc_cmx_ene(top_file1, top_file2, crd_file):

    if top_file1 == top_file2:
        psf = ForceField(top_file1)
    else:
        psf = ForceField(top_file1, top_file2)
    # Load the topology and coordinates
    pdb = PDBFile(crd_file)
    system = psf.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )
    integrator = VerletIntegrator(0.001 * picoseconds)  # 步长设为极小值（不影响单点能量计算）
    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    state = context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    return potential_energy


# load input para: --fraga --fragb --xyz
xyz = sys.argv[1]
# xyz="ETAM_HOH_00"

def write_pdb(ref_mol, mol, pdb_file):
    """
    
    """
    for i in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(i).GetSymbol() != ref_mol.GetAtomWithIdx(i).GetSymbol():
            raise ValueError(f"Atom symbol mismatch at index {i}: {mol.GetAtomWithIdx(i).GetSymbol()} != {ref_mol.GetAtomWithIdx(i).GetSymbol()}")
        ref_pos = mol.GetConformer().GetAtomPosition(i)
        ref_mol.GetConformer().SetAtomPosition(i, ref_pos)
    Chem.MolToPDBFile(ref_mol, pdb_file)

def write_cmx_pdb(mol1,mol2,pdb_file):
    length1 = mol1.GetNumAtoms()
    length2 = mol2.GetNumAtoms()
    combined = Chem.CombineMols(mol1, mol2)
    Chem.MolToPDBFile(combined, pdb_file)
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    with open(pdb_file, 'w') as f:
        for i,line in enumerate(lines):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if i < length1:
                    f.write(line[:25] + "1" + line[26:])
                elif length1 <= i < length1 + length2:
                    f.write(line[:25] + "2" + line[26:])
            else:
                f.write(line)


    
opls_para = "/pubhome/xtzhang/interaction/FF/opls/para"
energy_json_path = Path('/pubhome/xtzhang/interaction/data/pdbpairs/total_inteng_comb.json')
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)
# for key in qm_energy_data.keys():
ene = qm_energy_data[f"{xyz}.xyz"]
results = {}
with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
    lines = f.readlines()
sdfpath = next((line for line in lines if xyz in line), None)
sdfpath = sdfpath.strip()
sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
# Settings
import tempfile
import time
WORKDIR = tempfile.mkdtemp(prefix=f"{xyz}_")
total = len(sdf)
mol = sdf[0]
idx=0
# ref_fga_name, ref_fgb_name = mol.GetProp("FRAG_NAME").split()
# ref_file1 = f"{opls_para}/{ref_fga_name}.pdb"
# ref_file2 = f"{opls_para}/{ref_fgb_name}.pdb"
# ref_fga = Chem.MolFromPDBFile(ref_file1, removeHs=False)
# ref_fgb = Chem.MolFromPDBFile(ref_file2, removeHs=False)
    
for idx, mol in enumerate(sdf):

    # start_time = time.time()
    # with open(f"{WORKDIR}/log.txt", "a") as log_file:
    #     log_file.write(f"Processing conformer {idx + 1} / {total}\n")
    # if idx > 50:
    #     break
    if str(idx) not in ene.keys():
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
    conf_dir = os.path.join(WORKDIR, f"conf_{idx}")

    print(f"Processing conformer {idx + 1} / {total} at {conf_dir} at {datetime.datetime.now()}", flush=True)
    os.makedirs(conf_dir, exist_ok=True)
    os.chdir(conf_dir)

    fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    Chem.MolToPDBFile(fraga, f"{fga_name}_a.pdb")
    Chem.MolToPDBFile(fragb, f"{fgb_name}_b.pdb")
    e_a = calc_ene(f"{fga_name}_a.pdb")
    e_b = calc_ene(f"{fgb_name}_b.pdb")
    write_cmx_pdb(fraga, fragb, f"{fga_name}_{fgb_name}.pdb")
    e_complex = calc_ene(f"{fga_name}_{fgb_name}.pdb")
    interaction_energy = e_complex - (e_a + e_b)
    qm_ene = ene[str(idx)]
    results[str(idx)] = [round(interaction_energy, 2), round(qm_ene, 2)]

    os.chdir("../")
results_df = pd.DataFrame.from_dict(results, orient='index', columns=['FF_Energy', 'QM_Energy'])
results_df.index.name = 'Index'
if fgb_name[0] < fga_name[0]:
    pair = f"{fgb_name}_{fga_name}"
else:
    pair = f"{fga_name}_{fgb_name}"
output_file = Path(f'/pubhome/xtzhang/interaction/FF/opls/data/{pair}.csv')
results_df.to_csv(output_file, mode='a', header=not output_file.exists())
        