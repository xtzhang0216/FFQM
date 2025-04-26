#%%
import os
import shutil
from rdkit import Chem
from rdkit.Chem import rdmolfiles
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

#%%
class fragmol:
    def __init__(self, mol, resname, segid, het):
        self.mol = mol
        self.resname = resname
        self.segid = segid
        self.n_atoms = mol.GetNumAtoms()
        self.het = het
        self.ref_mol = Chem.MolFromMol2File(f"/pubhome/xtzhang/interaction/FF/charmm/para/{resname}.cgenff.mol2",removeHs=False)
        self.atom_names = []
        with open (f"/pubhome/xtzhang/interaction/FF/charmm/para/{resname}.cgenff.mol2", "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("@<TRIPOS>ATOM"):
                    start_line = lines.index(line) + 1
                if line.startswith("@<TRIPOS>BOND"):
                    end_line = lines.index(line)
                    break
            for line in lines[start_line:end_line]:
                atom_name = line.split()[1]
                self.atom_names.append(atom_name)
        # with open(f"/pubhome/xtzhang/interaction/FF/drude/ref_mols/{resname}.pdb") as f:
        #     self.drudename = f.readline().split()[3]
    def assign_crd(self):
        match = self.ref_mol.GetSubstructMatch(self.mol)
        for ref_idx, mol_idx in enumerate(match):
            ref_pos = self.mol.GetConformer().GetAtomPosition(ref_idx)
            self.ref_mol.GetConformer().SetAtomPosition(mol_idx, ref_pos)
    def write_crd(self, crd_file):
        self.assign_crd()
        atoms = [f"{self.n_atoms:>10}  EXT"]
        for i, atom in enumerate(self.ref_mol.GetAtoms(), start=1):
            pos = self.ref_mol.GetConformer().GetAtomPosition(atom.GetIdx())
            # atom_name = atom.GetMonomerInfo().GetName()
            atom_name = self.atom_names[i-1]
            line = (
                f"{i:>10}{1:>10}{self.resname:>6}{atom_name:>9}"
                f"{pos.x:>25.9f}{pos.y:>20.9f}{pos.z:>20.9f}  "
                f"{self.het}      1               0.0000000000"
            )
            atoms.append(line)
        with open(crd_file, "w") as f:
            f.write("\n".join(atoms))
    def write_pdb(self, pdb_file):
        self.assign_crd()
        Chem.MolToPDBFile(self.ref_mol, pdb_file)
        with open(pdb_file, "r") as f:
            lines = f.readlines()
        with open(pdb_file, "w") as f:
            for line in lines:
                if "UNL " in line:
                    line = line.replace("UNL ", f"{self.resname:<4}")
                f.write(line)
        # atoms = []
        # for i, atom in enumerate(self.ref_mol.GetAtoms(), start=1):
        #     pos = self.ref_mol.GetConformer().GetAtomPosition(atom.GetIdx())
        #     # atom_name = atom.GetMonomerInfo().GetName()
        #     atom_name = self.atom_names[i-1]
        #     line = (
        #         f"ATOM  {i:>5} {atom_name:<4} {self.resname:<3} {self.segid:<4}"
        #         f"{pos.x:>8.3f}{pos.y:>8.3f}{pos.z:>8.3f}  "
        #         f"1.00  0.00"
        #     )
        #     atoms.append(line)
        # with open(pdb_file, "w") as f:
        #     f.write("\n".join(atoms))

def assign_crd(ref_mol, mol):
        match = ref_mol.GetSubstructMatch(mol)
        for ref_idx, mol_idx in enumerate(match):
            ref_pos = mol.GetConformer().GetAtomPosition(ref_idx)
            ref_mol.GetConformer().SetAtomPosition(mol_idx, ref_pos)
        return ref_mol

def write_pdb0(mol, pdb_file):
    resname = pdb_file.split("_")[0]
    Chem.MolToPDBFile(mol, pdb_file, flavor=2)
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    with open(pdb_file, "w") as f:
        for line in lines:
            if "UNL " in line:
                line = line.replace("UNL ", f"{resname:<4}")
            f.write(line)

def write_pdb(mol2_file, mol, pdb_file, id):
    """
    combine atom label in mol2 file lines and coordinate form mol to write pdb file
    """
    resid = id
    resname = pdb_file.split("_")[0]
    pdb_lines = []
    with open(mol2_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                start_line = lines.index(line) + 1
            if line.startswith("@<TRIPOS>BOND"):
                end_line = lines.index(line)
                break
    for i,line in enumerate(lines[start_line:end_line]):
        atom_name = line.split()[1]
        atom_idx = i + 1
        pos = mol.GetConformer().GetAtomPosition(i)
        line = (
            f"HETATM{atom_idx:>5}  {atom_name:<3} {resname:<4}    {resid:<5}"
            f"{pos.x:>8.3f}{pos.y:>8.3f}{pos.z:>8.3f}  "
            f"1.00  0.00"
        )
        pdb_lines.append(line)
    with open(pdb_file, "w") as f:
        f.write("\n".join(pdb_lines) + "\n")


def write_cmx_pdb(mola, molb, pdb_file):
    combined = Chem.CombineMols(mola, molb)
    Chem.MolToPDBFile(combined, pdb_file)
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    with open(pdb_file, "w") as f:
        for line in lines:
            if "UNL " in line:
                line = line.replace("UNL ", "CMX ")
            f.write(line)
#%%
# xyz = sys.argv[1]
xyz="N1PA_N1PA"

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
#%%
# Settings
import tempfile
import time
from parmed import load_file, charmm

WORKDIR = tempfile.mkdtemp(prefix=f"{xyz}_")
total = len(sdf)
mol = sdf[0]
idx=0
ref_fga_name, ref_fgb_name = mol.GetProp("FRAG_NAME").split()
ref_fga = Chem.MolFromMol2File(f"/pubhome/xtzhang/interaction/FF/charmm/para/{ref_fga_name}.cgenff.mol2",removeHs=False)
ref_fgb = Chem.MolFromMol2File(f"/pubhome/xtzhang/interaction/FF/charmm/para/{ref_fgb_name}.cgenff.mol2",removeHs=False)
mol = sdf[0]

conf_dir = os.path.join(WORKDIR, f"conf_{idx}")
print(f"Processing conformer {idx + 1} / {total} at {conf_dir}")
os.makedirs(conf_dir, exist_ok=True)
os.chdir(conf_dir)
fga, fgb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()
if fga_name == ref_fga_name and fgb_name == ref_fgb_name:
    pass
elif fga_name == ref_fgb_name and fgb_name == ref_fga_name:
    ref_fga_name, ref_fgb_name = fgb_name, fga_name
    ref_fga, ref_fgb = ref_fgb, ref_fga
else:
    print("Different conformer, skip")
    raise ValueError("Different conformer, skip")

# fga = assign_crd(fga, ref_fga)
# fga = assign_crd(fgb, ref_fgb)
write_pdb(f"/pubhome/xtzhang/interaction/FF/charmm/para/{fga_name}.cgenff.mol2", fga, f"{fga_name}_1.pdb", 1)
# write_pdb(fgb, f"{fgb_name}_2.pdb")
# write_cmx_pdb(fga, fgb, f"{fga_name}_a.pdb", f"{fgb_name}_b.pdb")

shutil.copy("/pubhome/xtzhang/interaction/FF/charmm/para/toppar/top_all36_cgenff.rtf", conf_dir)
shutil.copy("/pubhome/xtzhang/interaction/FF/charmm/para/toppar/par_all36_cgenff.prm", conf_dir)
shutil.copy(f"/pubhome/xtzhang/interaction/FF/charmm/para/{fga_name}.str", conf_dir)
charmm_top = charmm.CharmmParameterSet(
    'top_all36_cgenff.rtf',     
    'par_all36_cgenff.prm',     
    f'{fga_name}.str'        
)
struct = load_file(f'{fga_name}_1.pdb')

psf = charmm.CharmmPsfFile.from_structure(struct)
psf.save(os.path.join(conf_dir, f"{fga_name}_1.psf"))
#%%
def calc_ene(top_file, crd_file):
    psf = charmm.CharmmPsfFile(top_file)
    pdb = PDBFile(crd_file)
    system = psf.createSystem(
        params=charmm_top,
        nonbondedMethod=NoCutoff,
        constraints=None,
        implicitSolvent=None
    )
    integrator = VerletIntegrator(0.001 * picoseconds)  # 步长设为极小值（不影响单点能量计算）
    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    state = context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    return potential_energy
ene_a = calc_ene(f"{fga_name}_1.psf", f"{fga_name}_1.pdb")
ene_b = calc_ene(f"{fgb_name}_b.psf", f"{fgb_name}_b.pdb")


# %%
# calculate energy with openmm

