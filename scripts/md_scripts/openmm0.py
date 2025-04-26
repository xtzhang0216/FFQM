#%%
from __future__ import print_function
import argparse
import sys
import os

from omm_readinputs import *
from omm_readparams import *
from omm_vfswitch import *

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

from rdkit import Chem
import pandas as pd
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="inpfile", help="Input parameter file", required=True)
args = parser.parse_args(
    "-i /pubhome/xtzhang/interaction/FF/Drude/run_drude/md_scripts/cmx.inp".split()
)


#%%


class fragmol:
    def __init__(self, mol, resname, segid, het):
        self.mol = mol
        self.resname = resname
        self.segid = segid
        self.n_atoms = mol.GetNumAtoms()
        self.het = het
        self.ref_mol = Chem.MolFromPDBFile(f"/pubhome/xtzhang/interaction/FF/drude/ref_mols/{resname}.pdb",removeHs=False)
        with open(f"/pubhome/xtzhang/interaction/FF/drude/ref_mols/{resname}.pdb") as f:
            self.drudename = f.readline().split()[3]

    def assign_crd(self):
        match = self.ref_mol.GetSubstructMatch(self.mol)
        for atom in match:
            pos = self.mol.GetConformer().GetAtomPosition(atom)
            self.ref_mol.GetConformer().SetAtomPosition(atom, pos)

    def write_pdb(self, pdb_file):
        self.assign_crd()
        Chem.MolToPDBFile(self.ref_mol, pdb_file)    

    def write_crd(self, crd_file):
        self.assign_crd()
        atoms = [f"{self.n_atoms:>10}  EXT"]
        for i, atom in enumerate(self.ref_mol.GetAtoms(), start=1):
            pos = self.ref_mol.GetConformer().GetAtomPosition(atom.GetIdx())
            atom_name = atom.GetMonomerInfo().GetName()
            line = (
                f"{i:>10}{1:>10}{self.drudename:>6}{atom_name:>9}"
                f"{pos.x:>25.9f}{pos.y:>20.9f}{pos.z:>20.9f}  "
                f"{self.het}      1               0.0000000000"
            )
            atoms.append(line)
        with open(crd_file, "w") as f:
            f.write("\n".join(atoms))
#%%
results = {}
pair="ETOH_MBZ"
sdfpath = "/pubhome/lzeng/data/pair25/del_CA_THR_ETOH/delCA_ETOH_MBZ_01_noproxim.sdf"
sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
# mol = sdf[0]
energy_json_path = Path('/pubhome/xtzhang/interaction/data/pdbpairs/total_inteng_comb.json')
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)
# for key in qm_energy_data.keys():
ene = qm_energy_data["N1PA_N1PA.xyz"]
def extract_ene_from_charmmout(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if "ENER>" in line:
            energy = float(line.split()[2])
            return energy
    return None
def extract_interene_from_charmmout(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if "INTE>" in line:
            energy = float(line.split()[2])
            return energy
    return None

def get_ene(inp_file, psf_file, pdb_file):
    inputs = read_inputs(inp_file)
    params = read_params(inputs.toppar)
    if inputs.ftype == "drude" or inputs.ftype == "charmm":
        psf = read_psf(psf_file)
        pdb = read_pdb(pdb_file)
        psf = gen_box_pdb(psf, pdb)

    # Build system
    if inputs.vdw == "Switch":
        system = psf.createSystem(
            params,
            nonbondedMethod=inputs.coulomb,
            nonbondedCutoff=inputs.r_off * nanometers,
            switchDistance=inputs.r_on * nanometers,
            constraints=None,
            ewaldErrorTolerance=inputs.ewald_Tol,
        )
    elif inputs.vdw == "Force-switch":
        system = psf.createSystem(
            params,
            nonbondedMethod=inputs.coulomb,
            nonbondedCutoff=inputs.r_off * nanometers,
            constraints=inputs.cons,
            ewaldErrorTolerance=inputs.ewald_Tol,
        )
        system = vfswitch(system, psf, inputs)
    # system.addForce(AndersenThermostat(inputs.temp*kelvin, 1/picosecond))

    if inputs.ftype == "drude":
        integrator = DrudeLangevinIntegrator(
            inputs.temp * kelvin,
            inputs.fric_coeff / picosecond,
            inputs.drude_temp * kelvin,
            inputs.drude_fric_coeff / picosecond,
            inputs.dt * picoseconds,
        )
        integrator.setMaxDrudeDistance(inputs.drude_hardwall)
    if inputs.ftype == "charmm":
        integrator = LangevinIntegrator(
            inputs.temp * kelvin, inputs.fric_coeff / picosecond, inputs.dt * picoseconds
        )

    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    # Drude VirtualSites
    if inputs.ftype == "drude":
        context.computeVirtualSites()
    ene = context.getState(getEnergy=True).getPotentialEnergy()
    return ene

def set_residue_info(mol, res_name, res_number):
    for atom in mol.GetAtoms():
        # 创建新的PDB残基信息对象
        pdb_info = Chem.AtomPDBResidueInfo()
        pdb_info.SetResidueName(res_name)  # PDB要求残基名为3字符
        pdb_info.SetResidueNumber(res_number)
        atom.SetMonomerInfo(pdb_info)
    return mol

def get_classicFF_ene(pdb_file):
    pdb = PDBFile(pdb_file)
    forcefield = ForceField('amber/ff14SB.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
    integrator = VerletIntegrator(0.001*picosecond)
    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    ene = context.getState(getEnergy=True).getPotentialEnergy()
    return ene

def assign_crd(ref_mol, mol):
    match = ref_mol.GetSubstructMatch(mol)
    for atom in match:
        pos = mol.GetConformer().GetAtomPosition(atom)
        ref_mol.GetConformer().SetAtomPosition(atom, pos)
    return ref_mol
#%%
fraga_name = "ETOH"
fragb_name = "MBZ"
fraga_ref_mol = Chem.MolFromPDBFile(f"/pubhome/xtzhang/interaction/FF/charmm/template/{fraga_name}.pdb",removeHs=False)
with open(f"/pubhome/xtzhang/interaction/FF/charmm/template/{fraga_name}.pdb") as f:
    ref_fraga_name = f.readline().split()[3]
fragb_ref_mol = Chem.MolFromPDBFile(f"/pubhome/xtzhang/interaction/FF/charmm/template/{fragb_name}.pdb",removeHs=False)
with open(f"/pubhome/xtzhang/interaction/FF/charmm/template/{fragb_name}.pdb") as f:
    ref_fragb_name = f.readline().split()[3]
for idx, mol in enumerate(sdf):
    if idx > 500:
        break
    if str(idx) not in ene.keys():
        continue
    fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)

    Chem.MolToPDBFile(fraga, "/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/fraga.pdb")
    Chem.MolToPDBFile(fragb, "/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/fragb.pdb")


    # # assign crd to ref
    # fraga = assign_crd(fraga_ref_mol, fraga)
    # fragb = assign_crd(fragb_ref_mol, fragb)


    # # # # save to pdb
    # tmp_fraga_path = "/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/fraga.pdb"
    # tmp_fragb_path = "/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/fragb.pdb"
    # Chem.MolToPDBFile(fraga, tmp_fraga_path)
    # Chem.MolToPDBFile(fragb, tmp_fragb_path)
    # import re
    # # 定义文件路径
    # # 读取文件内容
    # with open(tmp_fraga_path, 'r') as file:
    #     lines = file.readlines()
    # # 处理每一行，替换 18-21 列的内容
    # with open(tmp_fraga_path, 'w') as file:
    #     for line in lines:
    #         if not line.startswith("HETATM"):
    #             file.write(line)
    #             continue
    #         modified_line = re.sub(r'^(.{17}).{4}', r'\1' + ref_fraga_name, line)
    #         file.write(modified_line)
    # with open(tmp_fragb_path, 'r') as file:
    #     lines = file.readlines()
    # # 处理每一行，替换 18-23 列的内容
    # with open(tmp_fragb_path, 'w') as file:
    #     for line in lines:
    #         if not line.startswith("HETATM"):
    #             file.write(line)
    #             continue

    #         modified_line = re.sub(r'^(.{17}).{4}', r'\1' + ref_fragb_name, line)
    #         file.write(modified_line)

    # combined_mol = Chem.CombineMols(fraga, fragb)
    # Chem.MolToPDBFile(combined_mol, "/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/cmx.pdb")
    # with open("/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/cmx.pdb", 'r') as file:
    #     lines = file.readlines()
    # # 处理每一行，替换 18-23 列的内容
    # with open("/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/cmx.pdb", 'w') as file:
    #     for i, line in enumerate(lines):
    #         if not line.startswith("HETATM"):
    #             file.write(line)
    #             continue
    #         if i < 9 :
    #             modified_line = re.sub(r'^(.{17}).{4}', r'\1' + ref_fraga_name, line)
    #             file.write(modified_line)

    #         else:
    #             modified_line = re.sub(r'^(.{17}).{4}', r'\1' + ref_fragb_name, line)
    #             modified_line = modified_line[:25] + "2" + modified_line[26:]

    #             file.write(modified_line)
                

    fraga_energy = get_classicFF_ene("/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/fraga.pdb")
    fragb_energy = get_classicFF_ene("/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/fragb.pdb")
    cmx_energy = get_classicFF_ene("/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/cmx.pdb")
    classic_inter_ene = cmx_energy - fraga_energy - fragb_energy
    classic_inter_ene_kcal = classic_inter_ene.value_in_unit(kilocalories_per_mole)
    qm_ene = ene[str(idx)]

    results[str(idx)] = [classic_inter_ene_kcal, qm_ene]
#%%
# # save to pdb
# Chem.MolToPDBFile(mol, "/pubhome/xtzhang/interaction/FF/scripts/md_scripts/tmp/1.pdb")



#%%
# for idx, mol in enumerate(sdf):
#     if str(idx) not in ene.keys():
#         continue
#     print(f"Processing {idx}")
#     fraga_name, fragb_name = mol.GetProp("FRAG_NAME").split()
#     fga, fgb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
#     fraga = fragmol(fga, fraga_name, 1, "HETA")
#     fragb = fragmol(fgb, fragb_name, 2, "HETB")
#     fraga.write_crd(f"/pubhome/xtzhang/interaction/FF/drude/crd/step1_reader_seg1.crd")
#     fragb.write_crd(f"/pubhome/xtzhang/interaction/FF/drude/crd/step1_reader_seg2.crd")
#     os.chdir("/pubhome/xtzhang/interaction/FF/drude/crd")
#     os.system("charmm -i ene.inp > ene.out")
#     charmm_inter_ene = extract_interene_from_charmmout("ene.out")
#     # os.system("charmm -i fraga.inp > fraga.out")
#     # os.system("charmm -i fragb.inp > fragb.out")
#     # charmm_cmx_ene = extract_ene_from_charmmout("cmx.out")
#     # charmm_fraga_ene = extract_ene_from_charmmout("fraga.out")
#     # charmm_fragb_ene = extract_ene_from_charmmout("fragb.out")
#     # charmm_inter_ene = charmm_cmx_ene - charmm_fraga_ene - charmm_fragb_ene
    

#     # print(f"make crd file {idx} done")

#     # #%%
#     # # Load parameters
#     # psf_file = "/pubhome/xtzhang/interaction/FF/drude/crd/cmx.psf"
#     # pdb_file = "/pubhome/xtzhang/interaction/FF/drude/crd/cmx.pdb"
#     # fraga_psf_file = "/pubhome/xtzhang/interaction/FF/drude/crd/fraga.psf"
#     # fraga_pdb_file = "/pubhome/xtzhang/interaction/FF/drude/crd/fraga.pdb"
#     # fragb_psf_file = "/pubhome/xtzhang/interaction/FF/drude/crd/fragb.psf"
#     # fragb_pdb_file = "/pubhome/xtzhang/interaction/FF/drude/crd/fragb.pdb"
#     # os.chdir("/pubhome/xtzhang/interaction/FF/Drude/run_drude/produce_parameter/Classical2Drude/")
#     # print("Loading parameters")
#     # inp_file = args.inpfile
#     # cmx_ene = get_ene(inp_file, psf_file, pdb_file)
#     # fraga_ene = get_ene(inp_file, fraga_psf_file, fraga_pdb_file)
#     # fragb_ene = get_ene(inp_file, fragb_psf_file, fragb_pdb_file)


#     # # compute inter energy
#     # inter_coulob = cmx_ene - fraga_ene - fragb_ene
#     # inter_coulob_kcal = inter_coulob.value_in_unit(kilocalories_per_mole)

    # qm_ene = ene[str(idx)]
    # results[str(idx)] = [cmx_energy, qm_ene]


results_df = pd.DataFrame.from_dict(results, orient='index', columns=['charmm_inter_energy', 'QM_Energy'])
results_df.index.name = 'Index'
output_file = Path(f'/pubhome/xtzhang/interaction/FF/data/{pair}.csv')
results_df.to_csv(output_file, mode='a', header=not output_file.exists())
        
#     # print(f"Intermolecular Coulomb energy: {inter_coulob_kcal} kcal/mol")
#     # print(total_coulomb - solute_coulomb - solvent_coulomb)
#     # print(context.getState(getEnergy=True, groups={1}).getPotentialEnergy())

# #%%
# # for i, f in enumerate(system.getForces()):
# #     state = context.getState(getEnergy=True, groups={i})
# #     print(f.getName(), state.getPotentialEnergy())
# # %%
