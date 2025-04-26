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
ene = qm_energy_data["ETOH_MBZ_01.xyz"]
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



#%%
for idx, mol in enumerate(sdf):
    import time
    start_time = time.time()
    if str(idx) not in ene.keys():
        continue
    if idx > 1000:
        break
    print(f"Processing {idx}/{len(sdf)} at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    fraga_name, fragb_name = mol.GetProp("FRAG_NAME").split()
    fga, fgb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    fraga = fragmol(fga, fraga_name, 1, "HETA")
    fragb = fragmol(fgb, fragb_name, 2, "HETB")
    fraga.write_crd(f"/pubhome/xtzhang/interaction/FF/drude/crd/step1_reader_seg1.crd")
    fragb.write_crd(f"/pubhome/xtzhang/interaction/FF/drude/crd/step1_reader_seg2.crd")
    os.chdir("/pubhome/xtzhang/interaction/FF/drude/crd")
    os.system("charmm -i ene.inp > ene.out")
    # charmm_inter_ene = extract_interene_from_charmmout("ene.out")
    os.system("charmm -i fraga.inp > fraga.out")
    os.system("charmm -i fragb.inp > fragb.out")
    charmm_cmx_ene = extract_ene_from_charmmout("cmx.out")
    charmm_fraga_ene = extract_ene_from_charmmout("fraga.out")
    charmm_fragb_ene = extract_ene_from_charmmout("fragb.out")
    charmm_inter_ene = charmm_cmx_ene - charmm_fraga_ene - charmm_fragb_ene
    

    # print(f"make crd file {idx} done")

    # #%%
    # # Load parameters
    # psf_file = "/pubhome/xtzhang/interaction/FF/drude/crd/cmx.psf"
    # pdb_file = "/pubhome/xtzhang/interaction/FF/drude/crd/cmx.pdb"
    # fraga_psf_file = "/pubhome/xtzhang/interaction/FF/drude/crd/fraga.psf"
    # fraga_pdb_file = "/pubhome/xtzhang/interaction/FF/drude/crd/fraga.pdb"
    # fragb_psf_file = "/pubhome/xtzhang/interaction/FF/drude/crd/fragb.psf"
    # fragb_pdb_file = "/pubhome/xtzhang/interaction/FF/drude/crd/fragb.pdb"
    # os.chdir("/pubhome/xtzhang/interaction/FF/Drude/run_drude/produce_parameter/Classical2Drude/")
    # print("Loading parameters")
    # inp_file = args.inpfile
    # cmx_ene = get_ene(inp_file, psf_file, pdb_file)
    # fraga_ene = get_ene(inp_file, fraga_psf_file, fraga_pdb_file)
    # fragb_ene = get_ene(inp_file, fragb_psf_file, fragb_pdb_file)


    # # compute inter energy
    # inter_coulob = cmx_ene - fraga_ene - fragb_ene
    # inter_coulob_kcal = inter_coulob.value_in_unit(kilocalories_per_mole)

    qm_ene = ene[str(idx)]
    results[str(idx)] = [charmm_inter_ene, qm_ene]


results_df = pd.DataFrame.from_dict(results, orient='index', columns=['charmm_minimized_inter_energy', 'QM_Energy'])
results_df.index.name = 'Index'
output_file = Path(f'/pubhome/xtzhang/interaction/FF/drude/data/{pair}.csv')
results_df.to_csv(output_file, mode='a', header=not output_file.exists())
        
    # print(f"Intermolecular Coulomb energy: {inter_coulob_kcal} kcal/mol")
    # print(total_coulomb - solute_coulomb - solvent_coulomb)
    # print(context.getState(getEnergy=True, groups={1}).getPotentialEnergy())

#%%
# for i, f in enumerate(system.getForces()):
#     state = context.getState(getEnergy=True, groups={i})
#     print(f.getName(), state.getPotentialEnergy())
# %%
