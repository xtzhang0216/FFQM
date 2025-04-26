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
#%%
for idx, mol in enumerate(sdf):
    if str(idx) not in ene.keys():
        continue
    print(f"Processing {idx}")
    fraga_name, fragb_name = mol.GetProp("FRAG_NAME").split()
    fga, fgb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    fraga = fragmol(fga, fraga_name, 1, "HETA")
    fragb = fragmol(fgb, fragb_name, 2, "HETB")
    fraga.write_crd(f"/pubhome/xtzhang/interaction/FF/drude/crd/step1_reader_seg1.crd")
    fragb.write_crd(f"/pubhome/xtzhang/interaction/FF/drude/crd/step1_reader_seg2.crd")
    os.chdir("/pubhome/xtzhang/interaction/FF/drude/crd")
    os.system("charmm -i step2_drude.inp")
    print(f"make crd file {idx} done")

    #%%
    # Load parameters
    psf_file = "/pubhome/xtzhang/interaction/FF/drude/crd/step2_drude.psf"
    pdb_file = "/pubhome/xtzhang/interaction/FF/drude/crd/step2_drude.pdb"
    os.chdir("/pubhome/xtzhang/interaction/FF/Drude/run_drude/produce_parameter/Classical2Drude/")
    print("Loading parameters")
    inputs = read_inputs(args.inpfile)
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

    #%%
    # for i, f in enumerate(system.getForces()):
    #     f.setForceGroup(i)
    #%%
    # Set platform
    # platform = Platform.getPlatformByName("CUDA")
    # simulation = Simulation(psf.topology, system, integrator, platform)

    #%%
    # # Calculate initial system energy
    # print("\nInitial system energy")
    # print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

    solvent = set([a.index for a in pdb.topology.atoms() if a.residue.name in ('ETOH')])
    protein = set([a.index for a in pdb.topology.atoms() if a.index not in solvent])

#%%
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            force.setForceGroup(0)
            force.addGlobalParameter("solute_scale", 1)
            force.addGlobalParameter("solvent_scale", 1)
            for i in range(force.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(i)
                # Set the parameters to be 0 when the corresponding parameter is 0,
                # and to have their normal values when it is 1.
                param = "solute_scale" if i in protein else "solvent_scale"
                force.setParticleParameters(i, 0, 0, 0)
                force.addParticleParameterOffset(param, i, charge, sigma, epsilon)
            for i in range(force.getNumExceptions()):
                p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
                force.setExceptionParameters(i, p1, p2, 0, 0, 0)
        elif isinstance(force, CustomNonbondedForce):
            force.setForceGroup(1)
            force.addInteractionGroup(protein, solvent)
        else:
            force.setForceGroup(2)

#%%
    context = Context(system, integrator)
    context.setPositions(pdb.positions)

# simulation = Simulation(psf.topology, system, integrator)
# simulation.context.setPositions(pdb.positions)

# Drude VirtualSites
# if inputs.ftype == "drude":
#     simulation.context.computeVirtualSites()
    context.computeVirtualSites()


#%%
    def coulomb_energy(solute_scale, solvent_scale):
        context.setParameter("solute_scale", solute_scale)
        context.setParameter("solvent_scale", solvent_scale)
        return context.getState(getEnergy=True, groups={0}).getPotentialEnergy()

    total_coulomb = coulomb_energy(1, 1)
    solute_coulomb = coulomb_energy(1, 0)
    solvent_coulomb = coulomb_energy(0, 1)
    inter_coulob = total_coulomb - solute_coulomb - solvent_coulomb
    inter_coulob_kcal = inter_coulob.value_in_unit(unit.kilocalorie/unit.mole)  # Convert kJ/mol to kcal/mol
    qm_ene = ene[str(idx)]
    results[str(idx)] = [inter_coulob_kcal, qm_ene]
    # Save results to a CSV file
results_df = pd.DataFrame.from_dict(results, orient='index', columns=['FF_Energy', 'QM_Energy'])
results_df.index.name = 'Index'
output_file = Path(f'/pubhome/xtzhang/interaction/FF/data/drude_ene/{pair}.csv')
results_df.to_csv(output_file, mode='a', header=not output_file.exists())
        
    # print(f"Intermolecular Coulomb energy: {inter_coulob_kcal} kcal/mol")
    # print(total_coulomb - solute_coulomb - solvent_coulomb)
    # print(context.getState(getEnergy=True, groups={1}).getPotentialEnergy())

#%%
# for i, f in enumerate(system.getForces()):
#     state = context.getState(getEnergy=True, groups={i})
#     print(f.getName(), state.getPotentialEnergy())
# %%
