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

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="inpfile", help="Input parameter file", required=True)
args = parser.parse_args(
    "-i /pubhome/xtzhang/interaction/FF/Drude/run_drude/md_scripts/cmx.inp".split()
)
#%%
# Load parameters
os.chdir("/pubhome/xtzhang/interaction/FF/Drude/run_drude/produce_parameter/Classical2Drude/")
print("Loading parameters")
inputs = read_inputs(args.inpfile)
params = read_params(inputs.toppar)
if inputs.ftype == "drude" or inputs.ftype == "charmm":
    psf = read_psf(inputs.psffile)
    pdb = read_pdb(inputs.pdbfile)
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
print(total_coulomb - solute_coulomb - solvent_coulomb)
print(context.getState(getEnergy=True, groups={1}).getPotentialEnergy())

#%%
for i, f in enumerate(system.getForces()):
    state = context.getState(getEnergy=True, groups={i})
    print(f.getName(), state.getPotentialEnergy())
# %%
