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
context = Context(system, integrator)
context.setPositions(pdb.positions)
context.computeVirtualSites()

#%%
for i, f in enumerate(system.getForces()):
    state = context.getState(getEnergy=True, groups={i})
    print(f.getName(), state.getPotentialEnergy())
# %%
