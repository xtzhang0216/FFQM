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
args = parser.parse_args()

# Load parameters
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

# Set platform
platform = Platform.getPlatformByName("CUDA")
simulation = Simulation(psf.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

# Drude VirtualSites
if inputs.ftype == "drude":
    simulation.context.computeVirtualSites()

# Calculate initial system energy
print("\nInitial system energy")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Energy minimization
if inputs.start == 0:
    if inputs.mini_nstep > 0:
        print("\nEnergy minimization: %s steps" % inputs.mini_nstep)
        simulation.minimizeEnergy(
            tolerance=inputs.mini_Tol * kilojoule / mole,
            maxIterations=inputs.mini_nstep,
        )
        print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
        crd = simulation.context.getState(
            getPositions=True, enforcePeriodicBox=True
        ).getPositions()
        PDBFile.writeFile(psf.topology, crd, open("em_" + inputs.pdbid + ".pdb", "w"))
        simulation.saveState("em_" + inputs.pdbid + ".xml")

