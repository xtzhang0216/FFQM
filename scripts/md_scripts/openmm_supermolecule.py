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
#%%
def calc_energy(inp):
    os.chdir("/pubhome/xtzhang/interaction/FF/Drude/run_drude/produce_parameter/Classical2Drude/")
    inputs = read_inputs(inp)
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



    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    # Drude VirtualSites
    if inputs.ftype == "drude":
        context.computeVirtualSites()
    # Calculate initial system energy
    # print("\nInitial system energy")
    # print(context.getState(getEnergy=True).getPotentialEnergy())
    return context.getState(getEnergy=True).getPotentialEnergy()
# %%
total = calc_energy('/pubhome/xtzhang/interaction/FF/Drude/run_drude/md_scripts/cmx.inp')
fraga = calc_energy('/pubhome/xtzhang/interaction/FF/Drude/run_drude/md_scripts/fraga.inp')
fragb = calc_energy('/pubhome/xtzhang/interaction/FF/Drude/run_drude/md_scripts/fragb.inp')
inter = total - fraga - fragb
print(f"inter: {inter}")
# %%
