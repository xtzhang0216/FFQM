from __future__ import print_function
import argparse
import sys
import os

from omm_readinputs import *
from omm_readparams import *
from omm_vfswitch import *
from omm_barostat import *
from omm_restraints import *

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
    pdb = read_pdb("nvt_" + inputs.pdbid + ".pdb")
    psf = gen_box_pdb(psf, pdb)

# Build system
if inputs.vdw == "Switch":
    system = psf.createSystem(
        params,
        nonbondedMethod=inputs.coulomb,
        nonbondedCutoff=inputs.r_off * nanometers,
        switchDistance=inputs.r_on * nanometers,
        constraints=inputs.cons,
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
if inputs.pcouple == "yes":
    system = barostat(system, inputs)
system = restraints(system, pdb, inputs, 0, "NPT")

if inputs.ftype == "drude":
    integrator = DrudeLangevinIntegrator(
        inputs.temp * kelvin,
        inputs.fric_coeff / picosecond,
        inputs.drude_temp * kelvin,
        inputs.drude_fric_coeff / picosecond,
        inputs.dt * picoseconds,
    )
    integrator.setMaxDrudeDistance(inputs.drude_hardwall * nanometers)  # Drude Hardwall
if inputs.ftype == "charmm":
    integrator = LangevinIntegrator(
        inputs.temp * kelvin,
        inputs.fric_coeff / picosecond,
        inputs.nvt_dt * picoseconds,
    )

# Set platform
platform = Platform.getPlatformByName("CUDA")

# Build simulation context
simulation = Simulation(psf.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions.value_in_unit(angstroms) * angstroms)
if os.path.exists("nvt_" + inputs.pdbid + ".xml"):
    simulation.loadState("nvt_" + inputs.pdbid + ".xml")
elif os.path.exists("nvt_" + inputs.pdbid + ".chk"):
    simulation.loadCheckpoint("nvt_" + inputs.pdbid + ".chk")
else:
    raise "NO NVT File For Restart"

# Drude VirtualSites
if inputs.ftype == "drude":
    simulation.context.computeVirtualSites()

# Production
odcd = "npt_" + inputs.pdbid + ".dcd"
other = "npt_" + inputs.pdbid + ".csv"
ochk = "npt_" + inputs.pdbid + ".chk"
orst = "npt_" + inputs.pdbid + ".xml"
opdb = "npt_" + inputs.pdbid + ".pdb"

if inputs.npt_nstep > 0:
    print("\nNPT run: %s steps" % inputs.npt_nstep)
    simulation.reporters.append(
        StateDataReporter(
            other,
            inputs.npt_nstout,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=inputs.npt_nstep,
            volume=True,
            separator="\t",
        )
    )
    simulation.reporters.append(
        DCDReporter(odcd, inputs.npt_nstdcd, enforcePeriodicBox=True)
    )
    # Write restart file
    simulation.reporters.append(
        PDBReporter(opdb, inputs.npt_nstep, enforcePeriodicBox=True)
    )
    simulation.reporters.append(CheckpointReporter(ochk, inputs.npt_nstep))
    simulation.step(inputs.npt_nstep)
    simulation.saveState(orst)
