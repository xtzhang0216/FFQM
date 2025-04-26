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
parser.add_argument(
    "-t", dest="turn", help="Nth run for restraint reduce", required=True
)
args = parser.parse_args()
turn = int(args.turn)

# Load parameters
print("Loading parameters")
inputs = read_inputs(args.inpfile)
params = read_params(inputs.toppar)
if inputs.ftype == "drude" or inputs.ftype == "charmm":
    psf = read_psf(inputs.psffile)
    if turn == 1:
        pdb = read_pdb("npt_" + inputs.pdbid + ".pdb")
    else:
        pdb = read_pdb("reduce_restraint_npt_" + str(turn - 1) + ".pdb")
    psf = gen_box_pdb(psf, pdb)

# Build system
if inputs.vdw == "Switch":
    system = psf.createSystem(
        params,
        nonbondedMethod=inputs.coulomb,
        nonbondedCutoff=inputs.r_off * nanometers,
        switchDistance=inputs.r_on * nanometers,
        ewaldErrorTolerance=inputs.ewald_Tol,
    )
elif inputs.vdw == "Force-switch":
    system = psf.createSystem(
        params,
        nonbondedMethod=inputs.coulomb,
        nonbondedCutoff=inputs.r_off * nanometers,
        ewaldErrorTolerance=inputs.ewald_Tol,
    )
    system = vfswitch(system, psf, inputs)
if inputs.pcouple == "yes":
    system = barostat(system, inputs)

system = restraints(system, pdb, inputs, turn, "REDUCE")
if inputs.rest_lig == "yes":
    system = restraint_ligand(system, pdb, inputs)

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
if turn == 1:
    if os.path.exists("npt_" + inputs.pdbid + ".xml"):
        simulation.loadState("npt_" + inputs.pdbid + ".xml")
    elif os.path.exists("npt_" + inputs.pdbid + ".chk"):
        simulation.loadCheckpoint("npt_" + inputs.pdbid + ".chk")
    else:
        raise "NO File For Restart"
else:
    if os.path.exists("reduce_restraint_npt_" + str(turn - 1) + ".xml"):
        simulation.loadState("reduce_restraint_npt_" + str(turn - 1) + ".xml")
    elif os.path.exists("reduce_restraint_npt_" + str(turn - 1) + ".chk"):
        simulation.loadCheckpoint("reduce_restraint_npt_" + str(turn - 1) + ".chk")
    else:
        raise "NO File For Restart"

# Drude VirtualSites
if inputs.ftype == "drude":
    simulation.context.computeVirtualSites()

# Production
odcd = "reduce_restraint_npt_" + str(turn) + ".dcd"
other = "reduce_restraint_npt_" + str(turn) + ".csv"
ochk = "reduce_restraint_npt_" + str(turn) + ".chk"
orst = "reduce_restraint_npt_" + str(turn) + ".xml"
opdb = "reduce_restraint_npt_" + str(turn) + ".pdb"

if inputs.reduce_nstep > 0:
    print("\nTurn %s reduce restraint NPT run: %s steps" % (turn, inputs.reduce_nstep))
    simulation.reporters.append(
        StateDataReporter(
            other,
            inputs.reduce_nstout,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=inputs.reduce_nstep,
            volume=True,
            separator="\t",
        )
    )
    simulation.reporters.append(
        DCDReporter(odcd, inputs.reduce_nstdcd, enforcePeriodicBox=True)
    )
    # Write restart file
    simulation.reporters.append(
        PDBReporter(opdb, inputs.reduce_nstep, enforcePeriodicBox=True)
    )
    simulation.reporters.append(CheckpointReporter(ochk, inputs.reduce_nstep))
    simulation.step(inputs.reduce_nstep)
    simulation.saveState(orst)
