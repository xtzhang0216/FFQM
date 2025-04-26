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
    if inputs.start == 0:
        pdb = read_pdb("reduce_restraint_npt_" + str(inputs.turn - 1) + ".pdb")
    else:
        pdb = read_pdb("md_" + inputs.pdbid + "-" + str(inputs.start - 1) + ".pdb")
    psf = gen_box_pdb(psf, pdb)

# Build system
if inputs.vdw == "Switch":
    system = psf.createSystem(
        params,
        nonbondedMethod=inputs.coulomb,
        nonbondedCutoff=inputs.r_off * nanometers,
        switchDistance=inputs.r_on * nanometers,
        # constraints=inputs.cons,
        ewaldErrorTolerance=inputs.ewald_Tol,
    )
elif inputs.vdw == "Force-switch":
    system = psf.createSystem(
        params,
        nonbondedMethod=inputs.coulomb,
        nonbondedCutoff=inputs.r_off * nanometers,
        # constraints=inputs.cons,
        ewaldErrorTolerance=inputs.ewald_Tol,
    )
    system = vfswitch(system, psf, inputs)
if inputs.pcouple == "yes":
    system = barostat(system, inputs)
if inputs.rest == "yes":
    system = restraints(system, pdb, inputs, inputs.turn, "MD")
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
if inputs.start == 0:
    if inputs.turn != 0:
        if os.path.exists("reduce_restraint_npt_" + str(inputs.turn - 1) + ".xml"):
            simulation.loadState(
                "reduce_restraint_npt_" + str(inputs.turn - 1) + ".xml"
            )
        elif os.path.exists("reduce_restraint_npt_" + str(inputs.turn - 1) + ".chk"):
            simulation.loadCheckpoint(
                "reduce_restraint_npt_" + str(inputs.turn - 1) + ".chk"
            )
        else:
            raise "NO NPT File For Restart"
    else:
        if os.path.exists("npt_" + inputs.pdbid + ".xml"):
            simulation.loadState("npt_" + inputs.pdbid + ".xml")
        elif os.path.exists("npt_" + inputs.pdbid + ".chk"):
            simulation.loadCheckpoint("npt_" + inputs.pdbid + ".chk")
        raise "NO NPT File For Restart"
else:
    if os.path.exists("md_" + inputs.pdbid + "-" + str(inputs.start - 1) + ".xml"):
        simulation.loadState(
            "md_" + inputs.pdbid + "-" + str(inputs.start - 1) + ".xml"
        )
    elif os.path.exists("md_" + inputs.pdbid + "-" + str(inputs.start - 1) + ".chk"):
        simulation.loadCheckpoint(
            "md_" + inputs.pdbid + "-" + str(inputs.start - 1) + ".chk"
        )
    else:
        raise "NO MD File For Restart"

# Drude VirtualSites
if inputs.ftype == "drude":
    simulation.context.computeVirtualSites()

# Production
for i in range(inputs.start, inputs.end):
    if i - 2 > 0 and inputs.start == 0:
        os.system("rm md_%s-%s*" % (inputs.pdbid, str(i - 2)))
    print(">>>>> Simulation %d." % i)
    odcd = "md_" + inputs.pdbid + "-" + str(i) + ".dcd"
    other = "md_" + inputs.pdbid + "-" + str(i) + ".csv"
    ochk = "md_" + inputs.pdbid + "-" + str(i) + ".chk"
    orst = "md_" + inputs.pdbid + "-" + str(i) + ".xml"
    opdb = "md_" + inputs.pdbid + "-" + str(i) + ".pdb"

    if inputs.nstep > 0:
        print("\nMD run: %s steps" % inputs.nstep)
        simulation.reporters.append(
            StateDataReporter(
                other,
                inputs.nstout,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                progress=True,
                remainingTime=True,
                speed=True,
                totalSteps=inputs.nstep,
                volume=True,
                separator="\t",
            )
        )
        if inputs.nstdcd > 0:
            simulation.reporters.append(
                DCDReporter(odcd, inputs.nstdcd, enforcePeriodicBox=True)
            )
        simulation.reporters.append(
            PDBReporter(opdb, inputs.nstep, enforcePeriodicBox=True)
        )
        simulation.reporters.append(CheckpointReporter(ochk, inputs.nstep))
        simulation.step(inputs.nstep)
        simulation.saveState(orst)
