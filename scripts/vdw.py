#%%
import os
import shutil
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import AllChem
import subprocess
import pandas as pd
import json
from pathlib import Path
import sys
from openmm.app import AmberPrmtopFile, AmberInpcrdFile
from openmm import System, Context, Platform
from openmm.unit import *
from openmm.app import NoCutoff,PDBFile
from openmm import VerletIntegrator
from openmm.app import ForceField
import datetime


pdb_file = "/pubhome/xtzhang/interaction/FF/mol2/ACEM.pdb"
psf = ForceField('/pubhome/xtzhang/interaction/FF/opls/param/opls.zxt.xml','tip3p.xml') 

pdb = PDBFile(pdb_file)
for res in pdb.topology.residues():
    print(f"PDB Residue: {res.name}")
system = psf.createSystem(
    pdb.topology,
    nonbondedMethod=NoCutoff,
    constraints=None,
)
integrator = VerletIntegrator(0.001 * picoseconds)  # 步长设为极小值（不影响单点能量计算）
context = Context(system, integrator)
context.setPositions(pdb.positions)
state = context.getState(getEnergy=True)
potential_energy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
# %%
from openmm.app import *
from openmm import *
nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]

charges,sigmas,epsilons = [],[],[]
for i in range(system.getNumParticles()):
    charge, sigma, epsilon = nonbonded.getParticleParameters(i)
    charges.append(charge.value_in_unit(unit.elementary_charge))
    sigmas.append(sigma.value_in_unit(unit.nanometers))
    epsilons.append(epsilon.value_in_unit(unit.kilojoules_per_mole))
    
# %%
