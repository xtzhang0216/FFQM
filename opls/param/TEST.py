#%%
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
os.chdir("/pubhome/xtzhang/interaction/FF/opls/param")
pdb = PDBFile('NMA.pdb')
forcefield = ForceField('/pubhome/xtzhang/interaction/FF/opls/param/opls.zxt.xml','tip3p.xml')
# forcefield = ForceField('charmm36.xml')

system = forcefield.createSystem(pdb.topology,
                                 nonbondedMethod=NoCutoff,
                                 constraints=None)

integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

state = simulation.context.getState(getEnergy=True)
print('Energy:', state.getPotentialEnergy())


# %%
