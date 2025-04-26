#%%
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
os.chdir("/pubhome/xtzhang/interaction/FF/amoeba")
pdb = PDBFile('ala.pdb')
forcefield = ForceField('/pubhome/xtzhang/interaction/FF/amoeba/amoeba2018.zxt.xml')
# forcefield = ForceField('charmm36.xml')

system = forcefield.createSystem(pdb.topology,
                                 nonbondedMethod=NoCutoff,
                                #  polarisation='extrapolated',
                                #  mutualInducedTargetEpsilon=1e-5,
                                 constraints=None)

integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

state = simulation.context.getState(getEnergy=True)
print('Energy:', state.getPotentialEnergy())

# %%
