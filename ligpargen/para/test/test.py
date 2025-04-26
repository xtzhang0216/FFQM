#%%
from openmm.app import AmberPrmtopFile, AmberInpcrdFile
from openmm import System, Context, Platform
from openmm.unit import *
from openmm.app import NoCutoff,PDBFile
from openmm import VerletIntegrator
from openmm.app import ForceField
import os
#%%

topa = "N1PA.xml"
topb = "PRPA.xml"
pdba = "N1PA_a.pdb"
pdbb = "PRPA_b.pdb"
pdb_cmx = "N1PA_PRPA.pdb"
#%%
os.chdir("/pubhome/xtzhang/interaction/FF/opls/para/test")
psf_cmx = ForceField(topa, topb)

# Load the topology and coordinates
crd_cmx = PDBFile(pdb_cmx)
system = psf_cmx.createSystem(
    crd_cmx.topology,
    nonbondedMethod=NoCutoff,
    constraints=None,
)
# %%
psfb = ForceField(topb)
crdb = PDBFile(pdbb)
systemb = psfb.createSystem(
    crdb.topology,
    nonbondedMethod=NoCutoff,
    constraints=None,
)

# %%
