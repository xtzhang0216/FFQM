#%%
import pandas as pd
df = pd.read_csv("/pubhome/xtzhang/interaction/FF/outliar/outliar_head/ETAM_MSM/10.csv")
df = df[df["QM_Energy"]<0]
#%%
from rdkit import Chem
mols = Chem.SDMolSupplier("/pubhome/xtzhang/interaction/FF/outliar/outliar_head/ETAM_MSM/ETAM_MSM_outliar.charmm.sdf", removeHs=False)
writer = Chem.SDWriter("/pubhome/xtzhang/interaction/FF/outliar/outliar_head/ETAM_MSM/ETAM_MSM_outliar.favor.charmm.sdf")
for mol in mols:
    if int(mol.GetProp("Index")) in df["Index"].values and float(mol.GetProp("QM_Energy")) < 0:
        diff = float(mol.GetProp("FF_Energy")) - float(mol.GetProp("QM_Energy"))
        # if  diff > 10:
        writer.write(mol)
# %%