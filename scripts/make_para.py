#%%
import os
import subprocess
FragSmartsSmiles = {
    "ACEM": ("CC(N)=O", "CC(N)=O", 0),
    "ACET": ("CC([O])=O", "CC([O-])=O",-1),
    "MBZ": ("Cc1ccccc1", "Cc1ccccc1",0),
    "MIMD": ("[#6](-[#6]1-[#7](-[#6](-[#7]-[#6]-1-[H])-[H])-[H])(-[H])(-[H])-[H]", "[H]c1nc([H])n([H])c1C([H])([H])[H]",0),
    "MIME": ("[#6](-[#6]1-[#7]-[#6](-[#7](-[#6]-1-[H])-[H])-[H])(-[H])(-[H])-[H]", "[H]c1nc(c([H])([H])[H])c([H])n1[H]",0),
    "MIMM": ("[#6](-[#6]1-[#7](-[#6](-[#7](-[#6]-1-[H])-[H])-[H])-[H])(-[H])(-[H])-[H]", "n1c(C)c[nH+]c1",1),
    "MIND": ("Cc1cnc2ccccc12", "Cc1c[nH]c2ccccc12",0),
    "ETAM": ("CC[NH3+]", "CC[NH3+]",1),
    "ETOH": ("CCO", "CCO",0),
    "ETSH": ("CC[SH]", "CC[SH]",0),
    "MGDM": ("CNC(N)=[NH2]", "CNC(N)=[NH2+]",1),
    "MSM": ("CSC", "CSC",0),
    "NMA": ("CNC(~O)C", "CNC(=O)C",0),
    "MPHE": ("Cc1ccc(O)cc1", "Cc1ccc(O)cc1",0),
    "PRPA": ("C[CH1,CH2]C", "CCC",0),
    "N1PA": ("CC(~O)N1CCCC1", "CC(=O)N1CCCC1",0),
    "HOH": ("[#8](-[H])-[H]", "[H]O[H]",0)
}
#%%

os.chdir("/pubhome/xtzhang/interaction/FF/mol2")
for frag, (smarts, smiles,charge) in FragSmartsSmiles.items():
    
    with open(f"{frag}.smi", "w") as f:
        f.write(smiles)
    os.system(f"/pubhome/xtzhang/opt/unicon_1.4.2/unicon -i {frag}.smi -o {frag}.mol2 -g 3")
# %%
for frag, (smarts, smiles, netcharge) in FragSmartsSmiles.items():
    if frag != "ETSH":
        continue
    file = f"{frag}.sdf"
    subprocess.run([
        "/pubhome/xtzhang/mambaforge/envs/htmd/bin/antechamber", "-i", file, "-fi", "sdf",
        "-o", f"{frag}.mol2", "-fo", "mol2", "-c", "bcc", "-s", "2", "-nc", f"{netcharge}"
    ])
    subprocess.run([
        "/pubhome/xtzhang/mambaforge/envs/htmd/bin/parmchk2", "-i", f"{frag}.mol2", "-f", "mol2", "-o", f"{frag}.frcmod"
    ])
# %%
os.chdir("/pubhome/xtzhang/interaction/FF/mol2")
from rdkit import Chem
for frag, (smarts, smiles, netcharge) in FragSmartsSmiles.items():
    key = f"{frag}_{frag}"
    with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
        lines = f.readlines()
    mapped_line = next((line for line in lines if key in line), None)
    mapped_line = mapped_line.strip()

    mol = Chem.SDMolSupplier(mapped_line, removeHs=False)[0]
    fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    # with Chem.SDWriter(f"{frag}.sdf") as writer:
    #     writer.write(fraga)
    Chem.MolToPDBFile(fraga, f"{frag}.pdb")
    # os.system(f"/pubhome/xtzhang/opt/unicon_1.4.2/unicon -i {frag}.sdf -o {frag}.mol2")
       



# %%
# make param for charmm

os.chdir("/pubhome/xtzhang/interaction/FF/charmm/para")
for name in FragSmartsSmiles.keys():
    os.system(f"/pubhome/soft/silcsbio.2024.1/cgenff/cgenff -i mol2 /pubhome/xtzhang/interaction/FF/mol2/{name}.mol2 -w {name}.mol2 > {name}.str")
# %%
