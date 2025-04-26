from rdkit import Chem
import json
from pathlib import Path
energy_json_path = Path('/pubhome/xtzhang/interaction/data/pdbpairs/total_inteng_comb.json')
xyz = "ETAM_ETOH"
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)
ene = qm_energy_data[f"{xyz}.xyz"]
with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
    lines = f.readlines()
sdfpath = next((line for line in lines if xyz in line), None)
sdfpath = sdfpath.strip()
sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)

for idx, mol in enumerate(sdf):

    # start_time = time.time()
    # with open(f"{WORKDIR}/log.txt", "a") as log_file:
    #     log_file.write(f"Processing conformer {idx + 1} / {total}\n")
    if idx != 72:
        continue
    if str(idx) not in ene.keys():
        continue
    if mol is None:
        continue
    fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()
    fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    Chem.MolToPDBFile(fraga, f"/pubhome/xtzhang/interaction/FF/amoeba/test/{fga_name}_a.pdb")
    Chem.MolToPDBFile(fragb, f"/pubhome/xtzhang/interaction/FF/amoeba/test/{fgb_name}_b.pdb")
    Chem.MolToPDBFile(mol, f"/pubhome/xtzhang/interaction/FF/amoeba/test/{fga_name}_{fgb_name}.pdb")