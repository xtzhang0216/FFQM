#%%
from __future__ import print_function
import argparse
import sys
import os
import datetime
import time

# from omm_readinputs import *
# from omm_readparams import *
# from omm_vfswitch import *

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

from rdkit import Chem
import pandas as pd
import json
from pathlib import Path
import textwrap
import shutil


class fragmol:
    def __init__(self, mol, resname, segid, het):
        self.mol = mol
        self.resname = resname
        self.segid = segid
        self.n_atoms = mol.GetNumAtoms()
        self.het = het
        self.ref_mol = Chem.MolFromPDBFile(f"/pubhome/xtzhang/interaction/FF/drude/ref_mols/{resname}.pdb",removeHs=False)
        with open(f"/pubhome/xtzhang/interaction/FF/drude/ref_mols/{resname}.pdb") as f:
            self.drudename = f.readline().split()[3]
    def assign_crd(self):
        match = self.ref_mol.GetSubstructMatch(self.mol)
        for atom in match:
            pos = self.mol.GetConformer().GetAtomPosition(atom)
            self.ref_mol.GetConformer().SetAtomPosition(atom, pos)
    def write_crd(self, crd_file):
        self.assign_crd()
        atoms = [f"{self.n_atoms:>10}  EXT"]
        for i, atom in enumerate(self.ref_mol.GetAtoms(), start=1):
            pos = self.ref_mol.GetConformer().GetAtomPosition(atom.GetIdx())
            atom_name = atom.GetMonomerInfo().GetName()
            atom_name = atom_name.replace(" ", "")
            line = (
                f"{i:>10}{1:>10}  {self.drudename:<10}{atom_name:<8}"
                f"{pos.x:>20.9f}{pos.y:>20.9f}{pos.z:>20.9f}  "
                f"{self.het}      1               0.0000000000"
            )
            atoms.append(line)
        with open(crd_file, "w") as f:
            f.write("\n".join(atoms))
#%%
def write_frag_inp(resname, crd_file, out_file):
    inp_content = textwrap.dedent(f"""* Written by Xintong Zhang
    DIMENS CHSIZE 3000000 MAXRES 3000000

    stream drude_toppar_2023/toppar_drude_main_protein_2023a.str   ! master, includes protein, water and ions (must be read first)
    stream drude_toppar_2023/toppar_drude_model_2023a.str   ! model compounds

    open read unit 10 card name {crd_file}
    read sequence coor unit 10 RESI

    generate {resname} first none last none setup warn drude dmass 0.4 ! show

    !read coordinates
    rewind unit 10
    read coor card unit 10 append

    coor sdrude
    coor shake

    !SHAKE bonh param nofast -
    !      select  .not. type D*  end -
    !      select  .not. type D*  end

    nbonds atom vatom switch vswitch bycb -
           ctonnb 100.0 ctofnb 120.0 cutnb 16.0 -
           inbfrq -1 imgfrq -1 wmin 1.0 cdie eps 1.0

    energy
    !  minimize drudes
    CONS FIX SELECT .NOT. TYPE D* END
    mini ABNR nstep 2000 nprint 20
    CONS FIX SELE NONE END
    energy

    stop
    """)
    with open(out_file, 'w') as f:
        f.write(inp_content)

def write_complex_inp(fraga_crd, fragb_crd, out_file):
    inp_content = textwrap.dedent(f"""* Written by Xintong Zhang
    DIMENS CHSIZE 3000000 MAXRES 3000000

    stream drude_toppar_2023/toppar_drude_main_protein_2023a.str   ! master, includes protein, water and ions (must be read first)
    stream drude_toppar_2023/toppar_drude_model_2023a.str   ! model compounds

    open read unit 10 card name {fraga_crd}
    read sequence coor unit 10 RESI
    generate HETA first none last none setup warn drude dmass 0.4 ! show
    !read coordinates
    rewind unit 10
    read coor card unit 10 append

    open read unit 10 card name {fragb_crd}
    read sequence coor unit 10 RESI
    generate HETB first none last none setup warn drude dmass 0.4 ! show
    !read coordinates
    rewind unit 10
    read coor card unit 10 append

    coor sdrude
    coor shake

    !SHAKE bonh param nofast -
    !      select  .not. type D*  end -
    !      select  .not. type D*  end

    nbonds atom vatom switch vswitch bycb -
           ctonnb 100.0 ctofnb 120.0 cutnb 16.0 -
           inbfrq -1 imgfrq -1 wmin 1.0 cdie eps 1.0

    energy
    !  minimize drudes
    CONS FIX SELECT .NOT. TYPE D* END
    mini ABNR nstep 2000 nprint 20
    CONS FIX SELE NONE END
    energy

    DEFINE MOL1 SELECT SEGID HETA END
    DEFINE MOL2 SELECT SEGID HETB END
    INTEraction SELE MOL1 END SELE MOL2 END UNIT 6
    stop
    """)
    with open(out_file, 'w') as f:
        f.write(inp_content)



def extract_ene_from_charmmout(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    energy = None
    for line in lines:
        if "ENER>" in line:
            energy = float(line.split()[2])
    return energy
def extract_interene_from_charmmout(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    energy = None
    for line in lines:
        if "INTE>" in line:
            energy = float(line.split()[2])
    return energy


def check_nh_distance(fragA, fragB, fragA_name, idx):
    if fragA_name not in ["ACEM", "MIMM", "MIME", "MIMD"]:
        return True
    for atomA in fragA.GetAtoms():
        if atomA.GetSymbol() == 'N':
            for neighbor in atomA.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    posH = fragA.GetConformer().GetAtomPosition(neighbor.GetIdx())
                    for atomB in fragB.GetAtoms():
                        posB = fragB.GetConformer().GetAtomPosition(atomB.GetIdx())
                        dist = ((posH.x-posB.x)**2 + (posH.y-posB.y)**2 + (posH.z-posB.z)**2)**0.5
                        if dist < 1.5:
                            print(f"[Warning] {fragA_name}-N-H...atom distance < 1.5A at conformer {idx}, skip.")
                            return False
    return True

def main():
    parser = argparse.ArgumentParser(description='CHARMM Drude energy calculation')
    parser.add_argument('--xyz', type=str, default="ACEM_ACEM_00", help='xyz name, e.g. ACEM_HOH_00')
    parser.add_argument('--start', type=int, default=0, help='start index for sdf conformers')
    parser.add_argument('--end', type=int, default=100, help='end index for sdf conformers (exclusive)')
    parser.add_argument('--check_nh', type=bool, default=True, help='check N-H distance for specific molecules')
    args = parser.parse_args()

    xyz = args.xyz
    start = args.start
    end = args.end

    energy_json_path = Path('/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json')
    with open(energy_json_path, 'r') as energy_file:
        qm_energy_data = json.load(energy_file)
    ene = qm_energy_data[f"{xyz}.xyz"]
    results = {}
    with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
        lines = f.readlines()
    sdfpath = next((line for line in lines if xyz in line), None)
    sdfpath = sdfpath.strip()
    sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
    import tempfile
    import time
    WORKDIR = tempfile.mkdtemp(prefix=f"{xyz}_")
    total = len(sdf)
    for idx in range(start, end):
        start_time = time.time()
        if str(idx) not in ene.keys():
            continue
        mol = sdf[idx]
        fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()
        fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if fga_name == "ACEH" or fgb_name == "ACEH":
            continue
        if fga_name == "ACET" and fraga.GetNumAtoms() != 7:
            continue
        if fgb_name == "ACET" and fragb.GetNumAtoms() != 7:
            continue
        if args.check_nh:
            if not check_nh_distance(fraga, fragb, fga_name, idx) or not check_nh_distance(fragb, fraga, fgb_name, idx):
                continue
        conf_dir = os.path.join(WORKDIR, f"conf_{idx}")
        print(f"Processing conformer {idx} / {total} at {conf_dir} at {datetime.datetime.now()}", flush=True)
        os.makedirs(conf_dir, exist_ok=True)
        os.chdir(conf_dir)
        fraga = fragmol(fraga, fga_name, 1, "HETA")
        fragb = fragmol(fragb, fgb_name, 2, "HETB")
        fraga.write_crd(f"{fga_name.lower()}_a.crd")
        fragb.write_crd(f"{fgb_name.lower()}_b.crd")
        write_frag_inp("HETA", f"{fga_name.lower()}_a.crd", "fraga.inp")
        write_frag_inp("HETB", f"{fgb_name.lower()}_b.crd", "fragb.inp")
        write_complex_inp(f"{fga_name.lower()}_a.crd", f"{fgb_name.lower()}_b.crd", "cmx.inp")
        os.system("cp -r /pubhome/xtzhang/interaction/FF/drude/drude_toppar_2023 .")
        os.system("/pubhome/soft/charmm/c49b2/bin/charmm -i cmx.inp > cmx.out")
        os.system("/pubhome/soft/charmm/c49b2/bin/charmm -i fraga.inp > fraga.out")
        os.system("/pubhome/soft/charmm/c49b2/bin/charmm -i fragb.inp > fragb.out")
        charmm_cmx_ene = extract_ene_from_charmmout("cmx.out")
        charmm_fraga_ene_min = extract_ene_from_charmmout("fraga.out")
        charmm_fragb_ene_min = extract_ene_from_charmmout("fragb.out")
        charmm_inter_ene_inter = extract_interene_from_charmmout("cmx.out")
        charmm_inter_ene_min = charmm_cmx_ene - charmm_fraga_ene_min - charmm_fragb_ene_min
        qm_ene = ene[str(idx)]
        results[str(idx)] = [round(charmm_inter_ene_inter, 2), round(charmm_inter_ene_min, 2), round(qm_ene, 2)]
        Chem.MolToPDBFile(mol, f"{idx}.pdb")
        os.chdir("../")
        shutil.rmtree(conf_dir)
    os.chdir("/")
    shutil.rmtree(WORKDIR)
    results_df = pd.DataFrame.from_dict(results, orient='index', columns=['charmm_inter_energy', 'charmm_supermolecule_energy', 'QM_Energy'])
    results_df.index.name = 'Index'
    if fgb_name[0] < fga_name[0]:
        pair = f"{fgb_name}_{fga_name}"
    else:
        pair = f"{fga_name}_{fgb_name}"
    output_file = Path(f'/pubhome/xtzhang/interaction/FF/drude/data/{pair}.csv')
    results_df.to_csv(output_file, mode='a', header=not output_file.exists())

if __name__ == "__main__":
    main()
# ...existing code...
