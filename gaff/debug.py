#%%
import os
import shutil
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import subprocess
import pandas as pd
import json
from pathlib import Path
import sys
from openmm.app import AmberPrmtopFile, AmberInpcrdFile
from openmm import System, Context, Platform
from openmm.app import ForceField
from openmm.unit import *
from openmm.app import NoCutoff,PDBFile
from openmm import VerletIntegrator
#%%

def run_antechamber(file, prefix, netcharge=0):
    subprocess.run([
        "/pubhome/xtzhang/mambaforge/envs/htmd/bin/antechamber", "-i", file, "-fi", "sdf",
        "-o", f"{prefix}.mol2", "-fo", "mol2", "-c", "bcc", "-s", "2", "-nc", f"{netcharge}"
    ])  # generate mol2 with bcc charge
    subprocess.run([
        "/pubhome/xtzhang/mambaforge/envs/htmd/bin/parmchk2", "-i", f"{prefix}.mol2", "-f", "mol2", "-o", f"{prefix}.frcmod"
    ]) # generate mainly dihedral param

def write_leap(prefix):
    with open("leap.in", "w") as f:
        f.write("source leaprc.gaff\n")
        f.write(f"LIG = loadmol2 {prefix}.mol2\n")
        f.write(f"loadamberparams {prefix}.frcmod\n")
        f.write(f"saveamberparm LIG {prefix}.prmtop {prefix}.inpcrd\n")# generate topology and coordinate for specific conformation 
        f.write("quit\n")
        f.write("quit\n")
    subprocess.run(["/pubhome/xtzhang/mambaforge/envs/htmd/bin/tleap", "-f", "leap.in"])

def write_leap_cmx(fga,fgb):
# source leaprc.gaff
# loadamberparams etoh.frcmod
# loadamberparams tolu.frcmod
# ETOH = loadmol2 etoh.mol2
# TOLU = loadmol2 tolu.mol2
# COMPLEX = combine {ETOH TOLU}
# saveamberparm COMPLEX complex.prmtop complex.inpcrd
# quit
    with open("leap.in", "w") as f:
        f.write("source leaprc.gaff\n")
        f.write(f"loadamberparams {fga}_a.frcmod\n")
        f.write(f"loadamberparams {fgb}_b.frcmod\n")
        f.write(f"{fga}_a = loadmol2 {fga}_a.mol2\n")
        f.write(f"{fgb}_b = loadmol2 {fgb}_b.mol2\n")
        f.write(f"COMPLEX = combine {{{fga}_a {fgb}_b}}\n")
        f.write(f"saveamberparm COMPLEX {fga}_{fgb}.prmtop {fga}_{fgb}.inpcrd\n")
        f.write("quit\n")
    subprocess.run(["/pubhome/xtzhang/mambaforge/envs/htmd/bin/tleap", "-f", "leap.in"])
        
def mapping(frag, query,id):
    ref_mol_folder = "/pubhome/xtzhang/interaction/FF/gaff/parameter"

    with open(os.path.join(ref_mol_folder, f"{frag}.mol2"), "r") as f:
        line1 = f.readlines()
    # Modify coordinates in line1 based on uery):
    crd_start=8
    for i in range(query.GetNumAtoms()):
        atom = query.GetAtomWithIdx(i)
        query_symbol = atom.GetSymbol()
        ref_symbol = line1[crd_start + i].split()[1][0]
        if query_symbol != ref_symbol:
            raise ValueError(f"Atom symbol mismatch: {query_symbol} != {ref_symbol}")
        crd = query.GetConformer().GetAtomPosition(i)
        line_idx = crd_start + i
        line1[line_idx] = line1[line_idx][:16] + f"{crd[0]:11.4f}" + line1[line_idx][27:]
        line1[line_idx] = line1[line_idx][:27] + f"{crd[1]:11.4f}" + line1[line_idx][38:]
        line1[line_idx] = line1[line_idx][:38] + f"{crd[2]:11.4f}" + line1[line_idx][49:]
    with open(f"{frag}_{id}.mol2", "w") as f:
        f.writelines(line1)

def calc_ene(top_file, crd_file):
    # Load the topology and coordinates
    psf = AmberPrmtopFile(top_file)
    pdb = AmberInpcrdFile(crd_file)
    system = psf.createSystem(
        nonbondedMethod=NoCutoff,
        constraints=None,
        implicitSolvent=None
    )
    integrator = VerletIntegrator(0.001 * picoseconds)  # 步长设为极小值（不影响单点能量计算）
    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    state = context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    return potential_energy


# def write_min_input():
    with open("min.in", "w") as f:
        f.write("""Single point energy
&cntrl
imin=1,
maxcyc=0,
ntb=0,
cut=999.0,
igb=0,
/
""")

def run_sander(prefix):
    subprocess.run([
        "/pubhome/xtzhang/mambaforge/envs/htmd/bin/sander", "-O", "-i", "min.in", "-o", f"{prefix}.out",
        "-p", f"{prefix}.prmtop", "-c", f"{prefix}.inpcrd",
        "-r", f"{prefix}.rst", "-ref", f"{prefix}.inpcrd"
    ])
    with open(f"{prefix}.out") as f:
        for line in f:
            energy_line = None
            for line in f:
                if "ENERGY" in line:
                    energy_line = next(f, None)
            if energy_line:
                return float(energy_line.split()[1])
    return None


def calc_opls_ene( crd_file):
    # Load the topology and coordinates
    # psf = ForceField("/pubhome/xtzhang/interaction/FF/opls/opls36.zxt.xml")
    psf = ForceField('/pubhome/xtzhang/interaction/FF/opls/param/opls.zxt.xml','tip3p.xml') 

    pdb = PDBFile(crd_file)
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
    return potential_energy



def calc_opls_cmx_ene(top_file1, top_file2, crd_file):

    if top_file1 == top_file2:
        psf = ForceField(top_file1)
    else:
        psf = ForceField(top_file1, top_file2)
    # Load the topology and coordinates
    pdb = PDBFile(crd_file)
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
    return potential_energy

def write_pdb(ref_mol, mol, pdb_file):
    """
    
    """
    for i in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(i).GetSymbol() != ref_mol.GetAtomWithIdx(i).GetSymbol():
            raise ValueError(f"Atom symbol mismatch at index {i}: {mol.GetAtomWithIdx(i).GetSymbol()} != {ref_mol.GetAtomWithIdx(i).GetSymbol()}")
        ref_pos = mol.GetConformer().GetAtomPosition(i)
        ref_mol.GetConformer().SetAtomPosition(i, ref_pos)
    Chem.MolToPDBFile(ref_mol, pdb_file)

def write_cmx_pdb(mol1,mol2,pdb_file):
    length1 = mol1.GetNumAtoms()
    length2 = mol2.GetNumAtoms()
    combined = Chem.CombineMols(mol1, mol2)
    Chem.MolToPDBFile(combined, pdb_file)
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    with open(pdb_file, 'w') as f:
        for i,line in enumerate(lines):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if i < length1:
                    f.write(line[:25] + "1" + line[26:])
                elif length1 <= i < length1 + length2:
                    f.write(line[:25] + "2" + line[26:])
            else:
                f.write(line)


# load input para: --fraga --fragb --xyz
# xyz = sys.argv[1]
xyz="ACET_N1PA_00"

energy_json_path = Path('/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json')
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)
# for key in qm_energy_data.keys():
ene = qm_energy_data[f"{xyz}.xyz"]
results = {}
with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
    lines = f.readlines()
sdfpath = next((line for line in lines if xyz in line), None)
sdfpath = sdfpath.strip()
sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
# Settings
import tempfile
import time
WORKDIR = tempfile.mkdtemp(prefix=f"{xyz}_")
total = len(sdf)
for idx, mol in enumerate(sdf):
    # 
    # start_time = time.time()
    # with open(f"{WORKDIR}/log.txt", "a") as log_file:
    #     log_file.write(f"Processing conformer {idx + 1} / {total}\n")
    if idx > 5000:
        break
    if str(idx) not in ene.keys():
        continue
    if mol is None:
        continue
    fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()

    print(f"Processing conformer {idx + 1} / {total}")

    conf_dir = os.path.join(WORKDIR, f"conf_{idx}")
    os.makedirs(conf_dir, exist_ok=True)
    os.chdir(conf_dir)

    fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if fga_name == "ACEH" or fgb_name == "ACEH":
        continue
    if fga_name == "ACET" and fraga.GetNumAtoms() != 7:
        continue
    if fgb_name == "ACET" and fragb.GetNumAtoms() != 7:
        continue
    # get mol2
    mapping(fga_name, fraga,'a')
    mapping(fgb_name, fragb,'b')
    # get frcmod
    shutil.copy(f"/pubhome/xtzhang/interaction/FF/gaff/parameter/{fga_name}.frcmod", f"{conf_dir}/{fga_name}_a.frcmod")
    shutil.copy(f"/pubhome/xtzhang/interaction/FF/gaff/parameter/{fgb_name}.frcmod", f"{conf_dir}/{fgb_name}_b.frcmod")
    # with Chem.SDWriter(os.path.join(conf_dir, f"{fga_name}.sdf")) as writer:
    #     writer.write(fragb)
    # with Chem.SDWriter(os.path.join(conf_dir, f"{fgb_name}.sdf")) as writer:
    #     writer.write(fraga)
    # with Chem.SDWriter(os.path.join(conf_dir, f"{fga_name}_{fgb_name}.sdf")) as writer:
    #     writer.write(mol)
    # run_antechamber(f"{fga_name}.sdf", f"{fga_name}")
    write_leap(f"{fga_name}_a")
    # run_antechamber(f"{fgb_name}.sdf", f"{fgb_name}")
    write_leap(f"{fgb_name}_b")
    # complex
    # run_antechamber(f"{fga_name}_{fragb_name}.sdf",   f"{fga_name}_{fgb_name}")
    write_leap_cmx(fga_name, fgb_name)
    
    # Energy calculation
    # write_min_input()
    e_a = calc_ene(f"{fga_name}_a.prmtop", f"{fga_name}_a.inpcrd")
    e_b = calc_ene(f"{fgb_name}_b.prmtop", f"{fgb_name}_b.inpcrd")
    e_complex = calc_ene(f"{fga_name}_{fgb_name}.prmtop", f"{fga_name}_{fgb_name}.inpcrd")
    # e_a = run_sander(f"{fga_name}_a")
    # e_b = run_sander(f"{fgb_name}_b")
    # e_complex = run_sander(f"{fga_name}_{fgb_name}")

    interaction_energy = e_complex - (e_a + e_b)

    Chem.MolToPDBFile(fraga, f"{fga_name}_a.pdb")
    Chem.MolToPDBFile(fragb, f"{fgb_name}_b.pdb")
    e_a_opls = calc_opls_ene(f"{fga_name}_a.pdb")
    e_b_opls = calc_opls_ene(f"{fgb_name}_b.pdb")
    write_cmx_pdb(fraga, fragb, f"{fga_name}_{fgb_name}.pdb")
    e_complex_opls = calc_opls_ene(f"{fga_name}_{fgb_name}.pdb")
    interaction_energy_opls = e_complex_opls - (e_a_opls + e_b_opls)
         







    # interaction_energy_kcal = interaction_energy/4.184 
    qm_ene = ene[str(idx)]
    # save .2f
    results[str(idx)] = [round(interaction_energy, 2),
                         round(interaction_energy_opls, 2),
                         round(qm_ene, 2)]

    os.chdir("../")
results_df = pd.DataFrame.from_dict(results, orient='index', columns=['FF_Energy', 'FF_Energy_OPLS', 'QM_Energy'])
results_df.index.name = 'Index'
if fgb_name[0] < fga_name[0]:
    pair = f"{fgb_name}_{fga_name}"
else:
    pair = f"{fga_name}_{fgb_name}"
output_file = Path(f'/pubhome/xtzhang/interaction/FF/gaff/data/{pair}.csv')
results_df.to_csv(output_file, mode='a', header=not output_file.exists())


# %%
rot_sdf_path = next((path for path in Path('/pubhome/lzeng/data/pair25/rot_split/').glob(f'*{xyz}*')), None)
if rot_sdf_path:
    rot_sdf = Chem.SDMolSupplier(str(rot_sdf_path), removeHs=False)
else:
    print(f"No rotated SDF file found for {xyz}")
    rot_sdf = None

rot_ene = {k: v for k, v in ene.items() if k.startswith('rot')}
total = len(sdf)
if rot_ene:
    for idx, mol in enumerate(sdf):
        # 
        # start_time = time.time()
        # with open(f"{WORKDIR}/log.txt", "a") as log_file:
        #     log_file.write(f"Processing conformer {idx + 1} / {total}\n")
        if idx > 5000:
            break
        if f"rot-{idx}" not in rot_ene.keys():
            continue
        if mol is None:
            continue
        fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()

        print(f"Processing conformer {idx + 1} / {total}")

        conf_dir = os.path.join(WORKDIR, f"conf_{idx}")
        os.makedirs(conf_dir, exist_ok=True)
        os.chdir(conf_dir)

        fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if fga_name == "ACEH" or fgb_name == "ACEH":
            continue
        if fga_name == "ACET" and fraga.GetNumAtoms() != 7:
            continue
        if fgb_name == "ACET" and fragb.GetNumAtoms() != 7:
            continue
        # get mol2
        mapping(fga_name, fraga,'a')
        mapping(fgb_name, fragb,'b')
        # get frcmod
        shutil.copy(f"/pubhome/xtzhang/interaction/FF/gaff/parameter/{fga_name}.frcmod", f"{conf_dir}/{fga_name}_a.frcmod")
        shutil.copy(f"/pubhome/xtzhang/interaction/FF/gaff/parameter/{fgb_name}.frcmod", f"{conf_dir}/{fgb_name}_b.frcmod")
        # with Chem.SDWriter(os.path.join(conf_dir, f"{fga_name}.sdf")) as writer:
        #     writer.write(fragb)
        # with Chem.SDWriter(os.path.join(conf_dir, f"{fgb_name}.sdf")) as writer:
        #     writer.write(fraga)
        # with Chem.SDWriter(os.path.join(conf_dir, f"{fga_name}_{fgb_name}.sdf")) as writer:
        #     writer.write(mol)
        # run_antechamber(f"{fga_name}.sdf", f"{fga_name}")
        write_leap(f"{fga_name}_a")
        # run_antechamber(f"{fgb_name}.sdf", f"{fgb_name}")
        write_leap(f"{fgb_name}_b")
        # complex
        # run_antechamber(f"{fga_name}_{fragb_name}.sdf",   f"{fga_name}_{fgb_name}")
        write_leap_cmx(fga_name, fgb_name)
        
        # Energy calculation
        # write_min_input()
        e_a = calc_ene(f"{fga_name}_a.prmtop", f"{fga_name}_a.inpcrd")
        e_b = calc_ene(f"{fgb_name}_b.prmtop", f"{fgb_name}_b.inpcrd")
        e_complex = calc_ene(f"{fga_name}_{fgb_name}.prmtop", f"{fga_name}_{fgb_name}.inpcrd")
        # e_a = run_sander(f"{fga_name}_a")
        # e_b = run_sander(f"{fgb_name}_b")
        # e_complex = run_sander(f"{fga_name}_{fgb_name}")

        interaction_energy = e_complex - (e_a + e_b)

        Chem.MolToPDBFile(fraga, f"{fga_name}_a.pdb")
        Chem.MolToPDBFile(fragb, f"{fgb_name}_b.pdb")
        e_a_opls = calc_opls_ene(f"{fga_name}_a.pdb")
        e_b_opls = calc_opls_ene(f"{fgb_name}_b.pdb")
        write_cmx_pdb(fraga, fragb, f"{fga_name}_{fgb_name}.pdb")
        e_complex_opls = calc_opls_ene(f"{fga_name}_{fgb_name}.pdb")
        interaction_energy_opls = e_complex_opls - (e_a_opls + e_b_opls)
        


        # interaction_energy_kcal = interaction_energy/4.184 
        qm_ene = rot_ene[f"rot-{idx}"]
        # save .2f
        results[f"rot-{idx}"] = [round(interaction_energy, 2),
                                 round(interaction_energy_opls, 2),
                                 round(qm_ene, 2)]

        os.chdir("../")
results_df = pd.DataFrame.from_dict(results, orient='index', columns=['FF_Energy', 'FF_Energy_OPLS', 'QM_Energy'])
results_df.index.name = 'Index'
if fgb_name[0] < fga_name[0]:
    pair = f"{fgb_name}_{fga_name}"
else:
    pair = f"{fga_name}_{fgb_name}"
output_file = Path(f'/pubhome/xtzhang/interaction/FF/gaff/data/{pair}.csv')
results_df.to_csv(output_file, mode='a', header=not output_file.exists())
